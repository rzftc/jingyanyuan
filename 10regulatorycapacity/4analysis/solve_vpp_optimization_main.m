% solve_vpp_optimization_main.m
function [U_ac_optimal, U_ev_optimal, optimal_total_cost, exitflag, output] = solve_vpp_optimization_main( ...
    ac_params, ev_params, grid_params, ga_options)
% SOLVE_VPP_OPTIMIZATION_MAIN - 虚拟电厂全时域优化调度主函数
%
% 目标: 最小化总调度成本
% 约束: 
%   1. 每时刻满足电网调节需求
%   2. 全局SDCI指标改善
%   3. 全局Spearman Rho指标改善
%   4. 网络潮流不越限 (PTDF约束)
%
% 输入:
%   ac_params: 空调设备参数结构体
%     .num: (int) 空调数量
%     .P_potential: (num_ac x T) 每台空调在每个时间步的潜在调节功率
%     .Cost: (num_ac x 1) 每台空调的单位调节成本 (c_ACi)
%     .Location: (num_ac x 1) 每台空调所在的电网节点编号
%
%   ev_params: 电动汽车设备参数结构体
%     .num: (int) EV数量
%     .P_potential: (num_ev x T) 每台EV在每个时间步的潜在调节功率
%     .Cost: (num_ev x 1) 每台EV的单位调节成本 (c_EVj)
%     .Location: (num_ev x 1) 每台EV所在的电网节点编号
%
%   grid_params: 电网与需求参数结构体
%     .T: (int) 总时域长度 (时间步数)
%     .P_req: (T x 1) 每个时间步的电网调节需求 (Preq(t))
%     .SDCI_orig: (scalar) 原始或基准SDCI指标
%     .rho_orig: (scalar) 原始或基准Spearman Rho指标
%     .PTDF: (N_line x N_bus) 功率传输分布因子矩阵
%     .P_line_base: (N_line x T) 线路基础潮流
%     .P_line_max: (N_line x 1) 线路潮流上限
%     .N_bus: (int) 电网节点/母线数量
%     .N_line: (int) 电网线路数量
%
%   ga_options: (object) 遗传算法选项, 由 optimoptions('ga',...) 创建
%
% 输出:
%   U_ac_optimal: (num_ac x T) 空调的最优启停决策矩阵 (u_ACi(t))
%   U_ev_optimal: (num_ev x T) EV的最优启停决策矩阵 (u_EVj(t))
%   optimal_total_cost: (scalar) 优化后的最低总成本
%   exitflag, output: GA求解器的返回状态

    %% 1. 参数解构与准备
    num_ac = ac_params.num;
    num_ev = ev_params.num;
    T_horizon = grid_params.T;

    % 决策变量总数
    num_vars_ac = num_ac * T_horizon;
    num_vars_ev = num_ev * T_horizon;
    N_vars_total = num_vars_ac + num_vars_ev;
    
    if N_vars_total == 0
        warning('没有可调度的设备，跳过优化。');
        U_ac_optimal = zeros(num_ac, T_horizon);
        U_ev_optimal = zeros(num_ev, T_horizon);
        optimal_total_cost = 0;
        exitflag = -1;
        output = struct('message', '无设备可供优化');
        return;
    end

    %% 2. 定义适应度函数 (Fitness Function)
    % 适应度函数内部会计算成本和所有约束的惩罚项
    fitness_function = @(x) objective_function_with_penalties(x, ac_params, ev_params, grid_params);

    %% 3. 运行遗传算法
    fprintf('开始执行遗传算法优化，总决策变量数: %d...\n', N_vars_total);
    
    % GA 输入参数:
    % @(x) fitness_function(x) - 适应度函数句柄
    % N_vars_total              - 变量数量
    % [], [], [], []            - 线性不等式和等式约束 (A, b, Aeq, beq)，此处不用
    % zeros(1, N_vars_total)    - 变量下界 (0)
    % ones(1, N_vars_total)     - 变量上界 (1)
    % []                        - 非线性约束函数 (已在适应度函数中用罚函数法实现)
    % 1:N_vars_total            - 整数变量索引
    % ga_options                - GA配置
    
    [x_sol, fval, exitflag, output] = ga(fitness_function, ...
                                       N_vars_total, ...
                                       [], [], [], [], ...
                                       zeros(1, N_vars_total), ones(1, N_vars_total), ...
                                       [], 1:N_vars_total, ga_options);

    fprintf('遗传算法执行完毕，退出标志: %d。\n', exitflag);

    %% 4. 结果解析
    optimal_total_cost = fval; % 适应度函数返回的就是包含惩罚的成本

    % 将一维的决策变量向量 x_sol 转换回二维的决策矩阵 U
    if ~isempty(x_sol)
        % 解析空调决策
        if num_ac > 0
            u_ac_flat = x_sol(1:num_vars_ac);
            U_ac_optimal = reshape(logical(u_ac_flat), num_ac, T_horizon);
        else
            U_ac_optimal = zeros(0, T_horizon);
        end
        % 解析EV决策
        if num_ev > 0
            u_ev_flat = x_sol(num_vars_ac + 1 : end);
            U_ev_optimal = reshape(logical(u_ev_flat), num_ev, T_horizon);
        else
            U_ev_optimal = zeros(0, T_horizon);
        end
    else
        warning('GA 未返回有效解。');
        U_ac_optimal = zeros(num_ac, T_horizon);
        U_ev_optimal = zeros(num_ev, T_horizon);
    end
end

%% 嵌套函数: 适应度函数 (包含成本计算和罚函数)
function fitness = objective_function_with_penalties(x, ac_params, ev_params, grid_params)
    
    % --- 解码决策变量 ---
    T = grid_params.T;
    num_ac = ac_params.num;
    num_ev = ev_params.num;
    num_vars_ac = num_ac * T;

    U_ac = reshape(logical(x(1:num_vars_ac)), num_ac, T);
    U_ev = reshape(logical(x(num_vars_ac+1:end)), num_ev, T);

    [cite_start]% --- 1. 核心目标：计算总成本 [cite: 18] ---
    cost_ac = sum(U_ac .* ac_params.P_potential .* ac_params.Cost, 'all');
    cost_ev = sum(U_ev .* ev_params.P_potential .* ev_params.Cost, 'all');
    total_cost = cost_ac + cost_ev;

    % --- 2. 计算约束违反并施加惩罚 ---
    penalty = 0;
    
    [cite_start]% --- 约束1: 每时刻调节能力约束 [cite: 25] ---
    P_ac_dispatch = sum(U_ac .* ac_params.P_potential, 1);
    P_ev_dispatch = sum(U_ev .* ev_params.P_potential, 1);
    P_total_dispatch = P_ac_dispatch + P_ev_dispatch;
    
    power_shortage = max(0, grid_params.P_req' - P_total_dispatch);
    % 惩罚项与缺额的平方成正比
    penalty = penalty + 1e6 * sum(power_shortage.^2);

    [cite_start]% --- 约束2 & 3: SDCI 和 Rho 指标约束 [cite: 28, 31] ---
    n_ac_t = sum(U_ac, 1)';
    n_ev_t = sum(U_ev, 1)';
    avg_P_ac_t = (P_ac_dispatch ./ (n_ac_t' + eps))';
    avg_P_ev_t = (P_ev_dispatch ./ (n_ev_t' + eps))';
    
    sdci_opt = calculateSDCI(n_ac_t, n_ev_t, avg_P_ac_t, avg_P_ev_t);
    rho_opt = calculateSpearmanRho(n_ac_t, avg_P_ac_t, n_ev_t, avg_P_ev_t);

    sdci_violation = max(0, sdci_opt - grid_params.SDCI_orig);
    rho_violation = max(0, abs(rho_opt) - abs(grid_params.rho_orig)); % 约束rho的绝对值
    
    penalty = penalty + 1e5 * (sdci_violation^2 + rho_violation^2);

    [cite_start]% --- 约束4: 网络潮流约束 (PTDF) [cite: 34-45] ---
    if grid_params.N_line > 0
        % a. [cite_start]计算各节点净注入功率 [cite: 37]
        Delta_Pinj = zeros(grid_params.N_bus, T);
        for k = 1:grid_params.N_bus
            ac_at_bus_k = (ac_params.Location == k);
            ev_at_bus_k = (ev_params.Location == k);
            
            if any(ac_at_bus_k)
                Delta_Pinj(k,:) = Delta_Pinj(k,:) + sum(U_ac(ac_at_bus_k,:) .* ac_params.P_potential(ac_at_bus_k,:), 1);
            end
            if any(ev_at_bus_k)
                Delta_Pinj(k,:) = Delta_Pinj(k,:) + sum(U_ev(ev_at_bus_k,:) .* ev_params.P_potential(ev_at_bus_k,:), 1);
            end
        end
        
        % b. [cite_start]计算线路潮流变化 [cite: 41]
        Delta_P_line = grid_params.PTDF * Delta_Pinj;
        
        % c. 计算最终线路潮流
        P_line_final = grid_params.P_line_base + Delta_P_line;
        
        % d. [cite_start]计算越限惩罚 [cite: 45]
        line_overload = max(0, abs(P_line_final) - grid_params.P_line_max);
        penalty = penalty + 1e7 * sum(line_overload.^2, 'all');
    end

    % --- 最终适应度 ---
    fitness = total_cost + penalty;
end