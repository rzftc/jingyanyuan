%% *** 新增：局部函数: solve_hourly_dispatch_ptdf ***
function [u_ac_optimal, u_ev_optimal, optimal_cost, exitflag] = solve_hourly_dispatch_ptdf( ...
    num_ac, num_ev, p_ac_t, p_ev_t, c_ac, c_ev, P_req_t, ... % 基础参数
    n_ac_hourly_max, n_ev_hourly_max, ...                     % 上层GA给定的数量约束
    Location_AC, Location_EV, PTDF_matrix, ...                % 网络拓扑参数
    P_Line_Base_t, P_Line_Max, N_bus, N_line)                 % 网络潮流参数
    
    % p_ac_t, p_ev_t, c_ac, c_ev, Location_AC, Location_EV 必须是列向量
    % P_Line_Base_t, P_Line_Max 必须是列向量
    
    n_vars_ac = num_ac;
    n_vars_ev = num_ev;
    n_vars = n_vars_ac + n_vars_ev;
    
    u_ac_optimal = zeros(num_ac, 1); 
    u_ev_optimal = zeros(num_ev, 1);
    optimal_cost = 0;
    
    if n_vars == 0
        exitflag = -200; % 无设备
        return;
    end
    
    if P_req_t <= 1e-6 % 如果需求为0或负
        exitflag = 2; % 无需优化，需求已满足
        return;
    end

    % 1. 目标函数 (最小化成本) 
    f_obj_ac = p_ac_t .* c_ac; % (num_ac x 1)
    f_obj_ev = p_ev_t .* c_ev; % (num_ev x 1)
    f_obj = [f_obj_ac; f_obj_ev];

    % 2. 约束 (A_ineq * x <= b_ineq)
    A_ineq = [];
    b_ineq = [];

    % 约束 1: 功率需求 
    % sum(u_ac .* p_ac_t) + sum(u_ev .* p_ev_t) >= P_req_t
    % => -sum(...) <= -P_req_t
    A_ineq(1, :) = [-p_ac_t', -p_ev_t'];
    b_ineq(1, 1) = -P_req_t;

    % 约束 2: 设备数量约束 (来自上层GA)
    if n_vars_ac > 0
        A_ineq(end+1, 1:n_vars_ac) = 1;
        b_ineq(end+1, 1) = n_ac_hourly_max;
    end
    if n_vars_ev > 0
        A_ineq(end+1, n_vars_ac+1:n_vars) = 1;
        b_ineq(end+1, 1) = n_ev_hourly_max;
    end

    % 约束 3: 潮流约束 (PTDF) [cite: 34, 45]
    if N_line > 0 && N_bus > 0
        % 3.1 构造 A_ptdf 矩阵 (N_line x n_vars)
        % A_ptdf(l, i) = 线路l因设备i启动而产生的潮流变化
        A_ptdf = zeros(N_line, n_vars);
        
        % 填充AC部分 [cite: 37]
        if n_vars_ac > 0
            for i = 1:n_vars_ac
                 node_k = Location_AC(i);
                 if node_k > 0 && node_k <= N_bus
                    % 功率注入 Delta_P_inj = u_i * p_i [cite: 37]
                    % 潮流变化 Delta_PL = PTDF * Delta_P_inj [cite: 41]
                    % A_ineq*x 中的系数即为 PTDF(l, k) * p_i(t)
                    A_ptdf(:, i) = PTDF_matrix(:, node_k) * p_ac_t(i);
                 end
            end
        end
        % 填充EV部分 [cite: 37]
        if n_vars_ev > 0
            for j = 1:n_vars_ev
                 node_k = Location_EV(j);
                 if node_k > 0 && node_k <= N_bus
                    A_ptdf(:, n_vars_ac + j) = PTDF_matrix(:, node_k) * p_ev_t(j);
                 end
            end
        end

        % 3.2 添加约束 
        % P_Line_Base_t + A_ptdf * x <= P_Line_Max
        % => A_ptdf * x <= P_Line_Max - P_Line_Base_t
        A_ineq = [A_ineq; A_ptdf];
        b_ineq = [b_ineq; (P_Line_Max - P_Line_Base_t)];

        % -P_Line_Max <= P_Line_Base_t + A_ptdf * x
        % => -A_ptdf * x <= P_Line_Max + P_Line_Base_t
        A_ineq = [A_ineq; -A_ptdf];
        b_ineq = [b_ineq; (P_Line_Max + P_Line_Base_t)];
    end

    % 3. 边界和整数定义
    lb_vars = zeros(n_vars, 1);
    ub_vars = ones(n_vars, 1);
    intcon_vars = 1:n_vars;
    A_eq = [];
    b_eq = [];

    % 4. 求解
    options = optimoptions('intlinprog', 'Display', 'off');
    [x_sol, fval_sol, exitflag] = intlinprog(f_obj, intcon_vars, A_ineq, b_ineq, A_eq, b_eq, lb_vars, ub_vars, [], options);

    % 5. 提取结果
    if exitflag > 0
        if n_vars_ac > 0
            u_ac_optimal = x_sol(1:n_vars_ac);
        end
        if n_vars_ev > 0
            u_ev_optimal = x_sol(n_vars_ac+1:end);
        end
        optimal_cost = fval_sol;
    else
        % 求解失败
        optimal_cost = NaN;
    end
end
