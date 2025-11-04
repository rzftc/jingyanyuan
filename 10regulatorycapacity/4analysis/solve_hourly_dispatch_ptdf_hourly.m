% solve_hourly_dispatch_ptdf_hourly.m
function [u_ac, u_ev, cost, exitflag] = solve_hourly_dispatch_ptdf_hourly( ...
    num_ac, num_ev, ...
    P_ac_hourly, P_ev_hourly, ...
    c_ac, c_ev, P_req_hourly, ...
    n_ac_max_from_ga, n_ev_max_from_ga, ...
    Location_AC, Location_EV, PTDF_matrix, ...
    P_Line_Base_hourly, P_Line_Max, N_bus, N_line, h_global_for_warning)
% SOLVE_HOURLY_DISPATCH_PTDF_HOURLY (下层求解器 - 小时级)
%
% 使用 MILP (intlinprog) 求解 *整个小时* 的最小成本调度问题。
% 决策变量是选择哪些设备参与本小时（在所有时间步中保持不变），
% 以满足本小时内 *每个时间步* 的功率需求和网络约束。
%
% 输入:
%   P_ac_hourly (num_ac x steps_per_hour) - AC设备在该小时内每个时间步的潜力
%   P_ev_hourly (num_ev x steps_per_hour) - EV设备在该小时内每个时间步的潜力
%   P_req_hourly (steps_per_hour x 1)     - 该小时内每个时间步的需求
%   P_Line_Base_hourly (N_line x steps_per_hour) - 该小时内每个时间步的线路基础潮流
%   (其他参数为标量或向量)

    n_vars_ac = num_ac;
    n_vars_ev = num_ev;
    n_vars = n_vars_ac + n_vars_ev;
    steps_per_hour = size(P_ac_hourly, 2);

    if n_vars == 0 || all(P_req_hourly <= 1e-6)
        u_ac = zeros(num_ac, 1);
        u_ev = zeros(num_ev, 1);
        cost = 0;
        exitflag = 1;
        return;
    end

    % --- 1. 目标函数 (最小化整个小时的总成本) ---
    % f = [sum(P_ac_i(t)*c_ac_i) for t in hour; ...]
    f_obj_ac = sum(P_ac_hourly .* repmat(c_ac, 1, steps_per_hour), 2);
    f_obj_ev = sum(P_ev_hourly .* repmat(c_ev, 1, steps_per_hour), 2);
    f_obj = [f_obj_ac; f_obj_ev];

    intcon = 1:n_vars;
    lb = zeros(n_vars, 1);
    ub = ones(n_vars, 1);
    
    % --- 2. 约束 (A_ineq * x <= b_ineq) ---
    A_ineq = [];
    b_ineq = [];

    % --- 约束 2a: 上层GA的数量约束 (2个约束) ---
    % sum(u_ac) <= n_ac_max_from_ga
    % sum(u_ev) <= n_ev_max_from_ga
    A_ineq_count_ac = [ones(1, n_vars_ac), zeros(1, n_vars_ev)];
    b_ineq_count_ac = n_ac_max_from_ga;
    
    A_ineq_count_ev = [zeros(1, n_vars_ac), ones(1, n_vars_ev)];
    b_ineq_count_ev = n_ev_max_from_ga;
    
    A_ineq = [A_ineq; A_ineq_count_ac; A_ineq_count_ev];
    b_ineq = [b_ineq; b_ineq_count_ac; b_ineq_count_ev];

    
    % --- 循环构建每个时间步的约束 ---
    for t_local = 1:steps_per_hour
        
        P_req_t = P_req_hourly(t_local);
        if P_req_t <= 1e-6
            continue; % 如果当前时间步无需求，则跳过功率和网络约束
        end

        p_ac_t = P_ac_hourly(:, t_local);
        p_ev_t = P_ev_hourly(:, t_local);
        P_Line_Base_t = P_Line_Base_hourly(:, t_local);

        % --- 约束 2b: 功率平衡约束 (steps_per_hour 个约束) ---
        % -sum(u_ac .* p_ac_t) - sum(u_ev .* p_ev_t) <= -P_req_t
        A_ineq_power_t = -[p_ac_t', p_ev_t'];
        b_ineq_power_t = -P_req_t;
        A_ineq = [A_ineq; A_ineq_power_t];
        b_ineq = [b_ineq; b_ineq_power_t];

        % --- 约束 2c: PTDF 潮流约束 (2 * N_line * steps_per_hour 个约束) ---
        if N_line > 0 && N_bus > 0
            % a. 计算 *当前时间步t* 的节点注入功率矩阵 (N_bus x N_vars)
            P_inj_k_vec_t = zeros(N_bus, n_vars);
            for i = 1:n_vars_ac
                node_idx = Location_AC(i);
                P_inj_k_vec_t(node_idx, i) = p_ac_t(i); % 使用t时刻的功率
            end
            for j = 1:n_vars_ev
                node_idx = Location_EV(j);
                P_inj_k_vec_t(node_idx, n_vars_ac + j) = p_ev_t(j); % 使用t时刻的功率
            end
            
            % b. 计算 *当前时间步t* 的潮流变化矩阵 (N_line x N_vars)
            Delta_P_line_vec_t = PTDF_matrix * P_inj_k_vec_t;
            
            % c. 设置约束
            % 约束 (1): Delta_P_line_vec_t * x <= P_Line_Max - P_Line_Base_t
            A_ineq_ptdf_upper_t = Delta_P_line_vec_t; % (N_line x N_vars)
            b_ineq_ptdf_upper_t = P_Line_Max - P_Line_Base_t;
            
            % 约束 (2): -Delta_P_line_vec_t * x <= P_Line_Max + P_Line_Base_t
            A_ineq_ptdf_lower_t = -Delta_P_line_vec_t; % (N_line x N_vars)
            b_ineq_ptdf_lower_t = P_Line_Max + P_Line_Base_t;
            
            A_ineq = [A_ineq; A_ineq_ptdf_upper_t; A_ineq_ptdf_lower_t];
            b_ineq = [b_ineq; b_ineq_ptdf_upper_t; b_ineq_ptdf_lower_t];
        end
    end
    
    % --- 3. 求解 ---
    % 使用与 solve_hourly_dispatch_ptdf 相同的宽松设置
    options = optimoptions('intlinprog', ...
        'Display', 'off', ...
        'MaxTime', 120, ...              
        'RelativeGapTolerance', 0.05, ... % 放宽到 5%
        'AbsoluteGapTolerance', 10);     
    
    [x_sol, cost, exitflag] = intlinprog(f_obj, intcon, A_ineq, b_ineq, [], [], lb, ub, [], options);

    if exitflag > 0
        % (成功)
        u_ac = x_sol(1:n_vars_ac);
        u_ev = x_sol(n_vars_ac+1:end);
    else
        % (求解失败或超时)
        if exitflag == 0
            warning('intlinprog在第 %d 小时 达到最大时间(MaxTime)限制，未找到小时最优解。', h_global_for_warning);
        else
            warning('intlinprog在第 %d 小时 求解失败 (Exitflag: %d)。', h_global_for_warning, exitflag);
        end
        u_ac = zeros(num_ac, 1);
        u_ev = zeros(num_ev, 1);
        cost = NaN; % 标记失败
    end
end