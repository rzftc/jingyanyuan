% solve_hourly_dispatch_ptdf.m
function [u_ac, u_ev, cost, exitflag] = solve_hourly_dispatch_ptdf( ...
    num_ac, num_ev, ...
    p_ac_t, p_ev_t, ...
    c_ac_t, c_ev_t, P_req_t, ...
    n_ac_max_from_ga, n_ev_max_from_ga, ...
    Location_AC, Location_EV, PTDF_matrix, ...
    P_Line_Base_t, P_Line_Max, N_bus, N_line)
% SOLVE_HOURLY_DISPATCH_PTDF (下层求解器)
% 使用 MILP (intlinprog) 求解单时间步的最小成本调度问题，
% 同时满足功率需求、GA数量约束和PTDF网络潮流约束。

    % 决策变量 x = [u_ac_1, ..., u_ac_N; u_ev_1, ..., u_ev_M]
    n_vars_ac = num_ac;
    n_vars_ev = num_ev;
    n_vars = n_vars_ac + n_vars_ev;

    if n_vars == 0 || P_req_t <= 1e-6
        u_ac = zeros(num_ac, 1);
        u_ev = zeros(num_ev, 1);
        cost = 0;
        exitflag = 1;
        return;
    end

    % --- 1. 目标函数 (最小化成本) ---
    % f = [p_ac_i * c_ac_i; p_ev_j * c_ev_j]
    f_obj_ac = p_ac_t .* c_ac_t;
    f_obj_ev = p_ev_t .* c_ev_t;
    f_obj = [f_obj_ac; f_obj_ev];

    % 决策变量类型 (全部为整数)
    intcon = 1:n_vars;
    
    % 边界 (0 <= u <= 1)
    lb = zeros(n_vars, 1);
    ub = ones(n_vars, 1);
    
    A_ineq = [];
    b_ineq = [];

    % --- 2. 约束 ---
    
    % --- 约束 1: 功率平衡约束 ---
    % sum(u_ac .* p_ac) + sum(u_ev .* p_ev) >= P_req_t
    % 转换为: -sum(...) <= -P_req_t
    A_ineq_power = -[p_ac_t', p_ev_t'];
    b_ineq_power = -P_req_t;
    A_ineq = [A_ineq; A_ineq_power];
    b_ineq = [b_ineq; b_ineq_power];

    % --- 约束 2: 上层GA的数量约束 ---
    % sum(u_ac) <= n_ac_max_from_ga
    % sum(u_ev) <= n_ev_max_from_ga
    A_ineq_count_ac = [ones(1, n_vars_ac), zeros(1, n_vars_ev)];
    b_ineq_count_ac = n_ac_max_from_ga;
    
    A_ineq_count_ev = [zeros(1, n_vars_ac), ones(1, n_vars_ev)];
    b_ineq_count_ev = n_ev_max_from_ga;
    
    A_ineq = [A_ineq; A_ineq_count_ac; A_ineq_count_ev];
    b_ineq = [b_ineq; b_ineq_count_ac; b_ineq_count_ev];

    % --- 约束 3: PTDF 潮流约束 ---
    if N_line > 0 && N_bus > 0
        % a. 计算节点注入功率矩阵 (N_bus x N_vars)
        % P_inj_k_vec(k, i) = 当设备 i (AC) 启动时，节点 k 注入的功率
        P_inj_k_vec = zeros(N_bus, n_vars);
        for i = 1:n_vars_ac
            node_idx = Location_AC(i);
            P_inj_k_vec(node_idx, i) = p_ac_t(i);
        end
        for j = 1:n_vars_ev
            node_idx = Location_EV(j);
            P_inj_k_vec(node_idx, n_vars_ac + j) = p_ev_t(j);
        end
        
        % b. 计算潮流变化矩阵 (N_line x N_vars)
        % Delta_P_line_vec = PTDF * P_inj_k_vec
        % 每一列代表启动该设备(u=1)时，对所有线路的潮流贡献
        Delta_P_line_vec = PTDF_matrix * P_inj_k_vec;
        
        % c. 设置约束
        % P_Line_Base_t + Delta_P_line_vec * x <= P_Line_Max
        % -P_Line_Base_t - Delta_P_line_vec * x <= P_Line_Max
        
        % 约束 (1): Delta_P_line_vec * x <= P_Line_Max - P_Line_Base_t
        A_ineq_ptdf_upper = Delta_P_line_vec'; % (N_vars x N_line)
        b_ineq_ptdf_upper = P_Line_Max - P_Line_Base_t;
        
        % 约束 (2): -Delta_P_line_vec * x <= P_Line_Max + P_Line_Base_t
        A_ineq_ptdf_lower = -Delta_P_line_vec'; % (N_vars x N_line)
        b_ineq_ptdf_lower = P_Line_Max + P_Line_Base_t;

        % A_ineq 在 intlinprog 中是 A*x <= b, 
        % 我们需要转置 A 矩阵 (N_line x N_vars) -> (N_vars x N_line)
        % 不对，intlinprog 的 A 是 (n_constraints x n_vars)
        
        A_ineq_ptdf_upper = Delta_P_line_vec; % (N_line x N_vars)
        b_ineq_ptdf_upper = P_Line_Max - P_Line_Base_t;
        
        A_ineq_ptdf_lower = -Delta_P_line_vec; % (N_line x N_vars)
        b_ineq_ptdf_lower = P_Line_Max + P_Line_Base_t;
        
        A_ineq = [A_ineq; A_ineq_ptdf_upper; A_ineq_ptdf_lower];
        b_ineq = [b_ineq; b_ineq_ptdf_upper; b_ineq_ptdf_lower];
    end
    
    % --- 3. 求解 ---
   options = optimoptions('intlinprog', 'Display', 'off', 'MaxTime', 60);
    
    [x_sol, cost, exitflag] = intlinprog(f_obj, intcon, A_ineq, b_ineq, [], [], lb, ub, [], options);

    if exitflag > 0
        % (成功)
        u_ac = x_sol(1:n_vars_ac);
        u_ev = x_sol(n_vars_ac+1:end);
    else
        % (求解失败或超时)
        if exitflag == 0
            warning('intlinprog在时段 t=%d (非全局) 达到最大时间(MaxTime)限制，未找到最优解。', t_global_for_warning); % (t_global需要从主循环传入)
        end
        u_ac = zeros(num_ac, 1);
        u_ev = zeros(num_ev, 1);
        cost = NaN; % 标记失败
    end
end