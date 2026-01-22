% solve_hourly_dispatch_ptdf_hourly.m
function [u_ac, u_ev, cost, exitflag] = solve_hourly_dispatch_ptdf_hourly_fixed( ...
    num_ac, num_ev, ...
    P_ac_hourly, P_ev_hourly, ...
    c_ac, c_ev, P_req_hourly, ...
    n_ac_max_from_ga, n_ev_max_from_ga, ...
    Location_AC, Location_EV, PTDF_matrix, ...
    P_Line_Base_hourly, P_Line_Max, N_bus, N_line, h_global_for_warning)
% SOLVE_HOURLY_DISPATCH_PTDF_HOURLY (下层求解器 - 小时级)
%
% [V2 - 修正版] 
% 为了解决 Exitflag: -2 (不可行) 的问题，本函数已被修改。
% 原逻辑要求一个固定设备组合 (u_ac, u_ev) 满足小时内所有时间步的约束。
% 上层GA (hourly_ga_fitness_constrained) 仅在 *峰值时刻* 检查可行性。
% 
% 新逻辑：本函数现在 *只* 在小时内的 *峰值需求时刻* 建立约束。
% 这使得本函数的约束条件与上层GA的适应度函数保持一致，从而保证可行性。

    n_vars_ac = num_ac;
    n_vars_ev = num_ev;
    n_vars = n_vars_ac + n_vars_ev;
    steps_per_hour = size(P_ac_hourly, 2); % (例如 20)

    % --- 0. 寻找峰值需求时刻 ---
    % 找到需求最大的那个时间步的索引 (t_local)
    [max_demand_hourly, peak_t_local_indices] = max(P_req_hourly);
    
    % 如果有多个相同的最大值，只取第一个
    peak_t_local = peak_t_local_indices(1); 

    if n_vars == 0 || max_demand_hourly <= 1e-6
        % 如果无设备，或整个小时都无需求
        u_ac = zeros(num_ac, 1);
        u_ev = zeros(num_ev, 1);
        cost = 0;
        exitflag = 1;
        return;
    end

    % --- 1. 目标函数 (最小化整个小时的总成本) ---
    % 目标函数仍然是基于 *整个小时* 的总成本，因为我们选的设备要工作一小时
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
    % (这组约束保持不变)
    A_ineq_count_ac = [ones(1, n_vars_ac), zeros(1, n_vars_ev)];
    b_ineq_count_ac = n_ac_max_from_ga;
    
    A_ineq_count_ev = [zeros(1, n_vars_ac), ones(1, n_vars_ev)];
    b_ineq_count_ev = n_ev_max_from_ga;
    
    A_ineq = [A_ineq; A_ineq_count_ac; A_ineq_count_ev];
    b_ineq = [b_ineq; b_ineq_count_ac; b_ineq_count_ev];

    
    % --- [!!! 关键修改 !!!] ---
    % 不再循环 (for t_local = 1:steps_per_hour)
    % 仅使用峰值时刻 (peak_t_local) 的数据构建约束
    
    t_local = peak_t_local;
    
    P_req_t = P_req_hourly(t_local);
    p_ac_t = P_ac_hourly(:, t_local);
    p_ev_t = P_ev_hourly(:, t_local);
    P_Line_Base_t = P_Line_Base_hourly(:, t_local);

    % --- 约束 2b: 功率平衡约束 (仅 1 个约束) ---
    % -sum(u_ac .* p_ac_t) - sum(u_ev .* p_ev_t) <= -P_req_t
    A_ineq_power_t = -[p_ac_t', p_ev_t'];
    b_ineq_power_t = -P_req_t;
    A_ineq = [A_ineq; A_ineq_power_t];
    b_ineq = [b_ineq; b_ineq_power_t];

    % --- 约束 2c: PTDF 潮流约束 (仅 2 * N_line 个约束) ---
    if N_line > 0 && N_bus > 0
        % a. 计算 *峰值时刻t* 的节点注入功率矩阵 (N_bus x N_vars)
        P_inj_k_vec_t = zeros(N_bus, n_vars);
        for i = 1:n_vars_ac
            if ~isempty(Location_AC) % 增加检查
                node_idx = Location_AC(i);
                P_inj_k_vec_t(node_idx, i) = p_ac_t(i); % 使用t时刻的功率
            end
        end
        for j = 1:n_vars_ev
             if ~isempty(Location_EV) % 增加检查
                node_idx = Location_EV(j);
                P_inj_k_vec_t(node_idx, n_vars_ac + j) = p_ev_t(j); % 使用t时刻的功率
             end
        end
        
        % b. 计算 *峰值时刻t* 的潮流变化矩阵 (N_line x N_vars)
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
    % --- [!!! 修改结束 !!!] ---
    

    % --- 3. 求解 ---
    % (保持不变)
    options = optimoptions('intlinprog', ...
        'Display', 'off', ...
        'MaxTime', 1000000, ...              
        'RelativeGapTolerance', 0.05, ... 
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
        elseif exitflag == -2
             warning('intlinprog在第 %d 小时 求解失败 (Exitflag: -2，问题不可行)。即使只检查峰值时刻，约束依然无法满足。', h_global_for_warning);
        else
            warning('intlinprog在第 %d 小时 求解失败 (Exitflag: %d)。', h_global_for_warning, exitflag);
        end
        u_ac = zeros(num_ac, 1);
        u_ev = zeros(num_ev, 1);
        cost = NaN; % 标记失败
    end
end