% nonlinear_constraints_ga_ptdf.m
function [c_ineq, c_eq] = nonlinear_constraints_ga_ptdf(x, num_ac, num_ev, ...
                                P_ac_ind_t, P_ev_ind_t, P_req_t, T, ...
                                SDCI_orig, rho_orig, ...
                                Location_AC, Location_EV, PTDF_matrix, ...
                                P_Line_Base, P_Line_Max, N_bus, N_line)
    % --- 与 nonlinear_constraints_ga.m 相同的变量解码部分 ---
    num_vars_ac = num_ac * T;
    num_vars_ev = num_ev * T;
    N_vars_total = num_vars_ac + num_vars_ev;

    if isempty(x) || (N_vars_total > 0 && length(x) ~= N_vars_total)
        % (处理无效x的逻辑，同原函数)
        num_expected_constraints = T + 2 + 2 * N_line * T; % Power + SDCI + Rho + Line Limits (upper+lower)
        if T == 0; num_expected_constraints = 2 + 2*N_line*T; end
        c_ineq = ones(num_expected_constraints, 1) * 1e9;
        c_eq = [];
        return;
    end

    U_ac = zeros(num_ac, T);
    if num_ac > 0
        u_ac_flat = x(1:num_vars_ac);
        U_ac = reshape(logical(u_ac_flat), T, num_ac)';
    end

    U_ev = zeros(num_ev, T);
    if num_ev > 0
        u_ev_flat = x(num_vars_ac+1 : N_vars_total);
        U_ev = reshape(logical(u_ev_flat), T, num_ev)';
    end

    c_ineq = [];
    c_eq = [];
    eps_val = 1e-6;

    % --- 1. 计算原有的功率平衡、SDCI、Rho约束 ---
    c_ineq_power = zeros(T, 1);
    current_avg_P_ac_t = zeros(T, 1); current_n_ac_t = zeros(T, 1);
    current_avg_P_ev_t = zeros(T, 1); current_n_ev_t = zeros(T, 1);
    Delta_P_inj_k_all_t = zeros(N_bus, T); % *** 新增：存储所有节点的注入功率 ***

    for t_idx = 1:T
        % (计算 aggregated_ac_power_t, aggregated_ev_power_t,
        %  current_n_ac_t, current_avg_P_ac_t,
        %  current_n_ev_t, current_avg_P_ev_t 的逻辑同原函数)
        aggregated_ac_power_t = 0;
        if num_ac > 0
            aggregated_ac_power_t = sum(U_ac(:,t_idx) .* P_ac_ind_t(:,t_idx));
            current_n_ac_t(t_idx) = sum(U_ac(:,t_idx));
            current_avg_P_ac_t(t_idx) = aggregated_ac_power_t / (current_n_ac_t(t_idx) + eps_val);
        end
        aggregated_ev_power_t = 0;
        if num_ev > 0
            aggregated_ev_power_t = sum(U_ev(:,t_idx) .* P_ev_ind_t(:,t_idx));
            current_n_ev_t(t_idx) = sum(U_ev(:,t_idx));
            current_avg_P_ev_t(t_idx) = aggregated_ev_power_t / (current_n_ev_t(t_idx) + eps_val);
        end
        total_aggregated_power_t = aggregated_ac_power_t + aggregated_ev_power_t;
        c_ineq_power(t_idx) = P_req_t(t_idx) - total_aggregated_power_t; % P_req - P_supply <= 0

        % *** 新增：计算节点注入功率 Delta_Pinj_k(t) *** 
        Delta_P_inj_k_current_t = zeros(N_bus, 1);
        if num_ac > 0
            for k_bus = 1:N_bus
                 ac_at_bus_k_indices = find(Location_AC == k_bus);
                 if ~isempty(ac_at_bus_k_indices)
                    Delta_P_inj_k_current_t(k_bus) = Delta_P_inj_k_current_t(k_bus) + ...
                        sum(U_ac(ac_at_bus_k_indices, t_idx) .* P_ac_ind_t(ac_at_bus_k_indices, t_idx));
                 end
            end
        end
         if num_ev > 0
            for k_bus = 1:N_bus
                 ev_at_bus_k_indices = find(Location_EV == k_bus);
                 if ~isempty(ev_at_bus_k_indices)
                    Delta_P_inj_k_current_t(k_bus) = Delta_P_inj_k_current_t(k_bus) + ...
                        sum(U_ev(ev_at_bus_k_indices, t_idx) .* P_ev_ind_t(ev_at_bus_k_indices, t_idx));
                 end
            end
        end
        Delta_P_inj_k_all_t(:, t_idx) = Delta_P_inj_k_current_t;
        % *** 结束新增节点注入计算 ***
    end
    c_ineq = [c_ineq; c_ineq_power];

    SDCI_optimized = calculateSDCI(current_n_ac_t, current_n_ev_t, current_avg_P_ac_t, current_avg_P_ev_t);
    c_ineq_sdci = SDCI_optimized - SDCI_orig; % SDCI_opt - SDCI_target <= 0 [cite: 28]
    c_ineq = [c_ineq; c_ineq_sdci];

    rho_optimized = calculateSpearmanRho(current_n_ac_t, current_avg_P_ac_t, current_n_ev_t, current_avg_P_ev_t);
    c_ineq_rho = rho_optimized - rho_orig;   % rho_opt - rho_target <= 0 [cite: 31]
    c_ineq = [c_ineq; c_ineq_rho];

    % --- 2. *** 新增：计算网络潮流约束 *** ---
    if N_line > 0 && N_bus > 0
        c_ineq_line_upper = zeros(N_line * T, 1);
        c_ineq_line_lower = zeros(N_line * T, 1);
        idx_constraint_line = 1;

        for t_idx = 1:T
            % 计算线路潮流变化 Delta_PL_l(t) [cite: 41]
            Delta_P_Line_t = PTDF_matrix * Delta_P_inj_k_all_t(:, t_idx);

            % 计算最终线路潮流 P_L_l(t)
            P_Line_Final_t = P_Line_Base(:, t_idx) + Delta_P_Line_t;

            % 计算线路约束违反值 [cite: 45]
            % Upper bound: P_Line_Final_t(l) - P_Line_Max(l) <= 0
            violations_upper = P_Line_Final_t - P_Line_Max;
            % Lower bound: -P_Line_Final_t(l) - P_Line_Max(l) <= 0
            violations_lower = -P_Line_Final_t - P_Line_Max;

            % 存储约束违反值 (按时间顺序展开)
            c_ineq_line_upper(idx_constraint_line : idx_constraint_line + N_line - 1) = violations_upper;
            c_ineq_line_lower(idx_constraint_line : idx_constraint_line + N_line - 1) = violations_lower;
            idx_constraint_line = idx_constraint_line + N_line;
        end
        c_ineq = [c_ineq; c_ineq_line_upper; c_ineq_line_lower];
    end
    % --- 结束新增网络潮流约束 ---

end