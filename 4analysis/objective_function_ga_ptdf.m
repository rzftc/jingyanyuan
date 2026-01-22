% objective_function_ga_ptdf.m
function fitness = objective_function_ga_ptdf(x, num_ac, num_ev, ...
                                P_ac_ind_t, P_ev_ind_t, ...
                                C_ac_ind, C_ev_ind, P_req_t, T, ...
                                SDCI_orig, rho_orig, ...
                                Location_AC, Location_EV, PTDF_matrix, ... % 新增
                                P_Line_Base, P_Line_Max, N_bus, N_line)    % 新增
    % --- 与 objective_function_ga.m 相同的成本计算部分 ---
    num_vars_ac = num_ac * T;
    num_vars_ev = num_ev * T;
    N_vars_total = num_vars_ac + num_vars_ev;

    if isempty(x) || (N_vars_total > 0 && length(x) ~= N_vars_total)
        fitness = 1e20;
        return;
    end
    % (解码 U_ac, U_ev 的逻辑同原函数)
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

    total_cost = 0;
    if num_ac > 0
        cost_ac_matrix = U_ac .* P_ac_ind_t .* repmat(C_ac_ind(:), 1, T);
        total_cost = total_cost + sum(cost_ac_matrix(:));
    end
    if num_ev > 0
        cost_ev_matrix = U_ev .* P_ev_ind_t .* repmat(C_ev_ind(:), 1, T);
        total_cost = total_cost + sum(cost_ev_matrix(:));
    end
    fitness = total_cost; % 目标函数仍然是最小化成本 [cite: 18]

    % --- 调用新的约束函数计算约束违反和惩罚 ---
    [c_ineq_vals, ~] = nonlinear_constraints_ga_ptdf(x, num_ac, num_ev, ...
                                P_ac_ind_t, P_ev_ind_t, P_req_t, T, ...
                                SDCI_orig, rho_orig, ...
                                Location_AC, Location_EV, PTDF_matrix, ...
                                P_Line_Base, P_Line_Max, N_bus, N_line); % 传递新参数

    % --- 惩罚项计算 (同原函数，但需要考虑新增的约束) ---
    penalty_factor_power = max(1e6, abs(total_cost * 100) + 1e3);
    penalty_factor_sdci_rho = max(1e4, abs(total_cost * 10) + 1e2);
    penalty_factor_line = max(1e7, abs(total_cost * 200) + 1e4); % 对线路约束设置更高的惩罚因子

    % 功率平衡约束惩罚 [cite: 25]
    num_power_constraints = T; if T == 0; num_power_constraints = 0; end
    if length(c_ineq_vals) >= num_power_constraints && num_power_constraints > 0
        penalty_power = sum(max(0, c_ineq_vals(1:num_power_constraints)).^2);
        fitness = fitness + penalty_power * penalty_factor_power;
    end

    % SDCI 约束惩罚 [cite: 28]
    idx_sdci_constraint = num_power_constraints + 1;
    if length(c_ineq_vals) >= idx_sdci_constraint
        penalty_sdci = max(0, c_ineq_vals(idx_sdci_constraint))^2;
        fitness = fitness + penalty_sdci * penalty_factor_sdci_rho;
    end

    % Rho 约束惩罚 [cite: 31]
    idx_rho_constraint = num_power_constraints + 2;
    if length(c_ineq_vals) >= idx_rho_constraint
        penalty_rho = max(0, c_ineq_vals(idx_rho_constraint))^2;
        fitness = fitness + penalty_rho * penalty_factor_sdci_rho;
    end

    % *** 新增：线路潮流约束惩罚 *** [cite: 45]
    idx_line_start = num_power_constraints + 3;
    num_line_constraints = 2 * N_line * T; % 上下限
    if length(c_ineq_vals) >= idx_line_start + num_line_constraints - 1 && num_line_constraints > 0
        penalty_line = sum(max(0, c_ineq_vals(idx_line_start : idx_line_start + num_line_constraints - 1)).^2);
        fitness = fitness + penalty_line * penalty_factor_line;
    end
    % *** 结束新增惩罚 ***
end