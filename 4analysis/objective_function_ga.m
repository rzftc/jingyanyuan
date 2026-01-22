function fitness = objective_function_ga(x, num_ac, num_ev, P_ac_ind_t, P_ev_ind_t, C_ac_ind, C_ev_ind, P_req_t, T, SDCI_orig, rho_orig)
    % OBJECTIVE_FUNCTION_GA (精简版)
    % 适应度函数，计算成本并结合约束惩罚
    % 假设所有输入参数（包括x）均由调用者确保有效性

    num_vars_ac = num_ac * T;
    num_vars_ev = num_ev * T;
    N_vars_expected_here = num_vars_ac + num_vars_ev;

    % 核心检查：如果GA传递的 x 为空或长度与期望总变量数不匹配，则返回极差的适应度值
    % 这是功能性的，因为GA在某些情况下可能传递无效x，函数需要能处理
    if isempty(x) || (N_vars_expected_here > 0 && length(x) ~= N_vars_expected_here)
        fitness = 1e20; % 返回一个非常大的适应度值
        return;
    end

    % --- 安全地提取决策变量 ---
    if num_vars_ac > 0
        u_ac_flat = x(1:num_vars_ac);
    else
        u_ac_flat = [];
    end

    if num_vars_ev > 0
        start_idx_ev = num_vars_ac + 1;
        u_ev_flat = x(start_idx_ev : N_vars_expected_here); % 使用 N_vars_expected_here
    else
        u_ev_flat = [];
    end
    
    % --- 解码决策变量 ---
    U_ac = zeros(num_ac, T);
    if num_ac > 0 && ~isempty(u_ac_flat)
        U_ac = reshape(logical(u_ac_flat), T, num_ac)';
    end
    U_ev = zeros(num_ev, T);
    if num_ev > 0 && ~isempty(u_ev_flat)
        U_ev = reshape(logical(u_ev_flat), T, num_ev)';
    end

    % --- 计算总成本 ---
    total_cost = 0;
    if num_ac > 0 % 仅当有AC设备时计算AC成本
        cost_ac_matrix = U_ac .* P_ac_ind_t .* repmat(C_ac_ind(:), 1, T);
        total_cost = total_cost + sum(cost_ac_matrix(:));
    end
    if num_ev > 0 % 仅当有EV设备时计算EV成本
        cost_ev_matrix = U_ev .* P_ev_ind_t .* repmat(C_ev_ind(:), 1, T);
        total_cost = total_cost + sum(cost_ev_matrix(:));
    end
    
    fitness = total_cost;

    % --- 获取约束违反情况并施加惩罚 (调用nonlinear_constraints_ga) ---
    [c_ineq_vals, ~] = nonlinear_constraints_ga(x, num_ac, num_ev, P_ac_ind_t, P_ev_ind_t, P_req_t, T, SDCI_orig, rho_orig);

    penalty_factor_power = max(1e6, abs(total_cost * 100) + 1e3); 
    penalty_factor_sdci_rho = max(1e4, abs(total_cost * 10) + 1e2);

    num_power_constraints = T;
    if T == 0; num_power_constraints = 0; end

    % 功率平衡约束惩罚
    if length(c_ineq_vals) >= num_power_constraints && num_power_constraints > 0
        penalty_power = sum(max(0, c_ineq_vals(1:num_power_constraints)).^2); 
        fitness = fitness + penalty_power * penalty_factor_power;
    end

    % SDCI 约束惩罚
    idx_sdci_constraint = num_power_constraints + 1;
    if length(c_ineq_vals) >= idx_sdci_constraint 
        penalty_sdci = max(0, c_ineq_vals(idx_sdci_constraint))^2;
        fitness = fitness + penalty_sdci * penalty_factor_sdci_rho;
    end

    % Rho 约束惩罚
    idx_rho_constraint = num_power_constraints + 2;
    if length(c_ineq_vals) >= idx_rho_constraint
        penalty_rho = max(0, c_ineq_vals(idx_rho_constraint))^2;
        fitness = fitness + penalty_rho * penalty_factor_sdci_rho;
    end
end