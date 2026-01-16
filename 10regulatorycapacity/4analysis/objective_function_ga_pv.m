function fitness = objective_function_ga_pv(x, num_ac, num_ev, num_pv, ...
                                            P_ac_ind_t, P_ev_ind_t, P_pv_ind_t, ...
                                            C_ac_ind, C_ev_ind, C_pv_ind, ...
                                            P_req_t, T, SDCI_orig, rho_orig)
    % OBJECTIVE_FUNCTION_GA_PV
    % 适应度函数，计算成本并结合约束惩罚（支持光伏PV）
    
    num_vars_ac = num_ac * T;
    num_vars_ev = num_ev * T;
    num_vars_pv = num_pv * T;
    N_vars_expected_here = num_vars_ac + num_vars_ev + num_vars_pv;

    % 核心检查
    if isempty(x) || (N_vars_expected_here > 0 && length(x) ~= N_vars_expected_here)
        fitness = 1e20; 
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
        end_idx_ev = num_vars_ac + num_vars_ev;
        u_ev_flat = x(start_idx_ev : end_idx_ev);
    else
        u_ev_flat = [];
    end
    
    if num_vars_pv > 0
        start_idx_pv = num_vars_ac + num_vars_ev + 1;
        u_pv_flat = x(start_idx_pv : N_vars_expected_here);
    else
        u_pv_flat = [];
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
    U_pv = zeros(num_pv, T);
    if num_pv > 0 && ~isempty(u_pv_flat)
        U_pv = reshape(logical(u_pv_flat), T, num_pv)';
    end

    % --- 计算总成本 ---
    total_cost = 0;
    if num_ac > 0 
        cost_ac_matrix = U_ac .* P_ac_ind_t .* repmat(C_ac_ind(:), 1, T);
        total_cost = total_cost + sum(cost_ac_matrix(:));
    end
    if num_ev > 0 
        cost_ev_matrix = U_ev .* P_ev_ind_t .* repmat(C_ev_ind(:), 1, T);
        total_cost = total_cost + sum(cost_ev_matrix(:));
    end
    if num_pv > 0 % 新增 PV 成本
        cost_pv_matrix = U_pv .* P_pv_ind_t .* repmat(C_pv_ind(:), 1, T);
        total_cost = total_cost + sum(cost_pv_matrix(:));
    end
    
    fitness = total_cost;

    % --- 计算约束违反 (集成版，不依赖外部文件) ---
    [c_ineq_vals] = calculate_constraints_internal(U_ac, U_ev, U_pv, P_ac_ind_t, P_ev_ind_t, P_pv_ind_t, P_req_t, T, SDCI_orig, rho_orig);

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

function [c_ineq] = calculate_constraints_internal(U_ac, U_ev, U_pv, P_ac_t, P_ev_t, P_pv_t, P_req, T, SDCI_ref, rho_ref)
    % 内部子函数：计算约束违反值
    % 1. 功率平衡: P_total >= P_req
    % 2. SDCI: SDCI_new >= SDCI_ref (保持互补性不恶化)
    % 3. Rho: Rho_new <= Rho_ref (保持相关性不恶化，或者根据具体需求，此处假设维持原状)

    c_ineq = zeros(T + 2, 1);
    
    % --- 1. 功率约束 ---
    P_total_supplied = zeros(T, 1);
    if ~isempty(U_ac); P_total_supplied = P_total_supplied + sum(U_ac .* P_ac_t, 1)'; end
    if ~isempty(U_ev); P_total_supplied = P_total_supplied + sum(U_ev .* P_ev_t, 1)'; end
    if ~isempty(U_pv); P_total_supplied = P_total_supplied + sum(U_pv .* P_pv_t, 1)'; end % 加入PV
    
    % 违反条件：需求 > 供给
    c_ineq(1:T) = P_req - P_total_supplied;

    % --- 准备 SDCI/Rho 计算数据 (仅基于 AC 和 EV) ---
    n_ac_t = sum(U_ac, 1)'; 
    n_ev_t = sum(U_ev, 1)';
    
    eps_val = 1e-6;
    avg_P_ac_t = zeros(T, 1);
    if size(U_ac,1) > 0
        agg_ac = sum(U_ac .* P_ac_t, 1)';
        avg_P_ac_t = agg_ac ./ (n_ac_t + eps_val);
    end
    
    avg_P_ev_t = zeros(T, 1);
    if size(U_ev,1) > 0
        agg_ev = sum(U_ev .* P_ev_t, 1)';
        avg_P_ev_t = agg_ev ./ (n_ev_t + eps_val);
    end

    % --- 2. SDCI 约束 ---
    % 假设调用 calculateSDCI.m
    try
        sdci_curr = calculateSDCI(n_ac_t, n_ev_t, avg_P_ac_t, avg_P_ev_t);
        % 约束：新 SDCI 不应比参考值低太多 (SDCI_ref - sdci_curr <= 0)
        c_ineq(T+1) = SDCI_ref - sdci_curr; 
    catch
        c_ineq(T+1) = 0; % 出错则忽略
    end

    % --- 3. Rho 约束 ---
    % 假设调用 calculateSpearmanRho.m
    try
        rho_curr = calculateSpearmanRho(n_ac_t, avg_P_ac_t, n_ev_t, avg_P_ev_t);
        % 约束：通常希望相关性越低越好(负相关)，或者维持原状。
        % 此处假设：rho_curr - rho_ref <= 0 (不希望相关性异常升高)
        % 或者根据 calculateSpearmanRho 的具体逻辑 (通常互补性强意味着 rho 趋向 -1)
        % 假设这里的 rho_origin 是原始状态，我们不希望它变大（变正相关）
        c_ineq(T+2) = rho_curr - rho_ref; 
    catch
        c_ineq(T+2) = 0;
    end
end