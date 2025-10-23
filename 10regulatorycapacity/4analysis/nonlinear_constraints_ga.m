function [c_ineq, c_eq] = nonlinear_constraints_ga(x, num_ac, num_ev, P_ac_ind_t, P_ev_ind_t, P_req_t, T, SDCI_orig, rho_orig)
    % NONLINEAR_CONSTRAINTS_GA (精简版)
    % 计算非线性约束值，被适应度函数调用
    % 假设所有输入参数（包括x）均由调用者(objective_function_ga)确保有效性

    num_vars_ac = num_ac * T;
    num_vars_ev = num_ev * T;
    N_vars_expected_here = num_vars_ac + num_vars_ev;

    % 核心检查：如果x的长度与期望总变量数不匹配 (由适应度函数调用时应已保证)
    % 这个检查在这里作为双重保险，或者如果此函数可能被其他途径调用
    if isempty(x) || (N_vars_expected_here > 0 && length(x) ~= N_vars_expected_here)
        num_constraints = T + 2; if T == 0; num_constraints = 2; end; if num_constraints == 0; num_constraints=1; end
        c_ineq = ones(num_constraints, 1) * 1e9; % 返回大的约束违反值
        c_eq = [];
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
        u_ev_flat = x(start_idx_ev : N_vars_expected_here);
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
    
    c_ineq = [];
    c_eq = []; 
    eps_val = 1e-6;

    % --- 计算约束值 ---
    c_ineq_power = zeros(T, 1);
    current_avg_P_ac_t = zeros(T, 1); 
    current_n_ac_t = zeros(T, 1);     
    current_avg_P_ev_t = zeros(T, 1); 
    current_n_ev_t = zeros(T, 1);     

    for t_idx = 1:T 
        aggregated_ac_power_t = 0;
        if num_ac > 0 % 仅当有AC设备时计算
            aggregated_ac_power_t = sum(U_ac(:,t_idx) .* P_ac_ind_t(:,t_idx));
            current_n_ac_t(t_idx) = sum(U_ac(:,t_idx));
            if current_n_ac_t(t_idx) > 0
                current_avg_P_ac_t(t_idx) = aggregated_ac_power_t / current_n_ac_t(t_idx);
            % else avg_P_ac is 0 (initialized)
            end
        end

        aggregated_ev_power_t = 0;
        if num_ev > 0 % 仅当有EV设备时计算
            aggregated_ev_power_t = sum(U_ev(:,t_idx) .* P_ev_ind_t(:,t_idx));
            current_n_ev_t(t_idx) = sum(U_ev(:,t_idx));
            if current_n_ev_t(t_idx) > 0
                current_avg_P_ev_t(t_idx) = aggregated_ev_power_t / current_n_ev_t(t_idx);
            % else avg_P_ev is 0 (initialized)
            end
        end
        
        total_aggregated_power_t = aggregated_ac_power_t + aggregated_ev_power_t;
        c_ineq_power(t_idx) = P_req_t(t_idx) - total_aggregated_power_t; % P_req - P_supply <= 0
    end
    c_ineq = [c_ineq; c_ineq_power];

    % SDCI 约束
    SDCI_optimized = calculateSDCI(current_n_ac_t, current_n_ev_t, current_avg_P_ac_t, current_avg_P_ev_t);
    c_ineq_sdci = SDCI_optimized - SDCI_orig; % SDCI_opt - SDCI_target <= 0
    c_ineq = [c_ineq; c_ineq_sdci];
    
    % Rho 约束
    rho_optimized = calculateSpearmanRho(current_n_ac_t, current_avg_P_ac_t, current_n_ev_t, current_avg_P_ev_t);
    c_ineq_rho = rho_optimized - rho_orig;   % rho_opt - rho_target <= 0
    c_ineq = [c_ineq; c_ineq_rho];
end