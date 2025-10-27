function fitness = objectiveFunction_up(x, n_AC_max, n_EV_max, deltaP_AC, deltaP_EV, T, theta_1, theta_2, lambda_1, lambda_2_plus, Delta_P_grid_max)
    % 分解决策变量
    n_AC = x(1:T);          % 空调参与数量 [1×T]
    n_EV = x(T+1:2*T);      % EV参与数量 [1×T]
    
    % 确保deltaP为行向量（与n_AC/n_EV维度匹配）
    deltaP_AC = deltaP_AC(:).';  % [1×T]
    deltaP_EV = deltaP_EV(:).';
    
    % 计算调节能力
    p_AC_plus = n_AC .* deltaP_AC;  % [1×T]
    p_EV_plus = n_EV .* deltaP_EV;  % [1×T]
    
    % ==== 1. 计算归一化调节收益 ====
    total_power = sum(p_AC_plus) + sum(p_EV_plus);
    normalized_reward = total_power / (sum(n_AC_max*deltaP_AC)+sum(n_EV_max*deltaP_EV));
    
   % ==== 2. 动态SDCI惩罚（新增阈值0.2）====
    SDCI_plus = calculateSDCI(n_AC, n_EV, p_AC_plus, p_EV_plus);
    SDCI_penalty = lambda_1 * max(SDCI_plus - theta_1, 0);  % 硬约束：SDCI <= 0.2

    % ==== 3. 非对称相关性处理（新增rho约束）====
    rho_plus = calculateSpearmanRho(n_AC, p_AC_plus, n_EV, p_EV_plus);
    rho_threshold = theta_2;  % 硬约束：rho <= -0.5
    if rho_plus < rho_threshold
        rho_term = lambda_2_plus * abs(rho_plus);  % 达标奖励
    else
        rho_term = 1e6 * (rho_plus - rho_threshold);  % 未达标惩罚（极大值）
    end
    
    % ==== 综合目标函数 ====
    fitness = normalized_reward - SDCI_penalty - rho_term;
    
    % ==== 约束条件处理 ====
    if any(n_AC > n_AC_max) || any(n_EV > n_EV_max)
        fitness = -inf;
    end
end
