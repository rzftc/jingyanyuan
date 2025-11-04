function fitness = objectiveFunction_down_dt(x, n_AC_max, n_EV_max, deltaP_AC, deltaP_EV, T, ...
                                       theta_2, theta_3, lambda_1, lambda_2)
    % 分解决策变量（保持原样）
    n_AC = round(x(1:T));
    n_EV = round(x(T+1:2*T));
    
    % 初始化总调节能力
    p_AC = zeros(1,T);
    p_EV = zeros(1,T);
    
    % 使用parfor并行计算各时段功率
    parfor t = 1:T
        % 保持原有排序和求和逻辑
        sorted_AC = sort(deltaP_AC(:,t), 'descend');
        p_AC(t) = sum(sorted_AC(1:min(n_AC(t), end)));
        
        sorted_EV = sort(deltaP_EV(:,t), 'descend');
        p_EV(t) = sum(sorted_EV(1:min(n_EV(t), end)));
    end

    %% ===== 保持原有核心指标计算逻辑 =====
    total_power = sum(p_AC) + sum(p_EV);
    max_power_AC = sum(maxk(deltaP_AC(:), n_AC_max*T));
    max_power_EV = sum(maxk(deltaP_EV(:), n_EV_max*T));
    normalized_reward = total_power / (max_power_AC + max_power_EV + eps);
    
    avg_deltaP_AC = p_AC ./ (n_AC + eps);
    avg_deltaP_EV = p_EV ./ (n_EV + eps);
    SDCI_minus = calculateSDCI(n_AC', n_EV', avg_deltaP_AC', avg_deltaP_EV');
    rho_minus = calculateSpearmanRho(n_AC', avg_deltaP_AC', n_EV', avg_deltaP_EV');
    
    %% ===== 保持原有目标函数构建逻辑 =====
    base_fitness = normalized_reward - lambda_1*SDCI_minus - lambda_2*rho_minus;
    
    %% ===== 保持原有约束处理逻辑 =====
    penalty = 0;
    if SDCI_minus > theta_2
        penalty = penalty + 1e6*(SDCI_minus - theta_2);
    end
    if rho_minus > theta_3
        penalty = penalty + 1e6*(rho_minus - theta_3);
    end
    if any(n_AC > n_AC_max) || any(n_EV > n_EV_max)
        penalty = penalty + 1e6;
    end
    
    fitness = base_fitness - penalty;
end