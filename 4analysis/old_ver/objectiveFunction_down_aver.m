%% objectiveFunction_down2.m
function fitness = objectiveFunction_down_aver(x, n_AC_max, n_EV_max, deltaP_AC, deltaP_EV, T, ...
                                            theta_2, theta_3, lambda_1, lambda_2, Delta_P_grid_max)
    % 参数说明：
    % theta_2 - SDCI⁻上限阈值（典型值0.2）
    % theta_3 - rho上限阈值（典型值-0.3）

    % 分解决策变量
    n_AC = x(1:T);          % 空调数量
    n_EV = x(T+1:2*T);      % 电动汽车数量
    deltaP_AC = deltaP_AC(:).';  % 空调调节能力（1×T）
    deltaP_EV = deltaP_EV(:).';  % 电动汽车调节能力（1×T）

    % 计算实际调节能力（单台设备能力 × 数量）
    p_AC = n_AC .* deltaP_AC;
    p_EV = n_EV .* deltaP_EV;
    
    %% ===== 核心指标计算 =====
    % 调节收益计算
    total_power = sum(p_AC) + sum(p_EV);    % 总调节能力
    max_power = sum(n_AC_max * deltaP_AC) + sum(n_EV_max * deltaP_EV);  % 最大调节能力
    normalized_reward = total_power / (max_power + eps);  % 标准化调节收益
    
    % 计算互补性指标（注意传入单台设备能力）
    SDCI_minus = calculateSDCI(n_AC, n_EV, deltaP_AC, deltaP_EV);  % 计算SDCI⁻
    rho_minus = calculateSpearmanRho(n_AC, deltaP_AC, n_EV, deltaP_EV);  % 计算ρ值
    
    %% ===== 目标函数构建 =====
    base_fitness = normalized_reward ...                      % 标准化调节收益
                   - lambda_1 * SDCI_minus ...              % SDCI惩罚项
                   - lambda_2 * rho_minus;                 % 相关系数惩罚项
    
    %% ===== 硬约束处理 =====
    penalty = 0;  % 初始化惩罚项
    
    % 约束1：SDCI⁻ ≤ θ₂
    if SDCI_minus > theta_2
        penalty = penalty + 1e6 * (SDCI_minus - theta_2);  % SDCI⁻大于上限时惩罚
    end
    
    % 约束2：ρ ≤ θ₃（强制负相关）
    if rho_minus > theta_3
        penalty = penalty + 1e6 * (rho_minus - theta_3);  % ρ大于上限时惩罚
    end
    
    % 约束3：设备数量限制
    if any(n_AC > n_AC_max) || any(n_EV > n_EV_max)
        penalty = penalty + 1e6;  % 如果设备数量超出限制，进行惩罚
    end
    
    % 计算最终的适应度值
    fitness = base_fitness - penalty;
end
