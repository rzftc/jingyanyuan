function fitness = objectiveFunction_up_aver(x, n_AC_max, n_EV_max, deltaP_AC, deltaP_EV, T, ...
                        theta_1, theta_3, lambda_1, lambda_2)
    % 目标函数计算聚合商调节能力与互补性指标的综合适应度
    % 
    % 输入参数：
    %   x            - 决策变量向量 [2T×1]
    %   n_AC_max     - 空调最大可调度数量
    %   n_EV_max     - 电动汽车最大可调度数量
    %   deltaP_AC    - 单台空调调节能力向量 [1×T]
    %   deltaP_EV    - 单台电动汽车调节能力向量 [1×T]
    %   T            - 时间窗口长度
    %   theta_1      - SDCI⁺阈值上限（典型值0.2）
    %   theta_3      - Spearman相关系数阈值上限（典型值-0.3）
    %   lambda_1     - SDCI⁺惩罚系数
    %   lambda_2     - 相关系数惩罚系数
    %
    % 输出参数：
    %   fitness      - 综合适应度值

    %% 决策变量解析
    n_AC = x(1:T);              % 各时段空调调度数量
    n_EV = x(T+1:2*T);          % 各时段电动汽车调度数量
    
    % 确保调节能力向量为行向量
    deltaP_AC = deltaP_AC(:).';  
    deltaP_EV = deltaP_EV(:).';
    
    % 计算总调节能力
    p_AC = n_AC .* deltaP_AC;   % 空调总调节功率
    p_EV = n_EV .* deltaP_EV;   % 电动汽车总调节功率

    %% 核心指标计算
    % 标准化调节收益计算
    total_power = sum(p_AC) + sum(p_EV);
    max_power = sum(n_AC_max*deltaP_AC) + sum(n_EV_max*deltaP_EV);
    normalized_reward = total_power / (max_power + eps);  % 添加eps防止除零
    
    % 互补性指标计算
    SDCI_plus = calculateSDCI(n_AC, n_EV, deltaP_AC, deltaP_EV);
    rho_plus = calculateSpearmanRho(n_AC, deltaP_AC, n_EV, deltaP_EV);

    %% 目标函数构建
    base_fitness = normalized_reward ...                % 基础收益项
                   - lambda_1 * SDCI_plus ...           % SDCI⁺惩罚项
                   - lambda_2 * rho_plus;               % 相关系数惩罚项

    %% 硬约束处理
    penalty = 0;
    
    % 约束1: SDCI⁺ ≤ θ₁
    if SDCI_plus > theta_1
        penalty = penalty + 1e6 * (SDCI_plus - theta_1);
    end
    
    % 约束2: ρ ≤ θ₃（强制负相关性）
    if rho_plus > theta_3
        penalty = penalty + 1e6 * (rho_plus - theta_3);
    end
    
    % 约束3: 设备数量限制
    if any(n_AC > n_AC_max) || any(n_EV > n_EV_max)
        penalty = penalty + 1e6;
    end
    
    %% 综合适应度计算
    fitness = base_fitness - penalty;
end
