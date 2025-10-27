function fitness = objectiveFunction_down(x, n_AC_max, n_EV_max, deltaP_AC, deltaP_EV, T, theta_1, theta_2, lambda_1, lambda_2, Delta_P_grid_max)
    % 目标函数计算下调节场景的适应度值
    % 输入参数：
    % x - 决策变量向量
    % n_AC_max - 空调最大可调度数量
    % n_EV_max - 电动汽车最大可调度数量
    % deltaP_AC - 单台空调功率调节量（kW）
    % deltaP_EV - 单台电动汽车功率调节量（kW）
    % T - 时间步长总数
    % theta_1, theta_2 - 约束阈值参数
    % lambda_1, lambda_2 - 惩罚系数
    % Delta_P_grid_max - 电网功率变化上限

    % 解析决策变量
    n_AC = x(1:T);          % 各时段空调参与数量（T维行向量）
    n_EV = x(T+1:2*T);      % 各时段电动汽车参与数量（T维行向量）
    
    % 确保功率调节量为行向量以匹配维度
    deltaP_AC = deltaP_AC(:).';  % 转换为行向量 [1×T]
    deltaP_EV = deltaP_EV(:).';  % 转换为行向量 [1×T]
    
    % 计算总调节能力（下调方向）
    p_AC_minus = n_AC .* deltaP_AC;  % 空调总调节功率 [1×T]
    p_EV_minus = n_EV .* deltaP_EV;  % 电动汽车总调节功率 [1×T]
    
    % 1. 归一化调节收益计算
    total_power = sum(p_AC_minus) + sum(p_EV_minus);  % 总调节功率
    denominator = sum(n_AC_max * deltaP_AC) + sum(n_EV_max * deltaP_EV);  % 最大可能调节功率
    if denominator <= 0
        normalized_reward = 0;  % 避免除零错误
    else
        normalized_reward = total_power / denominator;  % 归一化收益
    end
    
    % 2. 动态SDCI约束惩罚项（硬约束SDCI <= theta_1）
    SDCI_minus = calculateSDCI(n_AC, n_EV, p_AC_minus, p_EV_minus);
    SDCI_penalty = lambda_1 * max(SDCI_minus - theta_1, 0);  % 超阈值惩罚

    % 3. 非对称相关性约束处理（硬约束rho <= theta_2）
    rho_minus = calculateSpearmanRho(n_AC, p_AC_minus, n_EV, p_EV_minus);
    rho_threshold = theta_2;  % 相关系数阈值
    if rho_minus < rho_threshold
        rho_term = lambda_2 * abs(rho_minus);  % 达标奖励项
    else
        rho_term = 1e6 * (rho_minus - rho_threshold);  % 未达标惩罚项
    end
    
    % 综合适应度计算
    fitness = normalized_reward - SDCI_penalty - rho_term;
    
    % 约束条件处理
    if any(n_AC > n_AC_max) || any(n_EV > n_EV_max)  % 设备数量越界检查
        fitness = -inf;  % 违反约束则标记为不可行解
    end
end
