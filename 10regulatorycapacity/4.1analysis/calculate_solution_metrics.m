% calculate_solution_metrics.m
function Metrics = calculate_solution_metrics(u_optimal, SimData, name)
    % 计算给定决策 u_optimal 下的聚合曲线和互补性指标
    
    % u_optimal 是 (nAC+nEV) x T (logical)
    % SimData.p_AC 是 (nAC x T)
    % SimData.p_EV 是 (nEV x T)

    % 1. 提取 AC 和 EV 各自的决策
    u_AC = u_optimal(1:SimData.nAC, :);
    u_EV = u_optimal(SimData.nAC + 1 : end, :);
    
    % 2. 计算聚合调节功率 P_AC_opt(t) 和 P_EV_opt(t)
    % P_AC_opt(t) = sum_i(u_AC(i,t) * p_AC(i,t))
    P_AC_agg = sum(u_AC .* SimData.p_AC, 1)'; % 聚合为 T x 1 列向量
    P_EV_agg = sum(u_EV .* SimData.p_EV, 1)'; % 聚合为 T x 1 列向量
    
    % 3. 准备调用库函数
    % 库函数需要 n(t) 和 deltaP(t)
    % 我们假设 deltaP_AC 和 deltaP_EV 为 1，n_AC 和 n_EV 等于聚合功率
    
    n_AC_dummy = P_AC_agg;
    deltaP_AC_dummy = ones(SimData.T, 1);
    n_EV_dummy = P_EV_agg;
    deltaP_EV_dummy = ones(SimData.T, 1);
    
    % 4. 调用库函数
    % 确保 calculateSDCI 和 calculateSpearmanRho 在 MATLAB 路径中
    try
        % (v2 修正：确保 calculateSDCI 和 calculateSpearmanRho 的输入参数顺序)
        % (假设您的库函数签名是: (n_AC, deltaP_AC, n_EV, deltaP_EV))
        Metrics.SDCI = calculateSDCI(n_AC_dummy, deltaP_AC_dummy, n_EV_dummy, deltaP_EV_dummy);
    catch ME_sdci
        fprintf('错误: 调用 calculateSDCI 失败: %s\n', ME_sdci.message);
        Metrics.SDCI = NaN;
    end
    
    try
        Metrics.Rho = calculateSpearmanRho(n_AC_dummy, deltaP_AC_dummy, n_EV_dummy, deltaP_EV_dummy);
    catch ME_rho
        fprintf('错误: 调用 calculateSpearmanRho 失败: %s\n', ME_rho.message);
        Metrics.Rho = NaN;
    end
    
    Metrics.P_AC_agg = P_AC_agg;
    Metrics.P_EV_agg = P_EV_agg;
    Metrics.Name = name;
end