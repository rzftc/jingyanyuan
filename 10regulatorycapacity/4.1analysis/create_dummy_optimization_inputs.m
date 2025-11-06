function OptInputs = create_dummy_optimization_inputs(SimData)
    % create_dummy_optimization_inputs (v5 - 使用预加载聚合数据)

    % --- 用户可调参数 ---
    demand_ratio = 0.4;      % 需求占总潜力的比例
    line_limit_ratio = 0.15; % 线路限值占系统峰值的比例
    % --------------------

    nAC = SimData.nAC; nEV = SimData.nEV;
    T = SimData.T; N_total = SimData.N_total;
    
    %% 1. 动态虚拟电网需求 P_req(t)
    % [修改] 直接使用加载的聚合数据 (保留了原始符号：上调+, 下调-)
    total_potential_t = SimData.P_AC_agg_loaded + SimData.P_EV_agg_loaded;
    
    % 需求曲线跟随总潜力的符号和趋势，并加入随机波动
    random_fluctuation = 0.8 + 0.4 * rand(1, T); 
    OptInputs.P_req = total_potential_t * demand_ratio .* random_fluctuation;
    
    % 计算系统能力的幅度峰值，用于设置正值的网络参数 (线路容量必须为正)
    peak_magnitude = max(abs(total_potential_t));
    if peak_magnitude < 1e-9, peak_magnitude = 100; end

    %% 2. 虚拟调节成本 (始终为正)
    OptInputs.DeviceCosts.c_AC = 0.2 + 0.1 * rand(nAC, 1); 
    OptInputs.DeviceCosts.c_EV = 0.15 + 0.1 * rand(nEV, 1); 
    
    %% 3. 虚拟网络拓扑
    N_nodes = 10; N_lines = 12;
    OptInputs.Network.DeviceNodes = randi(N_nodes, N_total, 1);
    OptInputs.Network.PTDF = (rand(N_lines, N_nodes) - 0.5) * 0.2; 
    
    % 线路容量 (使用幅度)
    limit_val = peak_magnitude * line_limit_ratio;
    OptInputs.Network.LineLimits = repmat(limit_val, N_lines, 1);
    
    % 基础潮流
    base_flow = 0.2 * sin(linspace(0, 2*pi, T)) + 0.1 * randn(1, T);
    OptInputs.Network.BaseFlow = OptInputs.Network.LineLimits .* base_flow;
end