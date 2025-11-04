% create_dummy_optimization_inputs.m
function OptInputs = create_dummy_optimization_inputs(SimData)
    % (v2) 生成基于 SimData 动态缩放的虚拟优化输入
    % !!! 用户必须用真实数据替换此文件中的内容 !!!

    % --- 用户可调参数 ---
    % 设置需求为总可用调节能力的百分比
    demand_ratio = 0.3; % (例如, 30%)
    % 设置线路容量为总调节能力峰值的百分比
    line_limit_ratio = 0.1; % (例如, 每条线路承载峰值的10%)
    % --------------------

    % 从 SimData 提取信息
    nAC = SimData.nAC;
    nEV = SimData.nEV;
    T = SimData.T;
    N_total = SimData.N_total;
    
    %% 1. 动态虚拟电网需求 P_req(t)
    % 计算每个时间步t的总可用调节潜力 (1 x T)
    total_available_potential_t = sum(SimData.p_AC, 1) + sum(SimData.p_EV, 1);
    
    % P_req(t) 是总可用潜力的 30%
    OptInputs.P_req = total_available_potential_t * demand_ratio; % 1xT
    
    % 找出系统潜力的峰值
    peak_system_potential = max(total_available_potential_t);
    if peak_system_potential == 0
        warning('加载的总调节潜力为0，电网需求将全为0。');
        peak_system_potential = 1; % 避免除零
    end

    %% 2. 虚拟调节成本 c_i
    % (这部分保持不变，因为成本与潜力无关)
    OptInputs.DeviceCosts.c_AC = 0.1 + 0.1 * rand(nAC, 1); % nAC x 1
    OptInputs.DeviceCosts.c_EV = 0.3 + 0.2 * rand(nEV, 1); % nEV x 1
    
    %% 3. 动态虚拟网络拓扑
    % (假设的节点和线路数量保持不变)
    N_nodes = 10;
    N_lines = 12;
    
    % 3a. 设备到节点的映射 (DeviceNodes)
    OptInputs.Network.DeviceNodes = randi(N_nodes, N_total, 1); % N_total x 1
    
    % 3b. PTDF 矩阵
    OptInputs.Network.PTDF = (rand(N_lines, N_nodes) - 0.5) * 0.2; % L x K
    
    % 3c. 动态线路潮流限值 (LineLimits)
    % 假设每条线路的限值是系统总潜力峰值的 10%
    line_limit_value = peak_system_potential * line_limit_ratio;
    OptInputs.Network.LineLimits = repmat(line_limit_value, N_lines, 1); % L x 1 (kW)
    
    % 3d. 动态基础潮流 (BaseFlow)
    % 假设基础潮流是 *新* 线路限值的 10% 到 40%
    t_hours = linspace(0, 24, T); % 假设时间轴
    base_load_profile = 0.1 + 0.3 * sin(t_hours * pi / 24).^2;
    OptInputs.Network.BaseFlow = OptInputs.Network.LineLimits .* base_load_profile; % L x T
    
end