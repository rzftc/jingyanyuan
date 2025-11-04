% solve_times_step_greedy.m
function [u_opt_t, P_total_t, cost_t, is_feasible] = ...
    solve_times_step_greedy(p_devices_t, c_devices, P_req_t, Network)
    % 使用贪心算法求解单时刻、网络约束的背包问题
    %
    % 输入:
    %   p_devices_t: (N x 1) 当前时刻 t 的所有设备可用潜力
    %   c_devices:   (N x 1) 所有设备的单位调节成本 (元/kW)
    %   P_req_t:     (scalar) 当前时刻 t 的总电网需求
    %   Network:     (struct) 包含 PTDF, LineLimits, BaseFlow, DeviceNodes
    %
    % 输出:
    %   u_opt_t:     (N x 1) 0/1 决策变量 (logical)
    %   P_total_t:   (scalar) 优化后的总功率
    %   cost_t:      (scalar) 优化后的总成本
    %   is_feasible: (bool) 是否满足了 P_req_t

    N = length(p_devices_t);
    u_opt_t = false(N, 1); % 使用 logical 提高效率
    P_total_t = 0;
    cost_t = 0;
    
    if P_req_t <= 0
        is_feasible = true;
        return;
    end

    % --- 贪心策略 ---
    % 1. 计算每个设备的 "总成本" = 潜力 * 单位成本
    total_cost_per_device = p_devices_t .* c_devices;
    
    % 2. 找出潜力 > 0 且成本有限的设备
    valid_indices = find(p_devices_t > 1e-3 & isfinite(total_cost_per_device));
    if isempty(valid_indices)
        is_feasible = false;
        return;
    end
    
    % 3. 按总成本升序排序
    [~, sorted_order] = sort(total_cost_per_device(valid_indices), 'ascend');
    sorted_indices = valid_indices(sorted_order);
    
    % --- 网络约束准备 ---
    L = size(Network.PTDF, 1); % 线路数量
    current_delta_flow = zeros(L, 1); % 调节引起的潮流变化
    
    % 计算潮流安全裕度
    margin_upper = Network.LineLimits - Network.BaseFlow;
    margin_lower = -Network.LineLimits - Network.BaseFlow;

    % 4. 迭代选择设备
    for idx = sorted_indices' % 确保是列向量迭代
        
        p_i = p_devices_t(idx);
        node_k = Network.DeviceNodes(idx);
        
        % --- 检查增加该设备是否会违反网络约束 ---
        delta_flow_i = Network.PTDF(:, node_k) * p_i;
        predicted_flow = current_delta_flow + delta_flow_i;
        
        if all(predicted_flow <= margin_upper) && all(predicted_flow >= margin_lower)
            % --- 未违反约束：接受该设备 ---
            u_opt_t(idx) = true;
            P_total_t = P_total_t + p_i;
            cost_t = cost_t + total_cost_per_device(idx);
            current_delta_flow = predicted_flow; % 更新当前潮流
            
            % --- 检查是否满足总需求 ---
            if P_total_t >= P_req_t
                break; % 需求已满足，停止迭代
            end
        else
            % --- 违反约束：拒绝该设备 ---
            % u_opt_t(idx) 保持为 false
        end
    end
    
    is_feasible = (P_total_t >= P_req_t);
end