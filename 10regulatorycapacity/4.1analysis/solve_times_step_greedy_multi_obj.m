% solve_times_step_greedy_multi_obj.m
function [u_opt_t, P_total_t, cost_t, is_feasible] = ...
    solve_times_step_greedy_multi_obj(p_devices_t, c_devices, P_req_t, Network, ...
                                 W_cost, W_complementarity, Ratio_AC_t, device_types)
    % (v3) 使用集成了成本和互补性(SDCI+Rho)启发式的多目标贪心算法求解
    %
    % 新增输入:
    %   W_cost:             (scalar) 成本权重
    %   W_complementarity: (scalar) 互补性权重 (同时用于SDCI和Rho)
    %   Ratio_AC_t:         (scalar) 此刻 AC潜力 / (AC潜力 + EV潜力)
    %   device_types:       (N x 1) 1=AC, 2=EV

    N = length(p_devices_t);
    u_opt_t = false(N, 1);
    P_total_t = 0;
    cost_t = 0;
    
    if P_req_t <= 0
        is_feasible = true;
        return;
    end

    % --- 多目标贪心策略 ---
    
    % 1. 计算原始成本 (用于最终成本核算)
    total_cost_per_device = p_devices_t .* c_devices;

    % 2. 找出潜力 > 0 且成本有限的设备
    valid_indices = find(p_devices_t > 1e-3 & isfinite(total_cost_per_device));
    if isempty(valid_indices)
        is_feasible = false;
        return;
    end
    
    % 3. 计算多目标得分 (仅针对有效设备)
    
    % 3a. 成本得分 (归一化)
    costs_valid = total_cost_per_device(valid_indices);
    max_cost = max(costs_valid);
    if max_cost <= 1e-6, max_cost = 1; end % 避免除零
    norm_cost_score = costs_valid / max_cost;
    
    % 3b. 互补性得分 (Complementarity Score)
    % 策略: 惩罚 (高分) 使用 "稀缺" 资源，鼓励 (低分) 使用 "丰富" 资源
    % 这将自动促进低 SDCI (min(P_AC, P_EV)) 和低 Rho (负相关)
    
    Ratio_EV_t = 1.0 - Ratio_AC_t;
    
    % 获取有效设备的类型
    types_valid = device_types(valid_indices);
    
    % 初始化得分
    complementarity_score = zeros(length(valid_indices), 1);
    
    % 如果是AC (类型1)，得分 = 稀缺资源的比例 (EV)
    % (如果EV稀缺(Ratio_EV_t小), AC得分就低(好), 优先用AC)
    complementarity_score(types_valid == 1) = Ratio_EV_t;
    
    % 如果是EV (类型2)，得分 = 稀缺资源的比例 (AC)
    % (如果AC稀缺(Ratio_AC_t小), EV得分就低(好), 优先用EV)
    complementarity_score(types_valid == 2) = Ratio_AC_t;
    
    % 3c. 计算综合得分
    combined_metric = (W_cost * norm_cost_score) + (W_complementarity * complementarity_score);

    % 4. 按“综合得分”升序排序 (优先选择得分最低的)
    [~, sorted_order] = sort(combined_metric, 'ascend');
    sorted_indices = valid_indices(sorted_order);
    
    % --- 网络约束准备 (同前) ---
    L = size(Network.PTDF, 1); 
    current_delta_flow = zeros(L, 1); 
    margin_upper = Network.LineLimits - Network.BaseFlow;
    margin_lower = -Network.LineLimits - Network.BaseFlow;

    % 5. 迭代选择设备 (同前)
    for idx = sorted_indices' 
        
        p_i = p_devices_t(idx);
        node_k = Network.DeviceNodes(idx);
        
        delta_flow_i = Network.PTDF(:, node_k) * p_i;
        predicted_flow = current_delta_flow + delta_flow_i;
        
        if all(predicted_flow <= margin_upper) && all(predicted_flow >= margin_lower)
            % --- 未违反约束：接受该设备 ---
            u_opt_t(idx) = true;
            P_total_t = P_total_t + p_i;
            % 关键：成本核算仍使用 *原始* 成本
            cost_t = cost_t + total_cost_per_device(idx); 
            current_delta_flow = predicted_flow; 
            
            % --- 检查是否满足总需求 ---
            if P_total_t >= P_req_t
                break; 
            end
        else
            % --- 违反约束：拒绝该设备 ---
        end
    end
    
    is_feasible = (P_total_t >= P_req_t);
end