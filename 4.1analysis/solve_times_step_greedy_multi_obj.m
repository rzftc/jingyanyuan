function [u_opt_t, P_total_t, cost_t, is_feasible] = ...
    solve_times_step_greedy_multi_obj(p_devices_t, c_devices, P_req_t, Network, ...
                                 W_cost, W_complementarity, Ratio_AC_t, device_types)
    % (v7 - 双向兼容贪心算法)
    N = length(p_devices_t);
    u_opt_t = false(N, 1); P_total_t = 0; cost_t = 0;
    
    % 0. 快速检查：需求为0，或需求方向与所有设备能力相反
    if abs(P_req_t) < 1e-6
         is_feasible = true; return;
    elseif (P_req_t > 0 && all(p_devices_t <= 1e-6)) || (P_req_t < 0 && all(p_devices_t >= -1e-6))
         is_feasible = false; return; 
    end

    % 1. 找出有效设备 (能力与需求同向，且幅度足够)
    if P_req_t > 0
        valid_idx = find(p_devices_t > 1e-4);
    else % 下调需求为负，找负值能力的设备
        valid_idx = find(p_devices_t < -1e-4);
    end
    
    if isempty(valid_idx)
        is_feasible = false; return;
    end

    % 2. 多目标评分
    % 2a. [修改] 成本得分 (必须使用功率绝对值来计算正的成本)
    costs_valid = abs(p_devices_t(valid_idx)) .* c_devices(valid_idx);
    max_c = max(costs_valid); if max_c < 1e-9, max_c = 1; end
    score_cost = costs_valid / max_c;
    
    % 2b. 互补性得分 (基于资源稀缺性)
    Ratio_EV_t = 1.0 - Ratio_AC_t;
    score_comp = zeros(length(valid_idx), 1);
    types_valid = device_types(valid_idx);
    score_comp(types_valid == 1) = Ratio_EV_t; % AC得分
    score_comp(types_valid == 2) = Ratio_AC_t; % EV得分
    
    % 2c. 综合得分 (越低越优先)
    [~, sort_order] = sort(W_cost * score_cost + W_complementarity * score_comp, 'ascend');
    sorted_idx = valid_idx(sort_order);
    
    % 3. 贪心选择 (带网络约束)
    margin_up = Network.LineLimits - Network.BaseFlow;
    margin_lo = -Network.LineLimits - Network.BaseFlow;
    curr_flow = zeros(size(Network.PTDF,1), 1);

    for i = 1:length(sorted_idx)
        idx = sorted_idx(i);
        % 检查网络约束
        flow_new = curr_flow + Network.PTDF(:, Network.DeviceNodes(idx)) * p_devices_t(idx);
        if all(flow_new <= margin_up + 1e-5) && all(flow_new >= margin_lo - 1e-5)
            u_opt_t(idx) = true;
            P_total_t = P_total_t + p_devices_t(idx);
            % [修改] 累加绝对值成本
            cost_t = cost_t + abs(p_devices_t(idx)) * c_devices(idx); 
            curr_flow = flow_new;
            
            % [修改] 检查是否满足需求 (考虑方向: 上调需>=, 下调需<=)
            if (P_req_t > 0 && P_total_t >= P_req_t) || (P_req_t < 0 && P_total_t <= P_req_t)
                break; 
            end
        end
    end
    
    % 4. 最终可行性 (考虑方向，允许微小误差)
    if P_req_t > 0
        is_feasible = (P_total_t >= P_req_t - 1e-3);
    else
        is_feasible = (P_total_t <= P_req_t + 1e-3);
    end
end