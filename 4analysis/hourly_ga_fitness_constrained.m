% hourly_ga_fitness_constrained.m (已按 Task 2 修改)
function fitness = hourly_ga_fitness_constrained(x, ...
    P_ac_hourly, P_ev_hourly, P_req_hourly, steps_per_hour, ...
    Location_AC, Location_EV, PTDF_matrix, P_Line_Base_hourly, P_Line_Max, N_bus, N_line)
    % *** 【Task 2 修正】: 增加了新的输入参数 (Location_AC, ..., N_line) ***

    n_ac = round(x(1));
    n_ev = round(x(2));
    
    large_penalty_factor = 1e9; % 功率不足的惩罚
    network_penalty_factor = 1e8; % 网络过载的惩罚

    % --- 1. 找到本小时的峰值需求 ---
    [max_demand_hourly, peak_idx] = max(P_req_hourly);
    if max_demand_hourly <= 1e-6
        fitness = 0; % 如果本小时无需求，则成本为0
        return;
    end

    num_avail_ac = size(P_ac_hourly, 1);
    num_avail_ev = size(P_ev_hourly, 1);
    
    actual_n_ac = min(n_ac, num_avail_ac);
    actual_n_ev = min(n_ev, num_avail_ev);

    % --- 2. 检查功率约束 (在峰值时刻) ---
    % 检查 *最优* 的 n 台设备是否满足峰值需求
    current_potential_power = 0;
    
    if actual_n_ac > 0 && num_avail_ac > 0
        sorted_p_ac = sort(P_ac_hourly(:, peak_idx), 'descend');
        current_potential_power = current_potential_power + sum(sorted_p_ac(1:actual_n_ac));
    end
    if actual_n_ev > 0 && num_avail_ev > 0
        sorted_p_ev = sort(P_ev_hourly(:, peak_idx), 'descend');
        current_potential_power = current_potential_power + sum(sorted_p_ev(1:actual_n_ev));
    end
    
    power_shortage = max(0, max_demand_hourly - current_potential_power);
    
    if power_shortage > 1e-3 
        fitness = power_shortage * large_penalty_factor; % 功率不足，返回大惩罚
        return;
    end

    % --- 3. ***【Task 2 修正】: 估算网络潮流约束 *** ---
    % 如果功率满足，我们估算这组最优设备在峰值时刻是否导致网络过载
    network_penalty = 0;
    if (actual_n_ac > 0 || actual_n_ev > 0) && N_line > 0
        
        P_dispatch_per_node = zeros(N_bus, 1);
        
        % 获取峰值时刻的功率和排序索引
        [sorted_p_ac, idx_ac] = sort(P_ac_hourly(:, peak_idx), 'descend');
        [sorted_p_ev, idx_ev] = sort(P_ev_hourly(:, peak_idx), 'descend');
        
        % 找出被选中的设备 *索引*
        selected_ac_indices = idx_ac(1:actual_n_ac);
        selected_ev_indices = idx_ev(1:actual_n_ev);
        
        % 找出这些设备的 *功率*
        P_dispatch_ac = sorted_p_ac(1:actual_n_ac);
        P_dispatch_ev = sorted_p_ev(1:actual_n_ev);

        % 找出这些设备的 *位置*
        Loc_dispatch_ac = Location_AC(selected_ac_indices);
        Loc_dispatch_ev = Location_EV(selected_ev_indices);
        
        % 使用 accumarray 计算每个节点的总注入功率
        if actual_n_ac > 0
            P_dispatch_per_node = accumarray(Loc_dispatch_ac, P_dispatch_ac, [N_bus 1]);
        end
        if actual_n_ev > 0
            P_dispatch_per_node = P_dispatch_per_node + accumarray(Loc_dispatch_ev, P_dispatch_ev, [N_bus 1]);
        end

        % 估算潮流变化
        Delta_P_line = PTDF_matrix * P_dispatch_per_node;
        P_base_t = P_Line_Base_hourly(:, peak_idx);
        P_line_final = P_base_t + Delta_P_line;
        
        % 计算过载惩罚
        line_overload = max(0, abs(P_line_final) - P_Line_Max);
        network_penalty = sum(line_overload) * network_penalty_factor;
        
        if network_penalty > 0
            fitness = network_penalty; % 网络过载，返回惩罚
            return;
        end
    end
    % ************************************************

    % --- 4. 如果功率和网络约束都满足，计算 SDCI 和 Rho 指标 ---
    % (此处的计算使用简化的 "top N" 功率，作为指标的代理)
    p_ac_agg = zeros(steps_per_hour, 1);
    p_ev_agg = zeros(steps_per_hour, 1);
    if actual_n_ac > 0 && num_avail_ac > 0
        for t = 1:steps_per_hour
            sorted_p_ac = sort(P_ac_hourly(:, t), 'descend');
            p_ac_agg(t) = sum(sorted_p_ac(1:actual_n_ac));
        end
    end
     if actual_n_ev > 0 && num_avail_ev > 0
        for t = 1:steps_per_hour
            sorted_p_ev = sort(P_ev_hourly(:, t), 'descend');
            p_ev_agg(t) = sum(sorted_p_ev(1:actual_n_ev));
        end
    end
    
    n_ac_vec = ones(steps_per_hour, 1) * actual_n_ac;
    n_ev_vec = ones(steps_per_hour, 1) * actual_n_ev;
    avg_p_ac = p_ac_agg ./ (actual_n_ac + 1e-9);
    avg_p_ev = p_ev_agg ./ (actual_n_ev + 1e-9);

    sdci_hourly = calculateSDCI(n_ac_vec, n_ev_vec, avg_p_ac, avg_p_ev);
    rho_hourly = calculateSpearmanRho(n_ac_vec, avg_p_ac, n_ev_vec, avg_p_ev);
    
    % 目标：最小化 SDCI 和 |Rho| 
    w_sdci = 0.5;
    w_rho = 0.5;
    
    fitness = w_sdci * sdci_hourly + w_rho * abs(rho_hourly); 
end