% hourly_ga_fitness_constrained_FIXED.m (V2 - 成本修正版)
function fitness = hourly_ga_fitness_constrained_FIXED(x, ...
    P_ac_hourly, P_ev_hourly, P_req_hourly, steps_per_hour, ...
    Location_AC, Location_EV, PTDF_matrix, P_Line_Base_hourly, P_Line_Max, N_bus, N_line, ...
    c_ac, c_ev) % <-- [!!! 关键修改 3a !!!] 接收新增的成本向量
    % *** (V2 修正): 此函数现在模拟MILP的 *成本最小化* 选择逻辑 ***

    n_ac = round(x(1));
    n_ev = round(x(2));
    
    large_penalty_factor = 1e9; % 功率不足的惩罚
    network_penalty_factor = 1e8; % 网络过载的惩罚

    % --- 1. 找到本小时的峰值需求 ---
    [max_demand_hourly, peak_idx_all] = max(P_req_hourly);
    peak_idx = peak_idx_all(1); % 只取第一个峰值点
    
    if max_demand_hourly <= 1e-6
        fitness = 0; % 如果本小时无需求，则成本为0
        return;
    end

    num_avail_ac = size(P_ac_hourly, 1);
    num_avail_ev = size(P_ev_hourly, 1);
    
    actual_n_ac = min(n_ac, num_avail_ac);
    actual_n_ev = min(n_ev, num_avail_ev);

    % --- [!!! 关键修改 3b !!!] ---
    % 我们不再按“功率”排序，而是按“成本”排序，以模拟MILP的行为。
    % MILP的目标函数是最小化 *小时总成本*。
    
    % --- 2. 检查约束 (在峰值时刻) ---
    % 检查 *成本最低* 的 n 台设备是否满足峰值需求
    
    % (2a. 计算用于排序的小时总成本)
    hourly_cost_ac = zeros(num_avail_ac, 1);
    if num_avail_ac > 0
        % (N_ac x T) .* (N_ac x 1) -> (N_ac x T)，然后按行求和
        hourly_cost_ac = sum(P_ac_hourly .* repmat(c_ac, 1, steps_per_hour), 2); 
    end
    
    hourly_cost_ev = zeros(num_avail_ev, 1);
    if num_avail_ev > 0
        hourly_cost_ev = sum(P_ev_hourly .* repmat(c_ev, 1, steps_per_hour), 2);
    end

    % (2b. 找出成本最低的 N 台设备的 *索引*)
    [~, sorted_cost_idx_ac] = sort(hourly_cost_ac, 'ascend');
    [~, sorted_cost_idx_ev] = sort(hourly_cost_ev, 'ascend');
    
    selected_ac_indices = [];
    if actual_n_ac > 0 && num_avail_ac > 0
        selected_ac_indices = sorted_cost_idx_ac(1:actual_n_ac);
    end
    
    selected_ev_indices = [];
    if actual_n_ev > 0 && num_avail_ev > 0
        selected_ev_indices = sorted_cost_idx_ev(1:actual_n_ev);
    end

    % (2c. 检查这组 *成本最低* 的设备在 *峰值时刻* 的功率)
    current_potential_power = 0;
    if ~isempty(selected_ac_indices)
        current_potential_power = current_potential_power + sum(P_ac_hourly(selected_ac_indices, peak_idx));
    end
    if ~isempty(selected_ev_indices)
        current_potential_power = current_potential_power + sum(P_ev_hourly(selected_ev_indices, peak_idx));
    end
    
    power_shortage = max(0, max_demand_hourly - current_potential_power);
    
    if power_shortage > 1e-3 
        fitness = power_shortage * large_penalty_factor; % 功率不足，返回大惩罚
        return;
    end
    % --- [!!! 修改结束 3b !!!] ---


    % --- 3. 估算网络潮流约束 ---
    % (此部分逻辑不变，但现在使用的是 *成本最低* 的设备集)
    network_penalty = 0;
    if (actual_n_ac > 0 || actual_n_ev > 0) && N_line > 0
        
        P_dispatch_per_node = zeros(N_bus, 1);
        
        % 找出这些被选中设备的 *峰值功率*
        P_dispatch_ac = [];
        if ~isempty(selected_ac_indices)
             P_dispatch_ac = P_ac_hourly(selected_ac_indices, peak_idx);
        end
        
        P_dispatch_ev = [];
        if ~isempty(selected_ev_indices)
            P_dispatch_ev = P_ev_hourly(selected_ev_indices, peak_idx);
        end

        % 找出这些设备的 *位置*
        Loc_dispatch_ac = [];
        if ~isempty(selected_ac_indices) && ~isempty(Location_AC)
            Loc_dispatch_ac = Location_AC(selected_ac_indices);
        end
        
        Loc_dispatch_ev = [];
         if ~isempty(selected_ev_indices) && ~isempty(Location_EV)
            Loc_dispatch_ev = Location_EV(selected_ev_indices);
         end
        
        % 使用 accumarray 计算每个节点的总注入功率
        if actual_n_ac > 0 && ~isempty(Loc_dispatch_ac)
            P_dispatch_per_node = accumarray(Loc_dispatch_ac, P_dispatch_ac, [N_bus 1]);
        end
        if actual_n_ev > 0 && ~isempty(Loc_dispatch_ev)
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

    % --- 4. 如果功率和网络约束都满足，计算 SDCI 和 Rho 指标 ---
    % (此部分逻辑不变，我们仍然使用 SDCI/Rho 作为最终的 *优化目标*)
    
    p_ac_agg = zeros(steps_per_hour, 1);
    p_ev_agg = zeros(steps_per_hour, 1);
    if ~isempty(selected_ac_indices)
        p_ac_agg = sum(P_ac_hourly(selected_ac_indices, :), 1)'; % (T_h x 1)
    end
     if ~isempty(selected_ev_indices)
        p_ev_agg = sum(P_ev_hourly(selected_ev_indices, :), 1)'; % (T_h x 1)
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