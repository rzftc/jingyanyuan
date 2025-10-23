function fitness = objectiveFunction_down_hourly(x, n_AC_max, n_EV_max, deltaP_AC, deltaP_EV, T, steps_per_hour, theta_2, theta_3, lambda_1, lambda_2)
    % x: 决策变量，长度 2T, 对应每个细分时间步的设备数量
    % deltaP_AC, deltaP_EV: 个体设备在每个细分时间步的调节能力 (N_devices x T)
    % T: 总的细分时间步数量
    % steps_per_hour: 每小时包含的细分时间步数量 (例如 1/0.05 = 20)

    % 分解决策变量 (每个细分时间步的设备数量)
    n_AC_fine_steps = round(x(1:T));
    n_EV_fine_steps = round(x(T+1:2*T));

    % 初始化每个细分时间步的总调节能力
    p_AC = zeros(1, T);
    p_EV = zeros(1, T);

    num_hours = T / steps_per_hour;
    if mod(T, steps_per_hour) ~= 0
        error('T (%d) 必须是 steps_per_hour (%d) 的整数倍', T, steps_per_hour);
    end

    for h = 1:num_hours % 外层循环按小时
        % 当前小时的第一个细分时间步的索引
        start_step_of_hour = (h-1) * steps_per_hour + 1;
        % 当前小时的最后一个细分时间步的索引
        end_step_of_hour = h * steps_per_hour;

        % 1. 在该小时的开始 (start_step_of_hour)，根据该时刻的设备参与数量
        %    和该时刻的设备能力进行排序，并选定参与设备集合
        
        % 获取当前小时第一个时间步计划参与的设备数量
        % 同时确保不超过最大可调度数，且不小于0
        num_ac_this_hour = max(0, min(n_AC_max, n_AC_fine_steps(start_step_of_hour)));
        num_ev_this_hour = max(0, min(n_EV_max, n_EV_fine_steps(start_step_of_hour)));

        % 空调设备选择 (基于小时初的能力排序)
        selected_ac_indices_for_hour = []; % 这个小时内参与的空调的索引
        if size(deltaP_AC,1) > 0 && num_ac_this_hour > 0
            % 使用聚合周期开始时的能力进行排序和选择
            [~, sorted_indices_ac] = sort(deltaP_AC(:, start_step_of_hour), 'descend');
            num_to_select_ac = min(num_ac_this_hour, size(deltaP_AC,1)); % 实际能选的数量
            selected_ac_indices_for_hour = sorted_indices_ac(1:num_to_select_ac);
        end
        
        % 电动汽车设备选择 (基于小时初的能力排序)
        selected_ev_indices_for_hour = []; % 这个小时内参与的EV的索引
        if size(deltaP_EV,1) > 0 && num_ev_this_hour > 0
             % 使用聚合周期开始时的能力进行排序和选择
            [~, sorted_indices_ev] = sort(deltaP_EV(:, start_step_of_hour), 'descend');
            num_to_select_ev = min(num_ev_this_hour, size(deltaP_EV,1)); % 实际能选的数量
            selected_ev_indices_for_hour = sorted_indices_ev(1:num_to_select_ev);
        end

        % 2. 在该小时内的每个细分时间步，使用这个小时初选定的设备集合，
        %    并用这些设备在当前细分时间步的实际能力来计算总调节功率
        for t_fine = start_step_of_hour:end_step_of_hour
            % 空调总调节功率
            if ~isempty(selected_ac_indices_for_hour)
                p_AC(t_fine) = sum(deltaP_AC(selected_ac_indices_for_hour, t_fine));
            else
                p_AC(t_fine) = 0;
            end
            
            % 电动汽车总调节功率
            if ~isempty(selected_ev_indices_for_hour)
                p_EV(t_fine) = sum(deltaP_EV(selected_ev_indices_for_hour, t_fine));
            else
                p_EV(t_fine) = 0;
            end
        end
    end

    %% ===== 核心指标计算逻辑 (保持不变) =====
    total_power = sum(p_AC) + sum(p_EV);
    
    % 计算最大潜力时，仍然按每个细分时间步考虑调度所有最大数量设备
    max_power_AC_fine_steps = zeros(1,T);
    max_power_EV_fine_steps = zeros(1,T);

    if size(deltaP_AC,1) > 0 && n_AC_max > 0
        for t_fine = 1:T
            sorted_AC_all = sort(deltaP_AC(:,t_fine), 'descend');
            max_power_AC_fine_steps(t_fine) = sum(sorted_AC_all(1:min(n_AC_max, size(deltaP_AC,1))));
        end
    end
    if size(deltaP_EV,1) > 0 && n_EV_max > 0
        for t_fine = 1:T
            sorted_EV_all = sort(deltaP_EV(:,t_fine), 'descend');
            max_power_EV_fine_steps(t_fine) = sum(sorted_EV_all(1:min(n_EV_max, size(deltaP_EV,1))));
        end
    end
    max_total_power = sum(max_power_AC_fine_steps) + sum(max_power_EV_fine_steps);
        
    normalized_reward = total_power / (max_total_power + eps);
    
    avg_deltaP_AC_fine_steps = p_AC ./ (n_AC_fine_steps + eps);
    avg_deltaP_EV_fine_steps = p_EV ./ (n_EV_fine_steps + eps);
    
    avg_deltaP_AC_fine_steps(isnan(avg_deltaP_AC_fine_steps) | isinf(avg_deltaP_AC_fine_steps)) = 0;
    avg_deltaP_EV_fine_steps(isnan(avg_deltaP_EV_fine_steps) | isinf(avg_deltaP_EV_fine_steps)) = 0;

    SDCI_minus = calculateSDCI(n_AC_fine_steps', n_EV_fine_steps', avg_deltaP_AC_fine_steps', avg_deltaP_EV_fine_steps');
    rho_minus = calculateSpearmanRho(n_AC_fine_steps', avg_deltaP_AC_fine_steps', n_EV_fine_steps', avg_deltaP_EV_fine_steps');
    
    %% ===== 目标函数构建逻辑 (保持不变) =====
    base_fitness = normalized_reward - lambda_1*SDCI_minus - lambda_2*rho_minus;
    
    %% ===== 约束处理逻辑 (保持不变) =====
    penalty = 0;
    if SDCI_minus > theta_2
        penalty = penalty + 1e6*(SDCI_minus - theta_2);
    end
    if rho_minus > theta_3 
        penalty = penalty + 1e6*(rho_minus - theta_3);
    end
    % 设备数量约束仍然是针对每个细分时间步的原始决策变量
    if any(n_AC_fine_steps > n_AC_max) || any(n_EV_fine_steps > n_EV_max)
        penalty = penalty + 1e6;
    end
    if any(n_AC_fine_steps < 0) || any(n_EV_fine_steps < 0) % 确保数量非负
        penalty = penalty + 1e6;
    end
    
    fitness = base_fitness - penalty;
end