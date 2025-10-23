function fitness = objectiveFunction_down_real_time(x, n_AC_max, n_EV_max, deltaP_AC, deltaP_EV, T, steps_per_hour, theta_2, theta_3, lambda_1, lambda_2)
    % x: 决策变量，长度 2T, 对应每个细分时间步的设备数量
    % deltaP_AC, deltaP_EV: 个体设备在每个细分时间步的调节能力 (N_devices x T)
    % T: 总的细分时间步数量
    % steps_per_hour: 每小时包含的细分时间步数量 (例如 1/0.05 = 20)

    n_AC_fine_steps = round(x(1:T));
    n_EV_fine_steps = round(x(T+1:2*T));

    p_AC = zeros(1, T);
    p_EV = zeros(1, T);

    num_hours = T / steps_per_hour;
    if mod(T, steps_per_hour) ~= 0
        error('T (%d) 必须是 steps_per_hour (%d) 的整数倍', T, steps_per_hour);
    end

    % --- 并行化外层的小时循环 ---
    % 由于每个小时的设备选择和后续功率计算是独立的，适合parfor
    % 收集每个小时的功率段
    p_AC_hourly_segments = cell(1, num_hours);
    p_EV_hourly_segments = cell(1, num_hours);

    parfor h = 1:num_hours % 外层循环按小时并行化
        start_step_of_hour = (h-1) * steps_per_hour + 1;
        end_step_of_hour = h * steps_per_hour;
        
        current_segment_len = end_step_of_hour - start_step_of_hour + 1;
        local_p_AC_segment = zeros(1, current_segment_len);
        local_p_EV_segment = zeros(1, current_segment_len);

        % 1. 在该小时的开始 (start_step_of_hour)，根据该时刻的计划设备参与数量
        %    和该时刻的设备能力进行排序，并选定参与设备集合
        
        num_ac_to_select_this_hour = max(0, min(n_AC_max, n_AC_fine_steps(start_step_of_hour)));
        num_ev_to_select_this_hour = max(0, min(n_EV_max, n_EV_fine_steps(start_step_of_hour)));

        selected_ac_indices_for_this_hour = [];
        if size(deltaP_AC,1) > 0 && num_ac_to_select_this_hour > 0
            % 使用小时开始时的能力进行排序和选择设备索引
            [~, sorted_device_indices_ac] = sort(deltaP_AC(:, start_step_of_hour), 'descend');
            actual_num_selected_ac = min(num_ac_to_select_this_hour, size(deltaP_AC,1));
            selected_ac_indices_for_this_hour = sorted_device_indices_ac(1:actual_num_selected_ac);
        end
        
        selected_ev_indices_for_this_hour = [];
        if size(deltaP_EV,1) > 0 && num_ev_to_select_this_hour > 0
            [~, sorted_device_indices_ev] = sort(deltaP_EV(:, start_step_of_hour), 'descend');
            actual_num_selected_ev = min(num_ev_to_select_this_hour, size(deltaP_EV,1));
            selected_ev_indices_for_this_hour = sorted_device_indices_ev(1:actual_num_selected_ev);
        end

        % 2. 在该小时内的每个细分时间步，使用这个小时初选定的设备集合，
        %    并用这些设备在当前细分时间步的实际能力来计算总调节功率
        segment_idx = 0;
        for t_fine = start_step_of_hour:end_step_of_hour
            segment_idx = segment_idx + 1;
            if ~isempty(selected_ac_indices_for_this_hour)
                local_p_AC_segment(segment_idx) = sum(deltaP_AC(selected_ac_indices_for_this_hour, t_fine));
            else
                local_p_AC_segment(segment_idx) = 0;
            end
            
            if ~isempty(selected_ev_indices_for_this_hour)
                local_p_EV_segment(segment_idx) = sum(deltaP_EV(selected_ev_indices_for_this_hour, t_fine));
            else
                local_p_EV_segment(segment_idx) = 0;
            end
        end
        p_AC_hourly_segments{h} = local_p_AC_segment;
        p_EV_hourly_segments{h} = local_p_EV_segment;
    end

    % 拼接并行计算的结果
    p_AC = horzcat(p_AC_hourly_segments{:});
    p_EV = horzcat(p_EV_hourly_segments{:});


    %% ===== 核心指标计算逻辑 (保持不变) =====
    total_power = sum(p_AC) + sum(p_EV);
    
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
    if any(n_AC_fine_steps > n_AC_max) || any(n_EV_fine_steps > n_EV_max)
        penalty = penalty + 1e6;
    end
    if any(n_AC_fine_steps < 0) || any(n_EV_fine_steps < 0) 
        penalty = penalty + 1e6;
    end
    
    fitness = base_fitness - penalty;
end