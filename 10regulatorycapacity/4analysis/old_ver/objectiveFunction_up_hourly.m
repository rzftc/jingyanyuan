function fitness = objectiveFunction_up_hourly(x, n_AC_max, n_EV_max, deltaP_AC, deltaP_EV, T, steps_per_hour, theta_1, theta_3, lambda_1, lambda_2)
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
        start_step_of_hour = (h-1) * steps_per_hour + 1;
        end_step_of_hour = h * steps_per_hour;

        num_ac_for_this_hour = max(0, min(n_AC_max, n_AC_fine_steps(start_step_of_hour)));
        num_ev_for_this_hour = max(0, min(n_EV_max, n_EV_fine_steps(start_step_of_hour)));

        selected_ac_indices_for_hour = [];
        if size(deltaP_AC,1) > 0 && num_ac_for_this_hour > 0
            [~, sorted_indices_ac] = sort(deltaP_AC(:, start_step_of_hour), 'descend');
            num_to_select_ac = min(num_ac_for_this_hour, size(deltaP_AC,1));
            selected_ac_indices_for_hour = sorted_indices_ac(1:num_to_select_ac);
        end
        
        selected_ev_indices_for_hour = [];
        if size(deltaP_EV,1) > 0 && num_ev_for_this_hour > 0
            [~, sorted_indices_ev] = sort(deltaP_EV(:, start_step_of_hour), 'descend');
            num_to_select_ev = min(num_ev_for_this_hour, size(deltaP_EV,1));
            selected_ev_indices_for_hour = sorted_indices_ev(1:num_to_select_ev);
        end

        for t_fine = start_step_of_hour:end_step_of_hour
            if ~isempty(selected_ac_indices_for_hour)
                p_AC(t_fine) = sum(deltaP_AC(selected_ac_indices_for_hour, t_fine));
            else
                p_AC(t_fine) = 0;
            end
            
            if ~isempty(selected_ev_indices_for_hour)
                p_EV(t_fine) = sum(deltaP_EV(selected_ev_indices_for_hour, t_fine));
            else
                p_EV(t_fine) = 0;
            end
        end
    end

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

    SDCI_plus = calculateSDCI(n_AC_fine_steps', n_EV_fine_steps', avg_deltaP_AC_fine_steps', avg_deltaP_EV_fine_steps');
    rho_plus = calculateSpearmanRho(n_AC_fine_steps', avg_deltaP_AC_fine_steps', n_EV_fine_steps', avg_deltaP_EV_fine_steps');
    
    %% ===== 目标函数构建逻辑 (保持不变) =====
    base_fitness = normalized_reward - lambda_1*SDCI_plus - lambda_2*rho_plus;
    
    %% ===== 约束处理逻辑 (保持不变) =====
    penalty = 0;
    if SDCI_plus > theta_1 % 使用上调的 theta_1
        penalty = penalty + 1e6*(SDCI_plus - theta_1);
    end
    if rho_plus > theta_3 
        penalty = penalty + 1e6*(rho_plus - theta_3);
    end
    if any(n_AC_fine_steps > n_AC_max) || any(n_EV_fine_steps > n_EV_max)
        penalty = penalty + 1e6;
    end
     if any(n_AC_fine_steps < 0) || any(n_EV_fine_steps < 0) 
        penalty = penalty + 1e6;
    end
    
    fitness = base_fitness - penalty;
end