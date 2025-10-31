%% 局部函数: 上层GA的目标函数 (带功率约束惩罚)
function fitness = hourly_ga_fitness_constrained(x, P_ac_hourly, P_ev_hourly, P_req_hourly, steps_per_hour)
    n_ac = round(x(1));
    n_ev = round(x(2));
    
    % --- 1. 检查功率约束 ---
    power_penalty = 0;
    large_penalty_factor = 1e9; % 巨大的惩罚因子
    max_demand_hourly = max(P_req_hourly); % 本小时最大需求
    min_potential_power_hourly = inf; % 初始化为无穷大

    num_avail_ac = size(P_ac_hourly, 1);
    num_avail_ev = size(P_ev_hourly, 1);
    
    actual_n_ac = min(n_ac, num_avail_ac); % 防止索引越界
    actual_n_ev = min(n_ev, num_avail_ev);

    for t = 1:steps_per_hour
        current_potential_power = 0;
        if actual_n_ac > 0 && num_avail_ac > 0
            sorted_p_ac = sort(P_ac_hourly(:, t), 'descend');
            current_potential_power = current_potential_power + sum(sorted_p_ac(1:actual_n_ac));
        end
        if actual_n_ev > 0 && num_avail_ev > 0
            sorted_p_ev = sort(P_ev_hourly(:, t), 'descend');
            current_potential_power = current_potential_power + sum(sorted_p_ev(1:actual_n_ev));
        end
        min_potential_power_hourly = min(min_potential_power_hourly, current_potential_power);
    end

    %% --- MODIFICATION START ---
    % 为上层GA增加 "冗余因子" (安全系数)
    % 迫使GA选择比严格需求更多的设备，为下层MILP提供灵活性以规避PTDF约束
    safety_factor = 1.2; % 增加20%的功率冗余 (您可以调整此系数)
    
    power_shortage = max(0, (max_demand_hourly * safety_factor) - min_potential_power_hourly);
    %% --- MODIFICATION END ---
    
    if power_shortage > 1e-3 % 如果功率不足 (允许一点点误差)
        power_penalty = power_shortage * large_penalty_factor;
        fitness = power_penalty; % 直接返回大惩罚，无需计算后续指标
        return;
    end

    % --- 2. 如果功率约束满足，计算 SDCI 和 Rho 指标 ---
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
    
    n_ac_vec = ones(steps_per_hour, 1) * actual_n_ac; % 使用实际选取的数量
    n_ev_vec = ones(steps_per_hour, 1) * actual_n_ev;
    avg_p_ac = p_ac_agg ./ (actual_n_ac + 1e-9);
    avg_p_ev = p_ev_agg ./ (actual_n_ev + 1e-9);

    sdci_hourly = calculateSDCI(n_ac_vec, n_ev_vec, avg_p_ac, avg_p_ev);
    rho_hourly = calculateSpearmanRho(n_ac_vec, avg_p_ac, n_ev_vec, avg_p_ev);
    
    % 目标：最小化 SDCI 和 |Rho| (如果希望相关性尽量小，无论正负)
    w_sdci = 0.5;
    w_rho = 0.5;
    
    fitness = w_sdci * sdci_hourly + w_rho * abs(rho_hourly); % 或者 w_rho * (rho_hourly + 1) 如果希望负相关

    % 增加一个惩罚项，如果总功率为0 (虽然理论上功率约束已经处理了)
    if sum(p_ac_agg) + sum(p_ev_agg) < 1e-6
        fitness = fitness + 1e6;
    end
end