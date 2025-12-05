function [AC_Up_Sum, AC_Down_Sum] = run_AC_simulation_MC(random_seed, acFileName)
    % run_AC_simulation_MC
    % 修改：强制系统按基线运行，不响应外部指令，以评估最大物理潜力。
    % [修改] 仿真时间调整为 6:00 到 30:00 (24小时，跨天)
    
    %% 1. 系统初始化
    if nargin > 0
        rng(random_seed, 'Threefry'); % 使用传入的种子控制随机性
    end
    
    % --- [修改] 时间参数调整 ---
    simulation_start_hour = 6;  % 仿真开始于 6:00
    simulation_end_hour   = 30; % 仿真结束于次日 6:00
    dt = 5/60;                  % 时间步长 (5分钟)
    time_points = simulation_start_hour:dt:simulation_end_hour; 
    T_steps_total = length(time_points);
    % -------------------------

    base_price = 30; 

    %% 2. 初始化 AC 参数
    acFile = acFileName;
    try
        ACs = initializeACsFromExcel(acFile);
    catch ME
        error('初始化 AC 失败: %s', ME.message);
    end
    num_AC = length(ACs);

    for i = 1:num_AC
        ACs(i).Tset_original = ACs(i).Tset;
        ACs(i).Tmax_original = ACs(i).Tmax;
        ACs(i).Tmin_original = ACs(i).Tmin;
        if ~isfield(ACs(i), 'p_incentive')
            ACs(i).p_incentive = round(50*rand(), 1);
        end
    end

    %% 3. 激励响应与预计算
    p_min = 15; p_max = 50; p_min_prime = 10; p_max_prime = 40; T_set_max = 3;
    current_p = 25.0;

    temp_ACs = ACs;
    max_Tset_all = max([ACs.Tset_original]);
    
    % 环境温度生成 (time_points 现在是 6..30, 15.0 对应下午3点峰值)
    T_ja_min_ambient = max_Tset_all + 0.1; 
    T_ja_peak_ambient = max_Tset_all + 3.0; 
    T_ja_mean = (T_ja_min_ambient + T_ja_peak_ambient) / 2;
    T_ja_amplitude = (T_ja_peak_ambient - T_ja_min_ambient) / 2;
    base_trend = T_ja_mean + T_ja_amplitude * cos(2*pi*(time_points - 15)/24); 
    
    window_size = 2 * round(1/dt); 
    noise_padding = ceil(window_size / 2);
    white_noise = randn(1, T_steps_total + 2 * noise_padding);
    fluctuations_raw = movmean(white_noise, window_size);
    fluctuations_centered = fluctuations_raw(noise_padding + 1 : noise_padding + T_steps_total);
    fluctuation_scale = T_ja_amplitude * 0.2; 
    base_ambient_temp_unified = max(base_trend + (fluctuations_centered / std(fluctuations_centered, 'omitnan')) * fluctuation_scale, T_ja_min_ambient); 

    parfor i = 1:num_AC
        participation = calculateParticipation(current_p, base_price);
        [~, ~, deltaT] = incentiveTempAC(current_p, p_min, p_max, p_min_prime, p_max_prime, T_set_max);
        temp_ACs(i).ptcp = (rand() < participation);

        if temp_ACs(i).ptcp
            temp_ACs(i).Tmax = temp_ACs(i).Tset_original + deltaT;
            temp_ACs(i).Tmin = temp_ACs(i).Tset_original - deltaT;
            temp_ACs(i).Tset = temp_ACs(i).Tmin + (temp_ACs(i).Tmax - temp_ACs(i).Tmin) * rand();
        end
        temp_ACs(i).T_ja = base_ambient_temp_unified;

        [alpha, beta, gamma] = calculateACABC_single(...
            temp_ACs(i).R, temp_ACs(i).C, temp_ACs(i).eta,...
            temp_ACs(i).Tmax, temp_ACs(i).Tmin, temp_ACs(i).Tset, dt);
        temp_ACs(i).alpha = alpha;
        temp_ACs(i).beta = beta;
        temp_ACs(i).gamma = gamma;
    end
    ACs = temp_ACs;

    ACs_participating = ACs([ACs.ptcp]);
    num_AC_participating = length(ACs_participating);
    
    if num_AC_participating == 0
        AC_Up_Sum = zeros(T_steps_total, 1); AC_Down_Sum = zeros(T_steps_total, 1); 
        return;
    end
    
    AggParams = calculateAggregatedACParams(ACs_participating);
    CURRENT_SOC_AC = [ACs_participating.SOC]'; 
    
    Agg_P_Potential_Up_History = zeros(T_steps_total, 1);
    Agg_P_Potential_Down_History = zeros(T_steps_total, 1);

    %% 5. 主时间循环 (基线模式)
    for t_idx = 1:T_steps_total
        SOC_agg_t = mean(CURRENT_SOC_AC, 'omitnan');

        % 【核心修正】: 强制调节指令为 0 (Delta_P = 0)
        % 模拟系统在无外部干预下的自然运行状态，以此评估最大调节潜力
        Delta_P_S_command = 0; 

        % 预测下一时刻的自然 SOC 演变
        SOC_target_next = AggParams.A * SOC_agg_t + AggParams.B * Delta_P_S_command + AggParams.C;
        SOC_target_next = max(0, min(1, SOC_target_next));

        temp_AC_Up_agg = 0;
        temp_AC_Down_agg = 0;
        temp_SOC_next = zeros(num_AC_participating, 1);

        parfor i = 1:num_AC_participating
            ac_i = ACs_participating(i);
            soc_curr = CURRENT_SOC_AC(i);
            
            T_ja_val = base_ambient_temp_unified(t_idx);
            P_base_i = ACbaseP_single(T_ja_val, ac_i.Tset, ac_i.R, ac_i.eta);
            
            % 计算该状态下的最大物理调节潜力
            [P_plus, P_minus] = calculateACAdjustmentPotentia(...
                P_base_i, 1.1*max(abs(ACbaseP_single(ac_i.T_ja, ac_i.Tset, ac_i.R, ac_i.eta))), 0, ...
                ac_i.alpha, ac_i.beta, ac_i.gamma, soc_curr, dt);

            temp_AC_Up_agg = temp_AC_Up_agg + P_plus;
            temp_AC_Down_agg = temp_AC_Down_agg + P_minus;

            % 状态更新：跟随基线 (Delta_Pj = 0)
            % 只有这样，系统才能始终保持在拥有最大潜力的“舒适区中心”
            soc_next = updateACSOC_single(soc_curr, 0, ac_i.alpha, ac_i.beta, ac_i.gamma);
            temp_SOC_next(i) = soc_next;
        end

        Agg_P_Potential_Up_History(t_idx) = temp_AC_Up_agg;
        Agg_P_Potential_Down_History(t_idx) = temp_AC_Down_agg;
        CURRENT_SOC_AC = temp_SOC_next;
    end
    
    AC_Up_Sum = Agg_P_Potential_Up_History;
    AC_Down_Sum = Agg_P_Potential_Down_History;
end