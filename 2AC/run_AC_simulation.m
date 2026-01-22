function [AC_Up_Sum, AC_Down_Sum, AC_Power_Sum] = run_AC_simulation(P_grid_command_24)
    % run_AC_simulation
    % 输入:
    %   P_grid_command_24: 24x1 的小时级电网**绝对功率**指令 (kW)
    %                      (注意：此处修改为接收绝对功率，而非调节偏差)
    % 输出:
    %   AC_Up_Sum:    (T_steps_total x 1) 所有单体上调潜力之和
    %   AC_Down_Sum:  (T_steps_total x 1) 所有单体下调潜力之和
    %   AC_Power_Sum: (T_steps_total x 1) 所有单体实际总功率之和

    %% 1. 系统初始化
    rng(2023, 'Threefry'); 
    T_total = 24; 
    dt = 5/60;    % 5分钟步长
    time_points = 0:dt:T_total; 
    T_steps_total = length(time_points);
    steps_per_hour = round(1/dt);
    
    % === 处理电网指令维度 (24 -> T_steps_total) ===
    % 检查输入是否为空或维度不对
    if length(P_grid_command_24) ~= 24
         warning('输入指令长度为 %d，期望为 24。正在尝试自动适配...', length(P_grid_command_24));
    end
    
    % 将小时级指令扩展到分钟级
    P_grid_command_series = repelem(P_grid_command_24, steps_per_hour);
    
    % 维度对齐
    current_len = length(P_grid_command_series);
    if current_len < T_steps_total
        padding = repmat(P_grid_command_series(end), T_steps_total - current_len, 1);
        P_grid_command_series = [P_grid_command_series(:); padding];
    elseif current_len > T_steps_total
        P_grid_command_series = P_grid_command_series(1:T_steps_total);
    end
    % ========================================================

    base_price = 30; 

    %% 2. 初始化 AC 参数
    acFile = 'AC_template2.xlsx';
    try
        ACs = initializeACsFromExcel(acFile);
    catch ME
        error('初始化 AC 失败: %s', ME.message);
    end
    num_AC = length(ACs);

    % 备份
    for i = 1:num_AC
        ACs(i).Tset_original = ACs(i).Tset;
        ACs(i).Tmax_original = ACs(i).Tmax;
        ACs(i).Tmin_original = ACs(i).Tmin;
        if ~isfield(ACs(i), 'p_incentive')
            ACs(i).p_incentive = round(60*rand(), 1);
        end
    end

    %% 3. 激励响应参数
    p_min = 15; p_max = 50; p_min_prime = 10; p_max_prime = 40; T_set_max = 3;
    current_p = 25.0;

    %% 4. 预计算
    temp_ACs = ACs;
    max_Tset_all = max([ACs.Tset_original]);
    
    % T_ja 生成逻辑
    T_ja_min_ambient = max_Tset_all + 0.1; 
    T_ja_peak_ambient = max_Tset_all + 6.0; 
    T_ja_mean = (T_ja_min_ambient + T_ja_peak_ambient) / 2;
    T_ja_amplitude = (T_ja_peak_ambient - T_ja_min_ambient) / 2;
    base_trend = T_ja_mean + T_ja_amplitude * cos(2*pi*(time_points - 15)/24); 
    
    % 添加噪声
    window_size = 2 * steps_per_hour; 
    noise_padding = ceil(window_size / 2);
    white_noise = randn(1, T_steps_total + 2 * noise_padding);
    fluctuations_raw = movmean(white_noise, window_size);
    fluctuations_centered = fluctuations_raw(noise_padding + 1 : noise_padding + T_steps_total);
    fluctuation_scale = T_ja_amplitude * 0.2; 
    scaled_fluctuations = (fluctuations_centered / std(fluctuations_centered, 'omitnan')) * fluctuation_scale;
    base_ambient_temp_unified = max(base_trend + scaled_fluctuations, T_ja_min_ambient); 

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

    % 4.2 计算聚合参数
    ACs_participating = ACs([ACs.ptcp]);
    num_AC_participating = length(ACs_participating);
    
    if num_AC_participating == 0
        AC_Up_Sum = zeros(T_steps_total, 1); 
        AC_Down_Sum = zeros(T_steps_total, 1); 
        AC_Power_Sum = zeros(T_steps_total, 1);
        warning('无 AC 参与，返回全零结果。');
        return;
    end
    
    AggParams = calculateAggregatedACParams(ACs_participating);

    % === [关键修改]：预计算聚合基线功率 ===
    % 目的：将输入的绝对功率指令 P_target 转换为调节指令 Delta_P
    % 公式：P_base_i(t) = (T_ja(t) - Tset_i) / (R_i * eta_i)
    % 聚合：P_base_agg(t) = sum(P_base_i(t))
    % 向量化计算：P_base_agg(t) = T_ja(t) * K1 - K2
    % 其中 K1 = sum(1/(R*eta)), K2 = sum(Tset/(R*eta))
    
    R_vec_p = [ACs_participating.R];
    eta_vec_p = [ACs_participating.eta];
    Tset_vec_p = [ACs_participating.Tset];
    
    % 避免除以零
    valid_params = (R_vec_p > 1e-6) & (eta_vec_p > 1e-6);
    if any(~valid_params)
        warning('发现 %d 台 AC 参数异常(R或eta接近0)，已排除基线计算。', sum(~valid_params));
        R_vec_p = R_vec_p(valid_params);
        eta_vec_p = eta_vec_p(valid_params);
        Tset_vec_p = Tset_vec_p(valid_params);
    end
    
    term_common = 1 ./ (R_vec_p .* eta_vec_p);
    K1 = sum(term_common);
    K2 = sum(Tset_vec_p .* term_common);
    
    % 计算全时段聚合基线功率 (T_steps x 1)
    Agg_Baseline_Power_Series = base_ambient_temp_unified(:) * K1 - K2;
    % 确保基线功率非负 (物理约束)
    Agg_Baseline_Power_Series = max(Agg_Baseline_Power_Series, 0);
    
    fprintf('  [INFO] 聚合基线功率范围: %.2f kW - %.2f kW\n', min(Agg_Baseline_Power_Series), max(Agg_Baseline_Power_Series));
    % ==========================================

    % 4.3 准备 T_ja 矩阵用于后续计算
    T_ja_matrix = repmat(base_ambient_temp_unified', 1, num_AC_participating); 
    Tset_vec = [ACs_participating.Tset];
    R_vec = [ACs_participating.R];
    eta_vec = [ACs_participating.eta];


    %% 5. 主时间循环 (状态化)
    CURRENT_SOC_AC = [ACs_participating.SOC]'; 

    % 结果存储
    Agg_P_Potential_Up_History = zeros(T_steps_total, 1);
    Agg_P_Potential_Down_History = zeros(T_steps_total, 1);
    Individual_Power_History = zeros(T_steps_total, num_AC_participating);

    for t_idx = 1:T_steps_total
        % 1. 聚合状态
        SOC_agg_t = mean(CURRENT_SOC_AC, 'omitnan');

        % 2. 获取并转换指令
        % 输入 P_grid_command_series 是绝对目标功率
        P_target_abs = P_grid_command_series(t_idx);
        P_base_curr = Agg_Baseline_Power_Series(t_idx);
        
        % === [核心逻辑] 计算功率差值 ===
        Delta_P_S_command = P_target_abs - P_base_curr;

        % 3. 预测目标 SOC
        SOC_target_next = AggParams.A * SOC_agg_t + AggParams.B * Delta_P_S_command + AggParams.C;
        SOC_target_next = max(0, min(1, SOC_target_next));

        % 4. 单体响应
        temp_AC_Up_agg = 0;
        temp_AC_Down_agg = 0;
        temp_SOC_next = zeros(num_AC_participating, 1);
        temp_P_achieved = zeros(num_AC_participating, 1);

        parfor i = 1:num_AC_participating
            ac_i = ACs_participating(i);
            soc_curr = CURRENT_SOC_AC(i);
            
            T_ja_val = base_ambient_temp_unified(t_idx);

            P_base_i = ACbaseP_single(T_ja_val, ac_i.Tset, ac_i.R, ac_i.eta);
            
            [P_plus, P_minus] = calculateACAdjustmentPotentia(...
                P_base_i, 2*abs(P_base_i), 0, ...
                ac_i.alpha, ac_i.beta, ac_i.gamma, soc_curr, dt);

            temp_AC_Up_agg = temp_AC_Up_agg + P_plus;
            temp_AC_Down_agg = temp_AC_Down_agg + P_minus;

            % 反解指令 (此时 Delta_Pj 是相对于基线的增量)
            delta_Pj = 0;
            if abs(ac_i.beta) > 1e-9
                delta_Pj = (SOC_target_next - ac_i.alpha * soc_curr - ac_i.gamma) / ac_i.beta;
            end
            delta_Pj_clipped = max(P_minus, min(P_plus, delta_Pj));

            % 更新 SOC
            soc_next = updateACSOC_single(soc_curr, delta_Pj_clipped, ...
                ac_i.alpha, ac_i.beta, ac_i.gamma);

            temp_SOC_next(i) = soc_next;
            temp_P_achieved(i) = delta_Pj_clipped;
        end

        Agg_P_Potential_Up_History(t_idx) = temp_AC_Up_agg;
        Agg_P_Potential_Down_History(t_idx) = temp_AC_Down_agg;
        Individual_Power_History(t_idx, :) = temp_P_achieved';
        CURRENT_SOC_AC = temp_SOC_next;
    end

    % 计算总功率 (基线 + 调节)
    % 基线功率 = (T_ja - Tset) / (R * eta)
    % 注意维度: T_ja 是 (T x 1), Tset 是 (1 x N)
    Baseline_Power = (repmat(base_ambient_temp_unified', 1, num_AC_participating) - repmat(Tset_vec, T_steps_total, 1)) ./ ...
                     (repmat(R_vec, T_steps_total, 1) .* repmat(eta_vec, T_steps_total, 1));
    
    Total_Power_Matrix = Baseline_Power + Individual_Power_History;
    P_standby = 0.05;
    Total_Power_Matrix(Total_Power_Matrix < P_standby) = P_standby; % 物理约束
    
    AC_Power_Sum = sum(Total_Power_Matrix, 2); % 按行求和，得到每个时间步的总功率
    
    % 输出结果
    AC_Up_Sum = Agg_P_Potential_Up_History;
    AC_Down_Sum = Agg_P_Potential_Down_History;
end