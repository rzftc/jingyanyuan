function [EV_Up_Sum, EV_Down_Sum, EV_Power_Sum, AggParams_EV, EV_E_Up_Cum, EV_E_Down_Cum] = run_EV_simulation_MC_soc_2bound(random_seed, evFileName)
    % run_EV_simulation_MC_soc
    % 修改说明：
    % 1. 输入参数改为 random_seed，用于控制蒙特卡洛模拟的随机性。
    % 2. 内部 P_grid_command_24 强制设为 0，计算自然状态下的调节潜力边界。
    % 3. 输出 AggParams_EV，包含时变聚合参数 A(t), B(t), C(t)。
    % 4. 输出 EV_E_Up_Cum, EV_E_Down_Cum，包含累积能量边界（能量包络）。
    % 5. 【本次修改】引入 main_potential_agg_ind_bound.m 中的偏差量随机生成逻辑 (actual_delta_E) 和调节能力计算函数 (calculateEVAdjustmentPotentia)。

    P_grid_command_24 = zeros(24, 1); % 强制基准指令为0

    %% 初始化参数
    if nargin > 0
        rng(random_seed); % 使用传入的随机种子
    end

    excelFile = evFileName;
    if ~exist(excelFile, 'file')
        % 如果没有文件，生成一个临时的
        generateEVParameters_real(excelFile, 1000, 1.0);
    end
    [EVs, t_sim, ~, ~, ~] = initializeFromExcel_8am(excelFile);
    
    %% 时间参数
    dt_short = 15;     % 短时间步长 (分钟)
    dt_long = 60;     % 长时间步长 (分钟)
    simulation_start_hour = 8;
    dt = dt_short / 60;      
    dt_minutes = dt_short;    
    dt_long_minutes = dt_long;
    t_adj = 15 / 60;   

    % 时间轴
    time_points_absolute = simulation_start_hour : dt : (simulation_start_hour + t_sim/60 - dt);
    num_time_points = length(time_points_absolute);
    
    num_long_steps = t_sim / dt_long; % 通常为 24
    num_short_per_long = dt_long / dt_short; % 通常为 12
    total_steps = num_long_steps * num_short_per_long;

    % === 处理电网指令 ===
    if length(P_grid_command_24) ~= num_long_steps
        if length(P_grid_command_24) > num_long_steps
            P_tar = P_grid_command_24(1:num_long_steps);
        else
            P_tar = [P_grid_command_24(:); zeros(num_long_steps - length(P_grid_command_24), 1)];
        end
    else
        P_tar = P_grid_command_24;
    end
    % =============================

    num_evs = length(EVs);
    
    % 输出容器
    results_EV_Up_Individual_Sum = zeros(1, total_steps);
    results_EV_Down_Individual_Sum = zeros(1, total_steps);
    results_P_agg = zeros(1, total_steps);
    
    % 能量边界输出容器
    results_EV_E_Up_Individual_Sum = zeros(1, total_steps);
    results_EV_E_Down_Individual_Sum = zeros(1, total_steps);

    %% 初始化参数 (向量化)
    base_Price_vec = 30 * ones(num_evs, 1);
    p_incentive_vec = [EVs.p_incentive]';
    participation_probabilities_vec = calculateParticipation(p_incentive_vec, base_Price_vec);
    ptcp_vec = (rand(num_evs, 1) < participation_probabilities_vec);

    E_tar_max_flex_vec = 0.2 * [EVs.C]';
    p_min = 15; p_max = 50; p_min_prime = 10; p_max_prime = 40;
    [deltaE_up_vec, deltaE_down_vec] = incentiveTempEV_updown(p_incentive_vec, p_min, p_max, p_min_prime, p_max_prime, E_tar_max_flex_vec);
    
    % === [修改点1]：随机生成 actual_delta_E 并更新 E_reg_min/max ===
    % 参考 main_potential_agg_ind_bound.m 的逻辑
    min_safe_r = 0.05;
    actual_delta_E_vec = zeros(num_evs, 1);
    
    % 仅对参与车辆生成随机偏差
    idx_ptcp = find(ptcp_vec);
    if ~isempty(idx_ptcp)
        d_up = deltaE_up_vec(idx_ptcp);
        d_down = deltaE_down_vec(idx_ptcp);
        % 生成随机实际偏差: actual = down + (up - down) * rand
        actual_delta_E_vec(idx_ptcp) = d_down + (d_up - d_down) .* rand(length(idx_ptcp), 1);
    end

    % 更新结构体
    for i=1:num_evs
        EVs(i).ptcp = ptcp_vec(i);
        EVs(i).E_tar_original = EVs(i).E_tar_set;
        EVs(i).SOC_original = 0;
        
        if EVs(i).ptcp
            % 获取计算好的随机偏差
            actual_delta_E = actual_delta_E_vec(i);
            
            % 更新 E_reg_min 和 E_reg_max (对称应用偏差)
            EVs(i).E_reg_min = EVs(i).E_tar_original - actual_delta_E;
            EVs(i).E_reg_max = EVs(i).E_tar_original + actual_delta_E;
            
            % 物理边界截断
            EVs(i).E_reg_min = max(EVs(i).E_reg_min, EVs(i).E_ini);
            EVs(i).E_reg_max = min(EVs(i).E_reg_max, EVs(i).C);
            
            % 修正 r (能量死区比例)
            EVs(i).r = max(actual_delta_E / EVs(i).C, min_safe_r);
            
            % 设定目标电量 (保持在合理范围内，这里保持原样或取 E_tar_original)
            % 在 bound 代码中通常保持 E_tar_original，此处不做额外随机化以保持一致
            EVs(i).E_tar_set = EVs(i).E_tar_original; 
        else
            EVs(i).E_reg_min = EVs(i).E_tar_original;
            EVs(i).E_reg_max = EVs(i).E_tar_original;
            EVs(i).r = 0.05; % 默认值
        end
        
        % 初始 E_tar 设为 max(E_tar_set, E_ini)
        EVs(i).E_tar = max(EVs(i).E_tar_set, EVs(i).E_ini);
        
        if EVs(i).P_N > 0
            EVs(i).tau_rem = 60 * (EVs(i).E_tar - EVs(i).E_ini) / EVs(i).P_N;
        else
            EVs(i).tau_rem = 0;
        end
        EVs(i).E_exp = EVs(i).E_ini;
        EVs(i).E_actual = EVs(i).E_ini;
        EVs(i).E_current = EVs(i).E_ini;
        EVs(i).P_current = 0;
        
        % 确保 r 存在 (上面已赋值，这里作为兜底)
        if ~isfield(EVs(i), 'r')
            EVs(i).r = 0.05;
        end
    end
    % ========================================================

    S_agg_current = 0;

    %% 预计算基线
    EVs_for_baseline = EVs;
    parfor i = 1:num_evs
        t_dep_h = EVs_for_baseline(i).t_dep / 60;
        t_in_h = EVs_for_baseline(i).t_in / 60;
            
        EVs_for_baseline(i).P_base_sequence = EVbaseP_Hysteresis(...
            EVs_for_baseline(i).C, EVs_for_baseline(i).eta,...
            EVs_for_baseline(i).E_tar, EVs_for_baseline(i).E_ini,...
            t_dep_h, t_in_h, dt, ...
            EVs_for_baseline(i).r, EVs_for_baseline(i).P_N, ...
            EVs_for_baseline(i).SOC_original, num_time_points, time_points_absolute);
        % EVs_for_baseline(i).P_base_sequence = EVbaseP_ChargeUntilFull(...
        %     EVs_for_baseline(i).C, EVs_for_baseline(i).eta,...
        %     EVs_for_baseline(i).E_tar, EVs_for_baseline(i).E_ini,...
        %     t_dep_h, t_in_h, dt, ...
        %     EVs_for_baseline(i).r, EVs_for_baseline(i).P_N, ...
        %     EVs_for_baseline(i).SOC_original, num_time_points, time_points_absolute);
       
    end
    for i = 1:num_evs
        EVs(i).P_base_sequence = EVs_for_baseline(i).P_base_sequence;
    end
    clear EVs_for_baseline;

    %% 预计算聚合参数 Kappa (向量化)
    C_vec = [EVs.C]';
    eta_vec = [EVs.eta]';
    r_vec = [EVs.r]';
    t_in_h_vec = [EVs.t_in]' / 60;
    t_dep_h_vec = [EVs.t_dep]' / 60;
    E_tar_vec = [EVs.E_tar]';
    E_in_vec = [EVs.E_ini]';
    ptcp_vec = [EVs.ptcp]'; 

    % 计算 Kappa1, Kappa2 
    Kappa1_all = - (C_vec .* r_vec) ./ (eta_vec .* dt);
    Kappa2_all =   (C_vec .* r_vec) ./ (eta_vec .* dt);

    % 计算 Kappa3 (P_req)
    duration_vec = t_dep_h_vec - t_in_h_vec;
    duration_vec(duration_vec < dt) = dt; 
    Kappa3_all = (E_tar_vec - E_in_vec) ./ (eta_vec .* duration_vec);

    Agg_A_trace = ones(total_steps, 1);
    Agg_B_trace = zeros(total_steps, 1);
    Agg_C_trace = zeros(total_steps, 1);

    %% 主仿真循环
    for long_idx = 1:num_long_steps
        t_long_start_minute = (long_idx - 1) * dt_long_minutes;

        % 长时间步：聚合与指令分配
        [lambda_star] = aggregateEVs_mont(EVs, P_tar(long_idx));
        [~, S_agg_next] = calculateVirtualSOC_agg(EVs, dt_long_minutes);

        % 短时间步：状态更新
        for short_idx = 1:num_short_per_long
            step_idx = (long_idx - 1) * num_short_per_long + short_idx;
            t_relative_minute = t_long_start_minute + (short_idx - 1) * dt_minutes; 
            t_current_minute_abs = t_relative_minute + simulation_start_hour * 60; 
            current_absolute_hour = time_points_absolute(step_idx);

            % 计算当前时刻的聚合参数
            is_connected = (t_in_h_vec <= current_absolute_hour) & (t_dep_h_vec > current_absolute_hour) & ptcp_vec;
            
            if any(is_connected)
                Sum_K1 = sum(Kappa1_all(is_connected));
                Sum_K2 = sum(Kappa2_all(is_connected));
                Sum_K3 = sum(Kappa3_all(is_connected));
                
                if abs(Sum_K1) > 1e-6
                    Agg_A_trace(step_idx) = - Sum_K2 / Sum_K1;
                    Agg_B_trace(step_idx) =   1 / Sum_K1;
                    Agg_C_trace(step_idx) = - Sum_K3 / Sum_K1;
                else
                    Agg_A_trace(step_idx) = 1; Agg_B_trace(step_idx) = 0; Agg_C_trace(step_idx) = 0;
                end
            else
                Agg_A_trace(step_idx) = 1; Agg_B_trace(step_idx) = 0; Agg_C_trace(step_idx) = 0;
            end

            temp_delta_p_plus_individual = zeros(num_evs, 1);
            temp_delta_p_minus_individual = zeros(num_evs, 1);
            temp_P_current = zeros(num_evs, 1);
            
            % 临时能量边界容器
            temp_e_up_individual = zeros(num_evs, 1);
            temp_e_down_individual = zeros(num_evs, 1);
            
            EVs_in_parfor = EVs;

            parfor i = 1:num_evs
                EV = EVs_in_parfor(i);
                EV = updateLockState_unlock(EV, t_current_minute_abs);

                EV_temp_with_handle = generateDemandCurve(EV);
                current_P_val = 0;
                if isfield(EV_temp_with_handle, 'demandCurve') && isa(EV_temp_with_handle.demandCurve, 'function_handle')
                    current_P_val = EV_temp_with_handle.demandCurve(lambda_star);
                else
                     switch EV.state
                        case 'LockON'
                            current_P_val = EV.P_N;
                        case {'LockOFF', 'OFF', 'ON'}
                            current_P_val = 0;
                        otherwise
                            current_P_val = 0;
                     end
                end
                EV.P_current = current_P_val;

                EV = calculateVirtualSOC_upgrade(EV, t_current_minute_abs, dt_minutes);

                temp_P_current(i) = EV.P_current;
                
                % 计算单体潜力
                is_online_h = (current_absolute_hour >= (EV.t_in / 60)) && (current_absolute_hour < (EV.t_dep / 60)); 
                
                if EV.ptcp && is_online_h 
                     P_base_i = EV.P_base_sequence(step_idx); 
                     t_dep_h = EV.t_dep / 60; 
                     t_in_h = EV.t_in / 60;

                     % === [修改点2]：使用 calculateEVAdjustmentPotentia 计算调节能力 ===
                     P_max_i = EV.P_N;
                     P_min_i = 0; 
                     
                     [DeltaP_plus_i, DeltaP_minus_i] = calculateEVAdjustmentPotentia(...
                         EV.C, EV.r, EV.eta, EV.E_tar, EV.E_ini, EV.E_actual, ...
                         t_dep_h, t_in_h, P_max_i, P_min_i, P_base_i, ...
                         EV.S_modified, t_adj);
                     % ==============================================================
                     
                     temp_delta_p_plus_individual(i) = DeltaP_plus_i;
                     temp_delta_p_minus_individual(i) = DeltaP_minus_i;
                     
                     % 2. 能量边界计算 (修改为基于 S_modified 在 [-1, 1] 的动态边界)
                     
                     % 计算当前时刻允许的能量上下限 (S_modified = +/- 1 对应的能量)
                     % 上限 E_upper = E_exp + C*r (且不超过电池容量 C)
                     E_bound_upper_dynamic = min(EV.C, EV.E_exp + EV.C * EV.r);
                     
                     % 下限 E_lower = E_exp - C*r (且不低于0)
                     E_bound_lower_dynamic = max(EV.E_ini, EV.E_exp - EV.C * EV.r);

                     % E_Up: 在满足 S>=-1 约束下，还能多充多少
                     temp_e_up_individual(i) = max(0, E_bound_upper_dynamic - EV.E_actual);
                     
                     % E_Down: 在满足 S<=1 约束下，还能多放多少
                     temp_e_down_individual(i) = max(0, EV.E_actual - E_bound_lower_dynamic);
                else
                     temp_delta_p_plus_individual(i) = 0;
                     temp_delta_p_minus_individual(i) = 0;
                     temp_e_up_individual(i) = 0;
                     temp_e_down_individual(i) = 0;
                end

                EVs_in_parfor(i) = EV;
            end
            EVs = EVs_in_parfor;
            
            % 记录当前步结果
            results_P_agg(step_idx) = sum(temp_P_current);
            results_EV_Up_Individual_Sum(step_idx) = sum(temp_delta_p_plus_individual);
            results_EV_Down_Individual_Sum(step_idx) = sum(temp_delta_p_minus_individual);
            
            % 记录能量边界结果
            results_EV_E_Up_Individual_Sum(step_idx) = sum(temp_e_up_individual);
            results_EV_E_Down_Individual_Sum(step_idx) = sum(temp_e_down_individual);
        end
        
        S_agg_current = S_agg_next;
    end

    % 输出结果
    EV_Up_Sum = results_EV_Up_Individual_Sum(:);
    EV_Down_Sum = results_EV_Down_Individual_Sum(:);
    EV_Power_Sum = results_P_agg(:);
    
    % 输出累积能量边界
    EV_E_Up_Cum = results_EV_E_Up_Individual_Sum(:);
    EV_E_Down_Cum = results_EV_E_Down_Individual_Sum(:);
    
    % 组装聚合参数输出
    AggParams_EV.A = Agg_A_trace;
    AggParams_EV.B = Agg_B_trace;
    AggParams_EV.C = Agg_C_trace;
    AggParams_EV.C_total = sum(C_vec(ptcp_vec)); 
end