%% 全功能虚拟电厂调节潜力分析系统 (理论最终修正版 + 在线车辆统计)
clear; close all; clc;
%% 1. 系统初始化
rng(2023, 'Threefry');

% --- 时间参数定义 ---
% 直接定义仿真的绝对开始和结束时间（单位：小时）
simulation_start_hour = 6;  % 第一天早上 6 点
simulation_end_hour   = 30; % 第二天早上 6 点 (24 + 6)

dt = 5/60; % 时间步长（小时），例如 5/60 代表5分钟
t_adj = 60/60; % 调节时长

% 基于开始和结束时间生成仿真时间点
time_points_absolute = simulation_start_hour:dt:simulation_end_hour;
T_total = simulation_end_hour - simulation_start_hour; % 计算得出总时长

% --- 其他初始化 ---
base_price = 30;

%% 2. 初始化参数
evFile = 'EV_residential.xlsx';
EVs = initializeEVsFromExcel(evFile);           
num_EV = length(EVs);                           

%% 4. 激励响应模块
p_min = 15; p_max = 50;                         
p_min_prime = 10; p_max_prime = 40;              
p_incentive_range = linspace(0, 50, 25);        
all_EV_Up = zeros(length(time_points_absolute), length(p_incentive_range));
all_EV_Down = zeros(length(time_points_absolute), length(p_incentive_range));
all_EV_Total_Charge_Power = zeros(length(time_points_absolute), length(p_incentive_range));

% ================ 数据存储初始化 ================
all_EV_Up_Individual   = zeros(num_EV, length(time_points_absolute), length(p_incentive_range));
all_EV_Down_Individual = zeros(num_EV, length(time_points_absolute), length(p_incentive_range));
all_SOC_Individual     = zeros(num_EV, length(time_points_absolute), length(p_incentive_range));
all_P_base_sequence    = zeros(num_EV, length(time_points_absolute), length(p_incentive_range));
all_E_current          = zeros(num_EV, length(time_points_absolute), length(p_incentive_range));
all_Online_EV_Count = zeros(length(time_points_absolute), length(p_incentive_range));
% ==============================================================

%% 参数扫描主循环
for p_idx = 1:length(p_incentive_range)
    current_p = p_incentive_range(p_idx);
    fprintf('\n== 分析激励电价 %.1f 元 (dt = %.2f 小时)==\n', current_p, dt);
    
    % --- 状态重置 ---
    for ev_k = 1:num_EV
        EVs(ev_k).E_current = EVs(ev_k).E_in;         
        EVs(ev_k).E_exp     = EVs(ev_k).E_in;        
        EVs(ev_k).P_current = 0;                      
        EVs(ev_k).SOC       = EVs(ev_k).SOC_original;
        EVs(ev_k).E_tar     = EVs(ev_k).E_tar_original;
    end
    
    stream = RandStream('Threefry', 'Seed', 2023+p_idx); 
    RandStream.setGlobalStream(stream); 
    
    % === 计算灵活性窗口 ===
    for i = 1:num_EV 
        EVs(i).p_incentive = current_p; 
        participation = calculateParticipation(current_p, base_price);
        EVs(i).ptcp = (rand() < participation);
        
        E_flex_max = 0.2 * EVs(i).C_EV;
        [deltaE_up, deltaE_down] = incentiveTempEV_updown(current_p, p_min, p_max, p_min_prime, p_max_prime, E_flex_max);
        
        if EVs(i).ptcp
            EVs(i).E_reg_min = EVs(i).E_tar_original - deltaE_down;
             if EVs(i).E_reg_min<=EVs(i).E_in
                 EVs(i).E_reg_min = EVs(i).E_in;
             end
            EVs(i).E_reg_max = EVs(i).E_tar_original + deltaE_up;
             if EVs(i).E_reg_max >= EVs(i).C_EV
                 EVs(i).E_reg_max = EVs(i).C_EV;
             end
        else
            EVs(i).E_reg_min = EVs(i).E_tar_original;
            EVs(i).E_reg_max = EVs(i).E_tar_original;
        end
    end
    
    num_time_points_for_baseline = length(time_points_absolute);
    current_p_idx_P_base_sequences = zeros(num_EV, length(time_points_absolute));

    for i = 1:num_EV
        % 注意：t_in 和 t_dep 是绝对时间，P_base_sequence 需要知道仿真的绝对时间点
        EVs(i).P_base_sequence = EVbaseP_ChargeUntilFull(...
            EVs(i).C_EV, EVs(i).eta,...
            EVs(i).E_tar_original, EVs(i).E_in,... 
            EVs(i).t_dep, EVs(i).t_in, dt, ... 
            EVs(i).r, EVs(i).p_on, EVs(i).SOC_original, num_time_points_for_baseline, time_points_absolute);
        current_p_idx_P_base_sequences(i, :) = EVs(i).P_base_sequence;
    end
    
    local_Up_for_p_idx = zeros(length(time_points_absolute),1); 
    local_Down_for_p_idx = zeros(length(time_points_absolute),1);
    local_Total_Charge_Power_for_p_idx = zeros(length(time_points_absolute),1);
    local_Online_EV_Count_for_p_idx = zeros(length(time_points_absolute), 1);
    local_SOC_for_p_idx = zeros(num_EV, length(time_points_absolute));
    local_EV_Up_Ind_for_p_idx   = zeros(num_EV, length(time_points_absolute));
    local_EV_Down_Ind_for_p_idx = zeros(num_EV, length(time_points_absolute));
    local_E_current_for_p_idx   = zeros(num_EV, length(time_points_absolute));
    
    for t_idx = 1:length(time_points_absolute)
        current_absolute_hour = time_points_absolute(t_idx); 
        
        if mod(t_idx-1, round(1/dt)) == 0 
            fprintf('\n仿真绝对时间 %.2f 小时, 激励电价 %.1f 元...', current_absolute_hour, current_p);
        end
        
        current_t_total_DeltaP_plus = 0;  
        current_t_total_DeltaP_minus = 0; 
        current_t_total_charge_power_this_step = 0;
        current_t_online_ev_count = 0;

        for i = 1:num_EV
            online = (current_absolute_hour >= EVs(i).t_in) && (current_absolute_hour < EVs(i).t_dep); 
           if EVs(i).ptcp
            if online
                current_t_online_ev_count = current_t_online_ev_count + 1;
                
                EVs(i).P_base = EVs(i).P_base_sequence(t_idx); 
                
                [~, ~, m3] = calculateEVABC_single(...
                    EVs(i).C_EV, EVs(i).eta,...
                    EVs(i).E_tar_original, EVs(i).E_in,...
                    EVs(i).t_dep, EVs(i).t_in, dt, EVs(i).r);
                
                [EVs(i).E_exp, EVs(i).E_current,...
                    EVs(i).P_current, EVs(i).SOC] = ...
                    calculateEVS_single(m3, EVs(i).E_exp,EVs(i).E_tar_original,...
                    EVs(i).eta, EVs(i).E_current,...
                    EVs(i).P_current, EVs(i).C_EV,...
                    EVs(i).r, EVs(i).p_on, dt,EVs(i).t_dep,current_absolute_hour);
                
                if EVs(i).P_current > 1e-3
                    current_t_total_charge_power_this_step = current_t_total_charge_power_this_step + EVs(i).P_current;
                end
      
                [DeltaP_plus_dt, DeltaP_minus_dt] = calculateEVAdjustmentPotentia_new(...
                    EVs(i).E_reg_min, EVs(i).E_reg_max, ...
                    EVs(i).E_current, EVs(i).t_dep, current_absolute_hour, ...
                    EVs(i).p_on, EVs(i).P_base, EVs(i).eta, t_adj);
                
                current_t_total_DeltaP_plus = current_t_total_DeltaP_plus + DeltaP_plus_dt;
                current_t_total_DeltaP_minus = current_t_total_DeltaP_minus + DeltaP_minus_dt;

                local_EV_Up_Ind_for_p_idx(i, t_idx)   = DeltaP_plus_dt;
                local_EV_Down_Ind_for_p_idx(i, t_idx) = DeltaP_minus_dt;
                local_SOC_for_p_idx(i, t_idx) = EVs(i).SOC;
                local_E_current_for_p_idx(i, t_idx) = EVs(i).E_current;
            else
                if t_idx > 1 
                    local_SOC_for_p_idx(i, t_idx) = local_SOC_for_p_idx(i, t_idx-1);
                    local_E_current_for_p_idx(i, t_idx) = local_E_current_for_p_idx(i, t_idx-1);
                else 
                    local_SOC_for_p_idx(i, t_idx) = EVs(i).SOC_original; 
                    local_E_current_for_p_idx(i, t_idx) = EVs(i).E_in;      
                end
                local_EV_Up_Ind_for_p_idx(i, t_idx)   = 0; 
                local_EV_Down_Ind_for_p_idx(i, t_idx) = 0; 
            end
           end
        end
        local_Up_for_p_idx(t_idx) = current_t_total_DeltaP_plus;
        local_Down_for_p_idx(t_idx) = current_t_total_DeltaP_minus;
        local_Total_Charge_Power_for_p_idx(t_idx) = current_t_total_charge_power_this_step;
        local_Online_EV_Count_for_p_idx(t_idx) = current_t_online_ev_count;
    end
    fprintf(' 完成.\n'); 

    all_EV_Up(:, p_idx) = local_Up_for_p_idx;  
    all_EV_Down(:, p_idx) = local_Down_for_p_idx;
    all_EV_Total_Charge_Power(:, p_idx) = local_Total_Charge_Power_for_p_idx;
    all_Online_EV_Count(:, p_idx) = local_Online_EV_Count_for_p_idx;

    all_EV_Up_Individual(:,:,p_idx)   = local_EV_Up_Ind_for_p_idx;
    all_EV_Down_Individual(:,:,p_idx) = local_EV_Down_Ind_for_p_idx;
    all_SOC_Individual(:,:,p_idx)     = local_SOC_for_p_idx;
    all_P_base_sequence(:,:,p_idx)    = current_p_idx_P_base_sequences; 
    all_E_current(:,:,p_idx)          = local_E_current_for_p_idx;
end

% 整理结果结构
results_3D.EV_Up = all_EV_Up;
results_3D.EV_Down = all_EV_Down;
results_3D.EV_Total_Charge_Power = all_EV_Total_Charge_Power;
results_3D.EV_Up_Individual = all_EV_Up_Individual;
results_3D.EV_Down_Individual = all_EV_Down_Individual;
results_3D.SOC_Individual = all_SOC_Individual;
results_3D.P_base_sequence = all_P_base_sequence;
results_3D.E_current = all_E_current;
results_3D.dt_actual = dt; 
results_3D.time_points_actual = time_points_absolute; 
results_3D.simulation_start_hour_actual = simulation_start_hour; 
results_3D.p_incentive_range_actual = p_incentive_range;
results_3D.Online_EV_Count = all_Online_EV_Count;