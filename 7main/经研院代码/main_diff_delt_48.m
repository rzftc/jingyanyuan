%% 全功能虚拟电厂调节潜力分析系统
clear; close all; clc;
%% 1. 系统初始化
rng(2023, 'Threefry');
T_total = 24;                                   % 仿真总时长（小时）

dt = 5/60; % <<<<<<------示例：根据需要修改此dt值
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
time_points = 0:dt:T_total; % 仿真的相对时间点 (0 到 T_total)                   
base_price = 30;            

simulation_start_offset_hours = 6.0; % 仿真 t=0 对应于实际的早上6点

%% 2. 初始化参数
evFile = 'EV_residential.xlsx'; % 在您的实际脚本中，这里可能是 EV_residential.xlsx
EVs = initializeEVsFromExcel(evFile);           
num_EV = length(EVs);                           

%% 4. 激励响应模块
p_min = 15; p_max = 50;                         
p_min_prime = 10; p_max_prime = 40;              
p_incentive_range = linspace(0, 50, 25);        
all_EV_Up = zeros(length(time_points), length(p_incentive_range));
all_EV_Down = zeros(length(time_points), length(p_incentive_range));
all_EV_Total_Charge_Power = zeros(length(time_points), length(p_incentive_range)); % 新增：存储总充电功率

% ================ 数据存储初始化 ================
all_EV_Up_Individual   = zeros(num_EV, length(time_points), length(p_incentive_range));
all_EV_Down_Individual = zeros(num_EV, length(time_points), length(p_incentive_range));
all_SOC_Individual     = zeros(num_EV, length(time_points), length(p_incentive_range));
all_P_base_sequence    = zeros(num_EV, length(time_points), length(p_incentive_range));
all_E_current          = zeros(num_EV, length(time_points), length(p_incentive_range));
% ==============================================================

%% 参数扫描主循环
for p_idx = 1:length(p_incentive_range)
    current_p = p_incentive_range(p_idx);
    fprintf('\n== 分析激励电价 %.1f 元 (dt = %.2f 小时)==\n', current_p, dt);
    
    for ev_k = 1:num_EV
        EVs(ev_k).E_tar     = EVs(ev_k).E_tar_original; 
        EVs(ev_k).E_current = EVs(ev_k).E_in;         
        EVs(ev_k).E_exp     = EVs(ev_k).E_in;        
        EVs(ev_k).P_current = 0;                      
        EVs(ev_k).SOC       = EVs(ev_k).SOC_original; 
    end
    
    stream = RandStream('Threefry', 'Seed', 2023+p_idx); 
    RandStream.setGlobalStream(stream); 
    
    E_tar_max_all_evs = zeros(1, num_EV); 
    for i_ev_init = 1:num_EV
        E_tar_max_all_evs(i_ev_init) = 0.2 * EVs(i_ev_init).C_EV;
    end
    
    for i = 1:num_EV 
        EVs(i).p_incentive = current_p; 
        participation = calculateParticipation(current_p, base_price);
        [~, ~, deltaE] = incentiveTempEV(current_p, p_min, p_max, p_min_prime, p_max_prime, E_tar_max_all_evs(i));
        EVs(i).ptcp = (rand() < participation);
        
        if EVs(i).ptcp
            adjusted_E_tar = EVs(i).E_tar_original - deltaE; 
            if adjusted_E_tar <= EVs(i).E_in
                EVs(i).E_tar = EVs(i).E_tar_original; 
            else
                EVs(i).E_tar = adjusted_E_tar;
            end
        else
            EVs(i).E_tar = EVs(i).E_tar_original; 
        end
    end
    
    % H = T_total; % 在您的代码中 H 是 24
    % H_steps = H / dt;  
    num_time_points_for_baseline = length(time_points); % 基线序列的长度与仿真窗口的时间点数一致

    current_p_idx_P_base_sequences = zeros(num_EV, length(time_points));

    for i = 1:num_EV
        t_in_sim_relative = EVs(i).t_in - simulation_start_offset_hours;
        t_dep_sim_relative = EVs(i).t_dep - simulation_start_offset_hours;

        EVs(i).P_base_sequence = EVbaseP_ChargeUntilFull(...
            EVs(i).C_EV, EVs(i).eta,...
            EVs(i).E_tar_original, EVs(i).E_in,... 
            t_dep_sim_relative, t_in_sim_relative, dt, ... 
            EVs(i).r, EVs(i).p_on, EVs(i).SOC_original, num_time_points_for_baseline); % 使用 length(time_points)
        current_p_idx_P_base_sequences(i, :) = EVs(i).P_base_sequence;
    end
    
    local_Up_for_p_idx = zeros(length(time_points),1); 
    local_Down_for_p_idx = zeros(length(time_points),1);
    local_Total_Charge_Power_for_p_idx = zeros(length(time_points),1); % 新增：当前p_idx的总充电功率时间序列

    local_SOC_for_p_idx = zeros(num_EV, length(time_points));
    local_EV_Up_Ind_for_p_idx   = zeros(num_EV, length(time_points));
    local_EV_Down_Ind_for_p_idx = zeros(num_EV, length(time_points));
    local_E_current_for_p_idx   = zeros(num_EV, length(time_points));
    
    for t_idx = 1:length(time_points)
        t_sim_relative = time_points(t_idx); 
        current_absolute_hour = t_sim_relative + simulation_start_offset_hours;
        
        if mod(t_idx-1, round(1/dt)) == 0 
            fprintf('\n仿真相对时间 %.2f 小时 (实际时间 %.2f), 激励电价 %.1f 元...', t_sim_relative, current_absolute_hour, current_p);
        end
        
        current_t_total_DeltaP_plus = 0;  
        current_t_total_DeltaP_minus = 0; 
        current_t_total_charge_power_this_step = 0; % 新增：当前时间步的总充电功率

        for i = 1:num_EV
            online = (current_absolute_hour >= EVs(i).t_in) && (current_absolute_hour < EVs(i).t_dep); 
           if EVs(i).ptcp
            if online
                EVs(i).P_base = EVs(i).P_base_sequence(t_idx); 
                
                [~, ~, m3] = calculateEVABC_single(...
                    EVs(i).C_EV, EVs(i).eta,...
                    EVs(i).E_tar, EVs(i).E_in,... 
                    EVs(i).t_dep, EVs(i).t_in, dt, EVs(i).r);
                
                [EVs(i).E_exp, EVs(i).E_current,...
                    EVs(i).P_current, EVs(i).SOC] = ...
                    calculateEVS_single(m3, EVs(i).E_exp,...
                    EVs(i).eta, EVs(i).E_current,...
                    EVs(i).P_current, EVs(i).C_EV,...
                    EVs(i).r, EVs(i).p_on, dt);
                
                % 新增：累加当前EV的充电功率 (P_current > 0 表示正在充电)
                if EVs(i).P_current > 1e-3 % 加一个小的阈值避免浮点误差
                    current_t_total_charge_power_this_step = current_t_total_charge_power_this_step + EVs(i).P_current;
                end
                % --- 结束新增 ---
      
                [DeltaP_plus_dt, DeltaP_minus_dt] = calculateEVAdjustmentPotentia(...
                    EVs(i).C_EV, EVs(i).r, EVs(i).eta,...
                    EVs(i).E_tar, EVs(i).E_in,... 
                    EVs(i).E_current, EVs(i).t_dep, EVs(i).t_in, ... 
                    EVs(i).p_on, 0,...
                    EVs(i).P_base, EVs(i).SOC, dt);
                
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
        local_Total_Charge_Power_for_p_idx(t_idx) = current_t_total_charge_power_this_step; % 新增：存储当前时间步的总充电功率
    end
    fprintf(' 完成.\n'); 

    all_EV_Up(:, p_idx) = local_Up_for_p_idx;  
    all_EV_Down(:, p_idx) = local_Down_for_p_idx;
    all_EV_Total_Charge_Power(:, p_idx) = local_Total_Charge_Power_for_p_idx; % 新增：存储当前激励价格下的总充电功率序列

    all_EV_Up_Individual(:,:,p_idx)   = local_EV_Up_Ind_for_p_idx;
    all_EV_Down_Individual(:,:,p_idx) = local_EV_Down_Ind_for_p_idx;
    all_SOC_Individual(:,:,p_idx)     = local_SOC_for_p_idx;
    all_P_base_sequence(:,:,p_idx)    = current_p_idx_P_base_sequences; 
    all_E_current(:,:,p_idx)          = local_E_current_for_p_idx;
end

% 整理结果结构
results_3D.EV_Up = all_EV_Up;
results_3D.EV_Down = all_EV_Down;
results_3D.EV_Total_Charge_Power = all_EV_Total_Charge_Power; % 新增：将总充电功率添加到结果结构体
results_3D.EV_Up_Individual = all_EV_Up_Individual;
results_3D.EV_Down_Individual = all_EV_Down_Individual;
results_3D.SOC_Individual = all_SOC_Individual;
results_3D.P_base_sequence = all_P_base_sequence;
results_3D.E_current = all_E_current;
results_3D.dt_actual = dt; 
results_3D.time_points_actual = time_points; 
results_3D.simulation_start_offset_hours = simulation_start_offset_hours; 
results_3D.p_incentive_range_actual = p_incentive_range; 
