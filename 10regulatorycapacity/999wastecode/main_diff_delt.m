%% 全功能虚拟电厂调节潜力分析系统
clear; close all; clc;
%% 参数生成应用实例
% generateEVData(100, 'residential')  % 生成100辆居民区EV
% generateEVData(50, 'workplace', {'特斯拉Model Y'}) % 生成50辆工作区特斯拉
%% 1. 系统初始化
rng(2023, 'Threefry');                         % 增强随机数生成器
T_total = 24;                                   % 总时长（小时）
dt = 0.15;                                      % 时间分辨率（小时）
time_points = 0:dt:T_total;                     
base_price = 30;            


%% 2. 初始化参数
% generateEVData(1000, 'residential', 'all');
evFile = 'EV_居民区.xlsx';                 % 电动汽车数据文件

%% 3. 读取设备参数
EVs = initializeEVsFromExcel(evFile);           % 从Excel导入电动汽车参数，包含 E_tar_original 和 SOC_original
num_EV = length(EVs);                           % 获取电动汽车总数

%% 4. 激励响应模块
%% 4.1 参数设定
p_min = 15; p_max = 50;                         % 原始电价范围
p_min_prime = 10; p_max_prime = 40;              % 调整后电价范围

%% 新增模块：激励电价参数扫描
p_incentive_range = linspace(0, 50, 25);        % 生成0-50之间的25个激励电价值
all_EV_Up = zeros(length(time_points), length(p_incentive_range));
all_EV_Down = zeros(length(time_points), length(p_incentive_range));

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
    fprintf('\n== 分析激励电价 %.1f 元 ==\n', current_p);
    
    % --- BEGIN EV状态重置 ---
    for ev_k = 1:num_EV
        EVs(ev_k).E_tar     = EVs(ev_k).E_tar_original; 
        EVs(ev_k).E_current = EVs(ev_k).E_in;         
        EVs(ev_k).E_exp     = EVs(ev_k).E_in;        
        EVs(ev_k).P_current = 0;                      
        EVs(ev_k).SOC       = EVs(ev_k).SOC_original; 
    end
    % --- END EV状态重置 ---
    
    %% 4.2 EV目标电量调整（动态版本）
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
    
    %% 5. 预计算模块
    H = 24;                                      
    H_steps = H / dt;                            
    current_p_idx_P_base_sequences = zeros(num_EV, length(time_points));

    for i = 1:num_EV
        EVs(i).P_base_sequence = EVbaseP_ChargeUntilFull(...
            EVs(i).C_EV, EVs(i).eta,...
            EVs(i).E_tar, EVs(i).E_in,... 
            EVs(i).t_dep, EVs(i).t_in, dt, ...
            EVs(i).r, EVs(i).p_on, EVs(i).SOC_original, H_steps+1); 
        current_p_idx_P_base_sequences(i, :) = EVs(i).P_base_sequence;
    end
    
    %% 6. 主时间循环
    local_Up_for_p_idx = zeros(length(time_points),1); 
    local_Down_for_p_idx = zeros(length(time_points),1);
    local_SOC_for_p_idx = zeros(num_EV, length(time_points));
    local_EV_Up_Ind_for_p_idx   = zeros(num_EV, length(time_points));
    local_EV_Down_Ind_for_p_idx = zeros(num_EV, length(time_points));
    local_E_current_for_p_idx   = zeros(num_EV, length(time_points));

    for t_idx = 1:length(time_points)
        t = time_points(t_idx);
        fprintf('\n时间 %.1f小时 [激励%.1f元]', t, current_p);
        
        current_t_total_DeltaP_plus = 0;  
        current_t_total_DeltaP_minus = 0; 

        for i = 1:num_EV
            online = (t >= EVs(i).t_in) && (t < EVs(i).t_dep); 
            
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
      
                [DeltaP_plus_dt, DeltaP_minus_dt] = calculateEVAdjustmentPotentia(...
                    EVs(i).C_EV, EVs(i).r, EVs(i).eta,...
                    EVs(i).E_tar, EVs(i).E_in,... 
                    EVs(i).E_current, EVs(i).t_dep,...
                    EVs(i).t_in, EVs(i).p_on, 0,...
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
        local_Up_for_p_idx(t_idx) = current_t_total_DeltaP_plus;
        local_Down_for_p_idx(t_idx) = current_t_total_DeltaP_minus;
    end

    % 存储当前 p_idx 激励价格下的聚合结果和原有单体结果
    all_EV_Up(:, p_idx) = local_Up_for_p_idx;  
    all_EV_Down(:, p_idx) = local_Down_for_p_idx;

    all_EV_Up_Individual(:,:,p_idx)   = local_EV_Up_Ind_for_p_idx;
    all_EV_Down_Individual(:,:,p_idx) = local_EV_Down_Ind_for_p_idx;
    all_SOC_Individual(:,:,p_idx)     = local_SOC_for_p_idx;
    all_P_base_sequence(:,:,p_idx)    = current_p_idx_P_base_sequences; 
    all_E_current(:,:,p_idx)          = local_E_current_for_p_idx;
    
end

% 整理结果结构
results_3D.EV_Up = all_EV_Up;
results_3D.EV_Down = all_EV_Down;
results_3D.EV_Up_Individual = all_EV_Up_Individual;
results_3D.EV_Down_Individual = all_EV_Down_Individual;
results_3D.SOC_Individual = all_SOC_Individual;
results_3D.P_base_sequence = all_P_base_sequence;
results_3D.E_current = all_E_current;
