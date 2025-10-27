%% 全功能虚拟电厂调节潜力分析系统
clear; close all; clc;
%% 参数生成应用实例
% generateEVData(100, 'residential')  % 生成100辆居民区EV
% generateEVData(50, 'workplace', {'特斯拉Model Y'}) % 生成50辆工作区特斯拉
%% 1. 系统初始化
rng(2023, 'Threefry');                         % 增强随机数生成器
T_total = 24;                                   % 总时长（小时）
dt = 0.05;                                      % 时间分辨率（小时）
time_points = 0:dt:T_total;                     
base_price = 30;            

num_cores = feature('numCores'); % 获取物理核心数
if isempty(gcp('nocreate'))
    parpool('local', num_cores); % 创建与核心数匹配的进程池
else
    disp('已有并行池运行中');
end

%% 2. 初始化参数
% generateEVData(1000, 'residential', 'all');
evFile = 'EV_居民区.xlsx';                 % 电动汽车数据文件

%% 3. 读取设备参数
EVs = initializeEVsFromExcel(evFile);           % 从Excel导入电动汽车参数
num_EV = length(EVs);                           % 获取电动汽车总数

%% 4. 激励响应模块
%% 4.1 参数设定
p_min = 15; p_max = 50;                         % 原始电价范围
p_min_prime = 10; p_max_prime = 40;              % 调整后电价范围

%% 新增模块：激励电价参数扫描
p_incentive_range = linspace(0, 50, 25);        % 生成0-50之间的25个激励电价值
all_EV_Up = zeros(length(time_points), length(p_incentive_range));
all_EV_Down = zeros(length(time_points), length(p_incentive_range));

% ================ 新增：单体数据存储初始化（含P_base_sequence和E_current） ================
% 维度为 [num_EV, 时间点数, 激励电价数]
all_EV_Up_Individual   = zeros(num_EV, length(time_points), length(p_incentive_range));
all_EV_Down_Individual = zeros(num_EV, length(time_points), length(p_incentive_range));
all_SOC_Individual     = zeros(num_EV, length(time_points), length(p_incentive_range));  % 单体SOC
all_P_base_sequence    = zeros(num_EV, length(time_points), length(p_incentive_range));  % 新增：单体基础功率序列
all_E_current          = zeros(num_EV, length(time_points), length(p_incentive_range));  % 新增：单体当前电量
% ==============================================================

%% 参数扫描主循环
parfor (p_idx = 1:length(p_incentive_range), 24)
    current_p = p_incentive_range(p_idx);
    fprintf('\n== 分析激励电价 %.1f 元 ==\n', current_p);
    
    %% 4.2 EV目标电量调整（动态版本）
    EVs_temp = copyEVStruct(EVs);              % 关键深拷贝操作
    
    %% 6.2 目标电量调整（独立随机数流）
    stream = RandStream('Threefry', 'Seed', 2023+p_idx);
    RandStream.setGlobalStream(stream);
    E_tar_max = 0.2 * [EVs_temp.C_EV];
    
    for i = 1:num_EV
        EVs_temp(i).p_incentive = current_p;
        participation = calculateParticipation(current_p, base_price);
        [~, ~, deltaE] = incentiveTempEV(current_p, p_min, p_max, p_min_prime, p_max_prime, E_tar_max(i));
        ptcp_result = (rand() < participation);
        EVs_temp(i).ptcp = ptcp_result;
        
        if EVs_temp(i).ptcp
            EVs_temp(i).E_tar = EVs_temp(i).E_tar - deltaE;
            if EVs_temp(i).E_tar <= EVs_temp(i).E_in
                EVs_temp(i).E_tar = EVs_temp(i).E_tar + deltaE;
            end
        end
    end
    
    %% 5. 预计算模块（完整保留，新增P_base_sequence记录）
    H = 24;                                      % 预测时域（小时）
    H_steps = H / dt;                            % 转换为时间步数
    local_P_base_sequence = zeros(num_EV, length(time_points));  % 本地临时存储P_base_sequence
    for i = 1:num_EV
        % [~, EVs_temp(i).P_base_sequence] = EVbaseP_single_longstep(...
        %     EVs_temp(i).C_EV, EVs_temp(i).eta,...
        %     EVs_temp(i).E_tar, EVs_temp(i).E_in,...
        %     EVs_temp(i).t_dep, EVs_temp(i).t_in, dt, EVs_temp(i).r,...
        %     EVs_temp(i).p_on, EVs_temp(i).SOC, H_steps+1);
        EVs_temp(i).P_base_sequence = EVbaseP_ChargeUntilFull(...
            EVs_temp(i).C_EV, EVs_temp(i).eta,...
            EVs_temp(i).E_tar, EVs_temp(i).E_in,...
            EVs_temp(i).t_dep, EVs_temp(i).t_in, dt, ...
            EVs_temp(i).r, EVs_temp(i).p_on, EVs_temp(i).SOC, H_steps+1);
        local_P_base_sequence(i, :) = EVs_temp(i).P_base_sequence;  % 记录基础功率序列
    end
    
    %% 6. 主时间循环（完整保留，新增E_current记录）
    local_Up = zeros(length(time_points),1);  % 本地临时变量
    local_Down = zeros(length(time_points),1);
    local_SOC = zeros(num_EV, length(time_points));  % 单体SOC
    local_EV_Up_Ind   = zeros(num_EV, length(time_points));  % 单体上调潜力
    local_EV_Down_Ind = zeros(num_EV, length(time_points));  % 单体下调潜力
    local_E_current   = zeros(num_EV, length(time_points));  % 新增：单体当前电量

    for t_idx = 1:length(time_points)
        t = time_points(t_idx);
        fprintf('\n时间 %.1f小时 [激励%.1f元]', t, current_p);
        
        for i = 1:num_EV
            online = (t >= EVs_temp(i).t_in) && (t <= EVs_temp(i).t_dep);
            if online
                EVs_temp(i).P_base = EVs_temp(i).P_base_sequence(t_idx);
                
                [~, ~, m3] = calculateEVABC_single(...
                    EVs_temp(i).C_EV, EVs_temp(i).eta,...
                    EVs_temp(i).E_tar, EVs_temp(i).E_in,...
                    EVs_temp(i).t_dep, EVs_temp(i).t_in, dt, EVs_temp(i).r);
                
                [EVs_temp(i).E_exp, EVs_temp(i).E_current,...
                    EVs_temp(i).P_current, EVs_temp(i).SOC] = ...
                    calculateEVS_single(m3, EVs_temp(i).E_exp,...
                    EVs_temp(i).eta, EVs_temp(i).E_current,...
                    EVs_temp(i).P_current, EVs_temp(i).C_EV,...
                    EVs_temp(i).r, EVs_temp(i).p_on, dt);
                
                [DeltaP_plus, DeltaP_minus] = calculateEVAdjustmentPotentia(...
                    EVs_temp(i).C_EV, EVs_temp(i).r, EVs_temp(i).eta,...
                    EVs_temp(i).E_tar, EVs_temp(i).E_in,...
                    EVs_temp(i).E_current, EVs_temp(i).t_dep,...
                    EVs_temp(i).t_in, EVs_temp(i).p_on, 0,...
                    EVs_temp(i).P_base, EVs_temp(i).SOC, dt);
                
                local_Up(t_idx) = sum(DeltaP_plus);  % 集群潜力累加
                local_Down(t_idx) = sum(DeltaP_minus);

                % 记录单体调节能力
                local_EV_Up_Ind(i, t_idx)   = DeltaP_plus;
                local_EV_Down_Ind(i, t_idx) = DeltaP_minus;
                % 记录单体SOC、当前电量（关键新增）
                local_SOC(i, t_idx) = EVs_temp(i).SOC;
                local_E_current(i, t_idx) = EVs_temp(i).E_current;  % 新增：记录当前电量
            else
                % 离线时保持最后SOC值（或根据需求设为初始值）
                local_SOC(i, t_idx) = EVs_temp(i).SOC;
                local_E_current(i, t_idx) = EVs_temp(i).E_current;  % 新增：离线时保持最后电量
                % 离线时调节能力设为0
                local_EV_Up_Ind(i, t_idx)   = 0;
                local_EV_Down_Ind(i, t_idx) = 0;
            end
        end
    end

    all_EV_Up(:, p_idx) = local_Up;  
    all_EV_Down(:, p_idx) = local_Down;

    % 存储单体数据（关键新增P_base_sequence和E_current）
    all_EV_Up_Individual(:,:,p_idx)   = local_EV_Up_Ind;
    all_EV_Down_Individual(:,:,p_idx) = local_EV_Down_Ind;
    all_SOC_Individual(:,:,p_idx)     = local_SOC;
    all_P_base_sequence(:,:,p_idx)    = local_P_base_sequence;  % 新增：存储基础功率序列
    all_E_current(:,:,p_idx)          = local_E_current;        % 新增：存储当前电量
end

% 整理结果结构（关键新增P_base_sequence和E_current字段）
results_3D.EV_Up = all_EV_Up;
results_3D.EV_Down = all_EV_Down;
results_3D.EV_Up_Individual = all_EV_Up_Individual;
results_3D.EV_Down_Individual = all_EV_Down_Individual;
results_3D.SOC_Individual = all_SOC_Individual;
results_3D.P_base_sequence = all_P_base_sequence;  % 新增：单体基础功率序列
results_3D.E_current = all_E_current;              % 新增：单体当前电量

