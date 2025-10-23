%% 全功能虚拟电厂调节潜力分析系统
clear; close all; clc;

%% 1. 系统初始化
rng(2023);                                      % 固定随机种子，保证结果可重复
T_total = 24;                                   % 总时长（小时）
dt = 5/60;                                      % 时间分辨率（小时）
time_points = 0:dt:T_total;                     % 生成时间序列
base_price = 30;                                % 基础电价（元/kWh）

%% 2. 初始化参数
% generateResidentialEV(1000);                     % 生成居民区电动汽车数据
% generateWorkplaceEV(1000);                       % 生成工作区电动汽车数据
evFile = 'EV_template.xlsx';                 % 电动汽车数据文件
% generateExampleExcel_real_24(100,100, 0.6);      
acFile = 'AC_template.xlsx';             
%% 3. 读取设备参数
ACs = initializeACsFromExcel(acFile);          
EVs = initializeEVsFromExcel(evFile);
num_AC = length(ACs);                          
num_EV = length(EVs);           

%% 4. 激励响应模块
%% 4.1 参数设定
p_min = 15; p_max = 50;                         % 原始电价范围
p_min_prime = 10; p_max_prime = 40; T_set_max = 3;                % 调整后电价范围

%% AC温度设定调整 - 修改为使用临时结构体数组
temp_ACs = ACs;  % 创建临时结构体数组
parfor i = 1:num_AC
    participation = calculateParticipation(temp_ACs(i).p_incentive, base_price);
    [~, ~, deltaT] = incentiveTempAC(...        
        temp_ACs(i).p_incentive, p_min, p_max, p_min_prime, p_max_prime, T_set_max);
    
    temp_ACs(i).ptcp = (rand() < participation);     
    
    if temp_ACs(i).ptcp
        temp_ACs(i).Tmax = temp_ACs(i).Tset + deltaT;
        temp_ACs(i).Tmin = temp_ACs(i).Tset - deltaT;
    end
    
    base_temp = temp_ACs(i).Tset + 4*sin(2*pi*time_points/24); 
    temp_range = temp_ACs(i).Tmax - temp_ACs(i).Tmin;
    noise = 0.2 * temp_range * randn(size(time_points));
    temp_ACs(i).T_ja = min(max(base_temp + noise, temp_ACs(i).Tmin), temp_ACs(i).Tmax);
end
ACs = temp_ACs;  % 更新原始结构体数组

%% 4.2 EV目标电量调整 - 修改为使用临时结构体数组
E_tar_max = 0.2 * [EVs.C_EV];                   % 计算最大目标电量变化范围
temp_EVs = EVs;  % 创建临时结构体数组
parfor i = 1:num_EV
    temp_EVs(i).p_incentive = 11;                     % 设置激励电价
    participation = calculateParticipation(temp_EVs(i).p_incentive, base_price);
    [~, ~, deltaE] = incentiveTempEV(temp_EVs(i).p_incentive, p_min, p_max, p_min_prime, p_max_prime, E_tar_max(i));
    ptcp_result = (rand() < participation);      % 生成参与决策
    temp_EVs(i).ptcp = ptcp_result;
    
    if temp_EVs(i).ptcp
        temp_EVs(i).E_tar = temp_EVs(i).E_tar - deltaE;
        % 防止目标电量低于初始电量
        if temp_EVs(i).E_tar <= temp_EVs(i).E_in
            temp_EVs(i).E_tar = temp_EVs(i).E_tar + deltaE;
        end
    end
end
EVs = temp_EVs;  % 更新原始结构体数组

%% 5. 预计算模块
%% 5.1 EV基线功率计算
H = 24;                                          % 预测时域（小时）
H_steps = H / dt;                                % 转换为时间步数
temp_EVs = EVs;  % 创建临时结构体数组
parfor i = 1:num_EV
    [~, temp_EVs(i).P_base_sequence] = EVbaseP_single_longstep(...
        temp_EVs(i).C_EV, temp_EVs(i).eta, temp_EVs(i).E_tar, temp_EVs(i).E_in,...
        temp_EVs(i).t_dep, temp_EVs(i).t_in, dt, temp_EVs(i).r, temp_EVs(i).p_on,...
        temp_EVs(i).SOC, H_steps+1);
end
EVs = temp_EVs;  % 更新原始结构体数组

temp_ACs = ACs;  % 创建临时结构体数组
parfor i = 1:num_AC
    [alpha, beta, gamma] = calculateACABC_single(...
        temp_ACs(i).R, temp_ACs(i).C, temp_ACs(i).eta,...
        temp_ACs(i).Tmax, temp_ACs(i).Tmin, temp_ACs(i).Tset, dt);
    temp_ACs(i).alpha = alpha;
    temp_ACs(i).beta = beta;
    temp_ACs(i).gamma = gamma;
end
ACs = temp_ACs;  % 更新原始结构体数组

%% 6. 主时间循环
%% 6.1 结果预分配
AC_Up = zeros(length(time_points),1);
AC_Down = zeros(length(time_points),1);
EV_Up = zeros(length(time_points),1);
EV_Down = zeros(length(time_points),1);
SOC_AC = zeros(num_AC, length(time_points));
SOC_EV = zeros(num_EV, length(time_points));
m3 = zeros(num_EV,1);

% 新增：存储每台设备的调节能力
AC_Up_Individual = zeros(num_AC, length(time_points));
AC_Down_Individual = zeros(num_AC, length(time_points));
EV_Up_Individual = zeros(num_EV, length(time_points));
EV_Down_Individual = zeros(num_EV, length(time_points));

%% 6.2 时间步进循环
for t_idx = 1:length(time_points)
    t = time_points(t_idx);
    fprintf('\n== 时间 %.1f小时 ==\n', t);
    
    %% 电动汽车集群处理 - 使用临时结构体数组
    temp_EVs = EVs;  % 创建临时结构体数组
    temp_EV_Up = 0;
    temp_EV_Down = 0;
    temp_SOC_EV = zeros(num_EV,1);
    temp_m3 = m3;
    
    parfor i = 1:num_EV
        online = (t >= temp_EVs(i).t_in) && (t <= temp_EVs(i).t_dep);
        if online
            %% 基线功率获取
            temp_EVs(i).P_base = temp_EVs(i).P_base_sequence(t_idx);
            
            %% 状态更新
            [~, ~, m3_val] = calculateEVABC_single(...
                temp_EVs(i).C_EV, temp_EVs(i).eta, temp_EVs(i).E_tar,...
                temp_EVs(i).E_in, temp_EVs(i).t_dep, temp_EVs(i).t_in, dt, temp_EVs(i).r);
            temp_m3(i) = m3_val;
            [temp_EVs(i).E_exp, temp_EVs(i).E_current, temp_EVs(i).P_current, temp_EVs(i).SOC] = ...
                calculateEVS_single(m3_val, temp_EVs(i).E_exp, temp_EVs(i).eta,...
                temp_EVs(i).E_current, temp_EVs(i).P_current, temp_EVs(i).C_EV,...
                temp_EVs(i).r, temp_EVs(i).p_on, dt);
            
            %% 调节潜力计算
            [DeltaP_plus, DeltaP_minus] = calculateEVAdjustmentPotentia(...
                temp_EVs(i).C_EV, temp_EVs(i).r, temp_EVs(i).eta, temp_EVs(i).E_tar,...
                temp_EVs(i).E_in, temp_EVs(i).E_current, temp_EVs(i).t_dep,...
                temp_EVs(i).t_in, temp_EVs(i).p_on, 0, temp_EVs(i).P_base,...
                temp_EVs(i).SOC, dt);
            
            temp_EV_Up = temp_EV_Up + DeltaP_plus;
            temp_EV_Down = temp_EV_Down + DeltaP_minus;
            
            % 新增：存储单台EV调节能力
            EV_Up_Individual(i, t_idx) = DeltaP_plus;
            EV_Down_Individual(i, t_idx) = DeltaP_minus;
            
            temp_SOC_EV(i) = temp_EVs(i).SOC;
        else
            temp_SOC_EV(i) = temp_EVs(i).SOC;
        end
    end
    
    EVs = temp_EVs;  % 更新原始结构体数组
    m3 = temp_m3;
    EV_Up(t_idx) = temp_EV_Up;
    EV_Down(t_idx) = temp_EV_Down;
    SOC_EV(:, t_idx) = temp_SOC_EV;
    
    %% 空调分析 - 使用临时结构体数组
    temp_ACs = ACs;  % 创建临时结构体数组
    temp_AC_Up = 0;
    temp_AC_Down = 0;
    temp_SOC_AC = zeros(num_AC,1);
    
    parfor i = 1:num_AC
        if temp_ACs(i).ptcp
            % 计算基线功率
            temp_ACs(i).P_base = ACbaseP_single(...
                temp_ACs(i).T_ja(t_idx), temp_ACs(i).Tset, temp_ACs(i).R, temp_ACs(i).eta);
            
            % 更新SOC
            temp_ACs(i).SOC = calculateACS_single(...
                temp_ACs(i).T_ja(t_idx), temp_ACs(i).Tmax, temp_ACs(i).Tmin);
            
            % 计算调节潜力
            [DeltaP_plus, DeltaP_minus] = calculateACAdjustmentPotentia(...
                temp_ACs(i).P_base, 2*abs(temp_ACs(i).P_base), 0,...
                temp_ACs(i).alpha, temp_ACs(i).beta, temp_ACs(i).gamma,...
                temp_ACs(i).SOC, dt);
            
            temp_AC_Up = temp_AC_Up + DeltaP_plus;
            temp_AC_Down = temp_AC_Down + DeltaP_minus;
            
            % 新增：存储单台AC调节能力
            AC_Up_Individual(i, t_idx) = DeltaP_plus;
            AC_Down_Individual(i, t_idx) = DeltaP_minus;
            
            temp_SOC_AC(i) = temp_ACs(i).SOC;
        end
    end
    
    ACs = temp_ACs;  % 更新原始结构体数组
    AC_Up(t_idx) = temp_AC_Up;
    AC_Down(t_idx) = temp_AC_Down;
    SOC_AC(:, t_idx) = temp_SOC_AC;
end

%% 将结果组装回结构体
results = struct(...
    'AC_Up', AC_Up,...
    'AC_Down', AC_Down,...
    'EV_Up', EV_Up,...
    'EV_Down', EV_Down,...
    'SOC_AC', SOC_AC,...
    'SOC_EV', SOC_EV,...
    'm3', m3,...
    'AC_Up_Individual', AC_Up_Individual,...  % 新增：每台AC上调能力
    'AC_Down_Individual', AC_Down_Individual,... % 新增：每台AC下调能力
    'EV_Up_Individual', EV_Up_Individual,... % 新增：每台EV上调能力
    'EV_Down_Individual', EV_Down_Individual); % 新增：每台EV下调能力
