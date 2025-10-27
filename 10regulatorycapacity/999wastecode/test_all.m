%% 全功能虚拟电厂调节潜力分析系统
clear; close all; clc;

%% 1. 系统初始化
rng(2023);                      % 固定随机种子
num_AC    = 20;                 % 空调数量
num_EV    = 10;                 % 电动汽车数量
T_total   = 24;                 % 总时长(小时)
dt        = 0.5;                % 时间分辨率
time_points = 0:dt:T_total;
base_price = 30;                % 基础电价(元/kWh)

%% 2. 设备参数初始化
% 2.1 空调集群初始化
ACs = initializeACs(num_AC, rng);

% 添加激励响应参数
for i = 1:num_AC
    ACs(i).Tset_original = ACs(i).Tset;  % 保存原始设定温度
    ACs(i).Tmax_original = 26;           % 初始温度上限
    ACs(i).Tmin_original = 18;           % 初始温度下限
end

% 2.2 电动汽车集群初始化
EVs = initializeEVs(num_EV, rng);

% 添加灵活性参数
for i = 1:num_EV
    EVs(i).E_tar_original = EVs(i).E_tar;  % 保存原始目标电量
end

%% 3. 激励响应模块
% incentive_price = 0.6 + 0.4*sin(2*pi*time_points/24);  % 分时激励电价

% 3.1 AC温度设定调整
% T_set_max = 3;  % 最大温度偏移
p_min = 10; p_max = 50; 
p_min_prime = 15; p_max_prime = 45;
T_set_max = 3;  % 最大温度偏移

for i = 1:num_AC
    % 计算参与度
    participation = calculateParticipation(ACs(i).p_incentive, base_price);
    [~, ~, deltaT] = incentiveTempAC(...
        ACs(i).p_incentive, p_min, p_max, p_min_prime, p_max_prime, T_set_max);
    ACs(i).ptcp = (rand() < participation);
    
    % 80%概率响应温度变化
    if ACs(i).ptcp
        ACs(i).Tmax = ACs(i).Tset + deltaT;
        ACs(i).Tmin = ACs(i).Tset - deltaT;
    end
    
    % 基础温度波动（含日夜变化）
    base_temp = ACs(i).Tset + 4*sin(2*pi*time_points/24); 
    
    % 添加随机扰动（幅度为温度范围的20%）
    temp_range = ACs(i).Tmax - ACs(i).Tmin;
    noise = 0.2 * temp_range * randn(size(time_points));
    
    % 温度限幅（严格控制在Tmin-Tmax之间）
    ACs(i).T_ja = min(max(base_temp + noise, ACs(i).Tmin), ACs(i).Tmax);
end

% 3.2 EV目标电量调整
E_tar_max = 0.2 * [EVs.C_EV];  % 计算最大目标电量变化

for i = 1:num_EV
    participation = calculateParticipation(EVs(i).p_incentive, base_price);
    [~, ~, deltaE] = incentiveTempEV(EVs(i).p_incentive, 5, 50, 10, 40, E_tar_max(i));
    ptcp_result = (rand() < participation); 
    EVs(i).ptcp = ptcp_result;
    
    if ptcp_result
        EVs(i).E_tar = EVs(i).E_tar - deltaE;
        if EVs(i).E_tar <= EVs(i).E_in
            EVs(i).E_tar = EVs(i).E_tar + deltaE;
        end
    end
end

%% 4. 预计算模块
% 4.1 环境温度生成
% 4.2 EV基线功率预计算
for i = 1:num_EV
    [m1, m2, m3] = calculateEVABC_single(EVs(i).C_EV, EVs(i).eta,...
        EVs(i).E_tar, EVs(i).E_in, EVs(i).t_dep, EVs(i).t_in, dt, EVs(i).r);
    [~, EVs(i).P_base_sequence] = EVbaseP_single_longstep(m1, m2, m3,...
        EVs(i).p_on, EVs(i).SOC, length(time_points));
end

%% 5. 主时间循环
% 预分配存储矩阵
results = struct(...
    'AC_Up',     zeros(length(time_points),1),...
    'AC_Down',   zeros(length(time_points),1),...
    'EV_Up',     zeros(length(time_points),1),...
    'EV_Down',   zeros(length(time_points),1),...
    'SOC_AC',    zeros(num_AC, length(time_points)),...
    'SOC_EV',    zeros(num_EV, length(time_points)));

for t_idx = 1:length(time_points)
    t = time_points(t_idx);
    fprintf('\n== 时间 %.1f小时 ==\n', t);
    
    %% 5.1 空调集群处理
    for i = 1:num_AC
        % 基线功率计算
        if ACs(i).ptcp
            ACs(i).P_base = ACbaseP_single(ACs(i).T_ja(t_idx), ACs(i).Tset,...
                ACs(i).R, ACs(i).eta);
            
            % 状态系数更新
            [alpha, beta, gamma] = calculateACABC_single(ACs(i).R, ACs(i).C,...
                ACs(i).eta, ACs(i).Tmax, ACs(i).Tmin, ACs(i).Tset, dt);
            ACs(i).alpha = alpha;
            ACs(i).beta = beta;
            ACs(i).gamma = gamma;
            
            % SOC更新（含不确定性因子）
            ACs(i).SOC = calculateACS_single(ACs(i).T_ja(t_idx), ACs(i).Tmax, ACs(i).Tmin);
            
            % 调节潜力计算（含设备老化因子）
            [DeltaP_plus, DeltaP_minus] = calculateACAdjustmentPotentia(...
                ACs(i).P_base, 2*abs(ACs(i).P_base), 0, alpha, beta, gamma,...
                ACs(i).SOC, dt * (0.95 + 0.1*rand()));
            
            results.AC_Up(t_idx)   = results.AC_Up(t_idx) + DeltaP_plus;
            results.AC_Down(t_idx) = results.AC_Down(t_idx) + DeltaP_minus;
            results.SOC_AC(i, t_idx) = ACs(i).SOC;
        end
    end
    
    %% 5.2 电动汽车集群处理
    for i = 1:num_EV
        online = (t >= EVs(i).t_in) && (t <= EVs(i).t_dep);
        if online
            % 获取基线功率
            EVs(i).P_base = EVs(i).P_base_sequence(t_idx);
            
            % 状态更新
            [~, ~, m3] = calculateEVABC_single(EVs(i).C_EV, EVs(i).eta,...
                EVs(i).E_tar, EVs(i).E_in, EVs(i).t_dep, EVs(i).t_in, dt, EVs(i).r);
            [EVs(i).E_exp, EVs(i).E_current, EVs(i).P_current, EVs(i).SOC] = ...
                calculateEVS_single(m3, EVs(i).E_exp, EVs(i).eta,...
                EVs(i).E_current, EVs(i).P_current, EVs(i).C_EV,...
                EVs(i).r, EVs(i).p_on, dt);
            
            % 调节潜力计算（含电池衰减因子）
            [DeltaP_plus, DeltaP_minus] = calculateEVAdjustmentPotentia(...
                EVs(i).C_EV, EVs(i).r,EVs(i).eta, EVs(i).E_tar, EVs(i).E_in,EVs(i).E_current,...
                EVs(i).t_dep, EVs(i).t_in, EVs(i).p_on, 0,...
                EVs(i).P_base, EVs(i).SOC, dt);
            
            results.EV_Up(t_idx)   = results.EV_Up(t_idx) + DeltaP_plus;
            results.EV_Down(t_idx) = results.EV_Down(t_idx) + DeltaP_minus;
            results.SOC_EV(i, t_idx) = EVs(i).SOC;
        else
            results.SOC_EV(i, t_idx) = EVs(i).SOC;
        end
    end
    
    %% 5.3 协同调度优化
    % 当前时刻可用调节能力
    P_avail_up   = results.AC_Up(t_idx) + results.EV_Up(t_idx);
    P_avail_down = results.AC_Down(t_idx) + results.EV_Down(t_idx);
    
    % 电网需求（示例：跟踪可再生能源波动）
    P_renewable = 1000 * sin(2*pi*t/24);  % 模拟可再生能源出力
    P_demand    = 1500 + 500*randn();     % 模拟实时负荷需求
    delta_P     = P_demand - P_renewable; % 需要调节的功率
    
    % 优化调度
    if delta_P > 0
        % 需要上调
        dispatched = min(P_avail_up, delta_P);
        ratio_AC = results.AC_Up(t_idx) / P_avail_up;
        ratio_EV = 1 - ratio_AC;
    else
        % 需要下调
        dispatched = max(P_avail_down, delta_P);
        ratio_AC = results.AC_Down(t_idx) / abs(P_avail_down);
        ratio_EV = 1 - ratio_AC;
    end
    
    % 记录调度结果
    results.Dispatched(t_idx) = dispatched;
    results.Ratio_AC(t_idx)   = ratio_AC;
    results.Ratio_EV(t_idx)   = ratio_EV;
end

%% 6. 可视化与结果分析
% 6.1 调节潜力曲线
figure('Position', [100 100 1200 800])
subplot(3,1,1)
plot(time_points, results.AC_Up, 'r', time_points, results.AC_Down, 'b')
title('空调集群调节潜力')
legend('上调潜力', '下调潜力'), grid on

subplot(3,1,2)
plot(time_points, results.EV_Up, 'r', time_points, results.EV_Down, 'b')
title('电动汽车集群调节潜力')
legend('上调潜力', '下调潜力'), grid on

subplot(3,1,3)
total_up   = results.AC_Up + results.EV_Up;
total_down = results.AC_Down + results.EV_Down;
plot(time_points, total_up, 'r', time_points, total_down, 'b')
title('虚拟电厂总调节潜力')
legend('总上调潜力', '总下调潜力'), grid on

% 6.2 SOC时空分布
figure('Position', [100 100 1400 600])
subplot(1,2,1)
imagesc(time_points, 1:num_AC, results.SOC_AC)
title('空调SOC时空分布'), xlabel('时间 (小时)'), ylabel('设备编号')
colorbar

subplot(1,2,2)
imagesc(time_points, 1:num_EV, results.SOC_EV)
title('电动汽车SOC时空分布'), xlabel('时间 (小时)'), ylabel('设备编号')
colorbar

% 6.3 调度效果分析
figure('Position', [100 100 1200 800])
subplot(2,1,1)
plot(time_points, results.Dispatched, 'k', 'LineWidth', 2)
hold on
plot(time_points, results.Ratio_AC*100, 'b--')
plot(time_points, results.Ratio_EV*100, 'r--')
title('调度执行情况')
legend('实际调度量', '空调调度占比', '电动汽车调度占比'), grid on

subplot(2,1,2)
histogram(results.Ratio_AC, 20)
title('空调调度占比分布'), xlabel('占比'), ylabel('出现频次')
