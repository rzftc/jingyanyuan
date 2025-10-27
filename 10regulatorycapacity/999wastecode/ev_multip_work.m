%% 全功能虚拟电厂调节潜力分析系统
clear; close all; clc;

%% 1. 系统初始化
rng(2023);                                      % 固定随机种子
T_total = 24;                                   % 总时长(小时)
dt = 0.05;                                       % 时间分辨率
time_points = 0:dt:T_total;
base_price = 30;                                % 基础电价(元/kWh)

%% 2. 初始化参数
generateResidentialEV(100);
generateWorkplaceEV(100);
evFile = 'EV_Residential.xlsx';

%% 3. 读取设备参数
EVs = initializeEVsFromExcel(evFile);
num_EV = length(EVs);                           % 电动汽车数量

%% 4. 激励响应模块
p_min = 15; p_max = 50; 
p_min_prime = 10; p_max_prime = 45;

% 4.2 EV目标电量调整
E_tar_max = 0.2 * [EVs.C_EV];                   % 计算最大目标电量变化
for i = 1:num_EV
    participation = calculateParticipation(EVs(i).p_incentive, base_price);
    [~, ~, deltaE] = incentiveTempEV(EVs(i).p_incentive, 5, 50, 10, 40, E_tar_max(i));
    ptcp_result = (rand() < participation); 
    EVs(i).ptcp = ptcp_result;
    
    if EVs(i).ptcp
        EVs(i).E_tar = EVs(i).E_tar - deltaE;
        % 防止目标电量低于初始电量
        if EVs(i).E_tar <= EVs(i).E_in
            EVs(i).E_tar = EVs(i).E_tar + deltaE;
        end
    end
end

%% 5. 预计算模块
% 5.1 EV基线功率计算
H = 24;  % 预测时域24小时（需转换为时间步数）
H_steps = H / dt;   % 实际时间步数
for i = 1:num_EV
    [~, EVs(i).P_base_sequence] = EVbaseP_single_longstep(EVs(i).C_EV, EVs(i).eta,...
        EVs(i).E_tar, EVs(i).E_in, EVs(i).t_dep, EVs(i).t_in, dt, EVs(i).r,...
        EVs(i).p_on, EVs(i).SOC, H_steps);
end

%% 6. 主时间循环
% 预分配存储矩阵
results = struct(...
    'EV_Up', zeros(length(time_points),1),...
    'EV_Down', zeros(length(time_points),1),...
    'SOC_EV', zeros(num_EV, length(time_points)));

for t_idx = 1:length(time_points)
    t = time_points(t_idx);
    fprintf('\n== 时间 %.1f小时 ==\n', t);
    
    %% 6.2 电动汽车集群处理
    for i = 1:num_EV
        online = (t >= EVs(i).t_in) && (t <= EVs(i).t_dep);
        if online
            % 获取基线功率
            EVs(i).P_base = EVs(i).P_base_sequence(t_idx);
            
            % 状态更新
            [~, ~, m3] = calculateEVABC_single(EVs(i).C_EV, EVs(i).eta,...
                EVs(i).E_tar, EVs(i).E_in, EVs(i).t_dep, EVs(i).t_in, dt, EVs(i).r);
            results.m3(i)=m3;
            [EVs(i).E_exp, EVs(i).E_current, EVs(i).P_current, EVs(i).SOC] = ...
                calculateEVS_single(m3, EVs(i).E_exp, EVs(i).eta,...
                EVs(i).E_current, EVs(i).P_current, EVs(i).C_EV,...
                EVs(i).r, EVs(i).p_on, dt);
            
            % 调节潜力计算（含电池衰减因子）
            [DeltaP_plus, DeltaP_minus] = calculateEVAdjustmentPotentia(...
                EVs(i).C_EV, EVs(i).r,EVs(i).eta, EVs(i).E_tar, EVs(i).E_in,EVs(i).E_current,...
                EVs(i).t_dep, EVs(i).t_in, EVs(i).p_on, 0,...
                EVs(i).P_base, EVs(i).SOC, dt);
            
            results.EV_Up(t_idx) = results.EV_Up(t_idx) + DeltaP_plus;
            results.EV_Down(t_idx) = results.EV_Down(t_idx) + DeltaP_minus;
            results.SOC_EV(i, t_idx) = EVs(i).SOC;
        else
            results.SOC_EV(i, t_idx) = EVs(i).SOC;
        end
    end
end
   
%% 7. 可视化与结果分析
% 7.1 调节潜力曲线
figure('Position', [100 100 1200 800])
subplot(3,1,2)
plot(time_points, results.EV_Up, 'r', time_points, results.EV_Down, 'b')
title('电动汽车集群调节潜力')
legend('上调潜力', '下调潜力'), grid on
% 7.2 SOC时空分布
figure('Position', [100 100 1400 600])

subplot(1,2,2)
imagesc(time_points, 1:num_EV, results.SOC_EV)
title('电动汽车SOC时空分布'), xlabel('时间 (小时)'), ylabel('设备编号')
colorbar

