%% 全功能虚拟电厂调节潜力分析系统（优化版）
clear; close all; clc;

%% 1. 系统初始化
rng(2023);                                      % 固定随机种子，保证结果可重复
T_total = 24;                                   % 总时长（小时）
dt = 0.05;                                      % 时间分辨率（小时）
time_points = 0:dt:T_total;                    % 生成时间序列
base_price = 30;                                % 基础电价（元/kWh）

%% 2. 初始化参数
generateResidentialEV(100);                    % 生成居民区电动汽车数据
generateWorkplaceEV(100);                      % 生成工作区电动汽车数据
evFile = 'EV_居民区.xlsx';                 % 电动汽车数据文件

%% 3. 读取设备参数
EVs = initializeEVsFromExcel(evFile);          % 从Excel导入电动汽车参数
num_EV = length(EVs);                          % 获取电动汽车总数

%% 4. 参数预提取与预处理（优化关键点）
% 提取结构体参数到数组
[C_EV, eta, E_tar, E_in, t_dep, t_in, r, p_on, SOC] = deal(...
    [EVs.C_EV], [EVs.eta], [EVs.E_tar], [EVs.E_in],...
    [EVs.t_dep], [EVs.t_in], [EVs.r], [EVs.p_on], [EVs.SOC]);

% 预计算在线状态矩阵（关键优化）
online_status = false(num_EV, length(time_points));
parfor i = 1:num_EV
    online_status(i,:) = (time_points >= t_in(i)) & (time_points <= t_dep(i));
end
% 初始化状态数组（替代结构体字段访问）
[E_exp, E_current, P_current] = deal([EVs.E_exp], [EVs.E_current], [EVs.P_current]);
%% 5. 预计算模块
%% 5.1 EV基线功率计算
H = 24;                                          % 预测时域（小时）
H_steps = H / dt;                                % 转换为时间步数
parfor i = 1:num_EV
    [~, EVs(i).P_base_sequence] = EVbaseP_single_longstep(...
        EVs(i).C_EV, EVs(i).eta, EVs(i).E_tar, EVs(i).E_in,...
        EVs(i).t_dep, EVs(i).t_in, dt, EVs(i).r, EVs(i).p_on,...
        EVs(i).SOC, H_steps+1);
end
% 预提取基线功率序列（优化内存访问）
P_base_sequence = zeros(num_EV, length(time_points));
for i = 1:num_EV
    P_base_sequence(i,:) = EVs(i).P_base_sequence;
end

%% 7.1 结果预分配优化
results = struct(...
    'EV_Up', zeros(length(time_points), 1),...
    'EV_Down', zeros(length(time_points), 1),...
    'SOC_EV', zeros(num_EV, length(time_points)));

%% 7.2 时间步进循环优化
for t_idx = 1:length(time_points)
    %% 并行处理电动汽车集群（关键优化）
    DeltaP_plus = 0;
    DeltaP_minus = 0;
    active_EVs = find(online_status(:,t_idx)); % 预查找在线EV索引
    
    % 批量处理在线EV（减少循环次数）
    for i = active_EVs'
        %% 状态更新（使用预存数组）
        [~, ~, m3] = calculateEVABC_single(...
            C_EV(i), eta(i), E_tar(i), E_in(i),...
            t_dep(i), t_in(i), dt, r(i));
        
        [E_exp_i, E_current_i, P_current_i, SOC_i] = ...
            calculateEVS_single(m3, E_exp(i), eta(i),...
            E_current(i), P_current(i), C_EV(i),...
            r(i), p_on(i), dt);
        
        %% 调节潜力计算（使用数组参数）
        [dPlus, dMinus] = calculateEVAdjustmentPotentia(...
            C_EV(i), r(i), eta(i), E_tar(i),...
            E_in(i), E_current_i, t_dep(i),...
            t_in(i), p_on(i), 0,...
            P_base_sequence(i,t_idx), SOC_i, dt);
        
        %% 结果累加
        DeltaP_plus = DeltaP_plus + dPlus;
        DeltaP_minus = DeltaP_minus + dMinus;
        
        %% 状态更新回数组
        E_exp(i) = E_exp_i;
        E_current(i) = E_current_i;
        P_current(i) = P_current_i;
        results.SOC_EV(i, t_idx) = SOC_i;
    end
    
    %% 存储调节潜力结果
    results.EV_Up(t_idx) = DeltaP_plus;
    results.EV_Down(t_idx) = DeltaP_minus;
    
    %% 状态回写结构体（每5%进度更新一次）
    if mod(t_idx, round(length(time_points)/20)) == 0
        for i = active_EVs'
            EVs(i).E_exp = E_exp(i);
            EVs(i).E_current = E_current(i);
            EVs(i).P_current = P_current(i);
            EVs(i).SOC = results.SOC_EV(i, t_idx);
        end
    end
end

%% 8. 最终状态回写
for i = 1:num_EV
    EVs(i).E_exp = E_exp(i);
    EVs(i).E_current = E_current(i);
    EVs(i).P_current = P_current(i);
    EVs(i).SOC = results.SOC_EV(i, end);
end
%% 7. 可视化与结果分析
%% 7.1 调节潜力曲线可视化
figure('Position', [100 100 1200 800])
subplot(3, 1, 2)
plot(time_points, results.EV_Up, 'r', time_points, results.EV_Down, 'b')
title('电动汽车集群调节潜力'), grid on
legend('上调潜力', '下调潜力', 'Location', 'best')

%% 7.2 SOC时空分布可视化
figure('Position', [100 100 1400 600])
subplot(1, 2, 2)
imagesc(time_points, 1:num_EV, results.SOC_EV)
title('电动汽车SOC时空分布')
xlabel('时间 (小时)'), ylabel('设备编号')
colorbar
