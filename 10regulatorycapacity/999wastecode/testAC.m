%% 基于全流程的AC集群测试（结构体数组版本）
clear; close all; clc;

%% 1. 参数初始化（使用initializeACs）
num_AC = 10;     % AC数量
rngSeed = 2023;  % 统一随机种子
ACs = initializeACs(num_AC, rngSeed); % 调用初始化函数

% 电价参数
base_price = 30;               
T_total = 24;    % 总时间（小时）
dt = 0.5;        % 时间步长（小时）
p_incentive = linspace(0, 60, T_total/dt+1); 
time_points = 0:dt:T_total;

% 存储历史数据
history = struct();
history.time = 0:dt:T_total;
history.SOC = zeros(num_AC, length(history.time));
history.P_adj_plus = zeros(num_AC, length(history.time));
history.P_adj_minus = zeros(num_AC, length(history.time));

%% 阶段2: 更新温度设定
p_min = 10; p_max = 50; 
p_min_prime = 15; p_max_prime = 45;
T_set_max = 3; % 最大温度偏移

% 为每个AC单独计算温度变化
for i = 1:num_AC
    participation = calculateParticipation(ACs(i).p_incentive, base_price);
    [~, ~, deltaT] = incentiveTempAC(...
        ACs(i).p_incentive, p_min, p_max, p_min_prime, p_max_prime, T_set_max);
    ACs(i).ptcp=(rand() < participation);
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
%% 2. 主时间循环（保持核心逻辑）
for t_idx = 1:length(history.time)
    t = history.time(t_idx);
    fprintf('\n=== 时间 %.1f小时 ===\n', t);
    
    %% 阶段3: 计算基准功率（保持原有逻辑）
    for i = 1:num_AC
        if ACs(i).ptcp
            ACs(i).P_base = ACbaseP_single(ACs(i).T_ja(t_idx), ACs(i).Tset, ACs(i).R, ACs(i).eta);
        
        %% 阶段4: 更新状态参数
            [alpha, beta, gamma] = calculateACABC_single(... % 保持原有参数
                ACs(i).R, ACs(i).C, ACs(i).eta, ACs(i).Tmax, ACs(i).Tmin, ACs(i).Tset, dt);
            ACs(i).alpha = alpha;
            ACs(i).beta = beta;
            ACs(i).gamma = gamma;
        
        %% 阶段5: SOC更新与归一化
            ACs(i).SOC = calculateACS_single(ACs(i).T_ja(t_idx),ACs(i).Tmax, ACs(i).Tmin);
        
        %% 阶段6: 计算调节潜力（保持原有逻辑）
            [DeltaP_plus, DeltaP_minus] = calculateACAdjustmentPotentia(... % 参数不变
                ACs(i).P_base, abs(ACs(i).P_base)*2, 0, ACs(i).alpha, ACs(i).beta, ACs(i).gamma, ACs(i).SOC, dt);
            history.P_adj_plus(i, t_idx) = DeltaP_plus;
            history.P_adj_minus(i, t_idx) = DeltaP_minus;
        end
    end
    
    %% 存储SOC数据
    history.SOC(:, t_idx) = [ACs.SOC];
end

%% 3. 可视化（使用黑色不同线型，加大坐标轴名称字号）
% 定义线型数组

lineStyles = {'-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--'}; % 为10条曲线准备不同线型
% 图4：总调节能力
% 计算集群总潜力
total_plus = sum(history.P_adj_plus, 1);
total_minus = sum(history.P_adj_minus, 1);
% 处理全零情况
auto_ylim = @(data) [min(data)*1.1, max(data)*1.1]; % 自动适应正负范围
% 图4：总上调潜力独立展示
figure('Position', [100, 100, 600, 400]);
plot(history.time, total_plus, 'r-', 'LineWidth', 2);
title('集群总上调调节潜力', 'FontSize', 14);
xlabel('时间(小时)', 'FontSize', 12);
ylabel('调节能力 (kW)', 'FontSize', 12);
grid on;
ylim(auto_ylim(total_plus)); % 使用安全范围函数
print('-dpng', '-r600', 'Total_Up_Adjustment.png');

% 图5：总下调潜力独立展示
figure('Position', [100, 100, 600, 400]);
plot(history.time, total_minus, 'b-', 'LineWidth', 2);
title('集群总下调调节潜力', 'FontSize', 14);
xlabel('时间(小时)', 'FontSize', 12);
ylabel('调节能力 (kW)', 'FontSize', 12);
grid on;
ylim(auto_ylim(total_minus));
print('-dpng', '-r600', 'Total_Down_Adjustment.png');
% 图1：SOC演化
figure('Position', [100, 100, 600, 400]);
hold on;
for i = 1:num_AC
    plot(history.time, history.SOC(i,:), 'Color', 'k', 'LineStyle', lineStyles{i}, 'LineWidth', 1);
end
hold off;
title('归一化SOC演化', 'FontSize', 14);
xlabel('时间(小时)', 'FontSize', 14);
ylabel('SOC', 'FontSize', 14);
ylim([0 1]);
grid on;
print('-dpng', '-r600', 'SOC_Evolution.png');

% 图2：上调潜力
figure('Position', [100, 100, 600, 400]);
hold on;
for i = 1:num_AC
    plot(history.time, history.P_adj_plus(i,:), 'Color', 'k', 'LineStyle', lineStyles{i}, 'LineWidth', 1);
end
hold off;
title('上调潜力', 'FontSize', 14);
xlabel('时间(小时)', 'FontSize', 14);
ylabel('ΔP+ (kW)', 'FontSize', 14);
grid on;
print('-dpng', '-r600', 'Power_Adjustment_Plus.png');

% 图3：下调潜力
figure('Position', [100, 100, 600, 400]);
hold on;
for i = 1:num_AC
    plot(history.time, history.P_adj_minus(i,:), 'Color', 'k', 'LineStyle', lineStyles{i}, 'LineWidth', 1);
end
hold off;
title('下调潜力', 'FontSize', 14);
xlabel('时间(小时)', 'FontSize', 14);
ylabel('ΔP- (kW)', 'FontSize', 14);
grid on;
print('-dpng', '-r600', 'Power_Adjustment_Minus.png');