%% plot_reliable_regulation_domain.m
% 功能：读取可靠调节域数据文件并绘制可视化图表

clear; close all; clc;

%% 1. 加载数据
data_file = 'reliable_regulation_domain_1000_mix_01.mat';

if ~exist(data_file, 'file')
    error(['数据文件 %s 不存在。\n' ...
           '请先运行 main_scenario_generation_diff_mix.m 生成数据，' ...
           '或确认文件路径正确。'], data_file);
end

fprintf('正在加载数据: %s ...\n', data_file);
load(data_file);
fprintf('加载完成。包含 %d 个场景，时间步长 %d 个。\n', num_scenarios, length(time_points));

%% 2. 准备绘图参数
% 从加载的 time_points 自动推导仿真起止时间
simulation_start_hour = time_points(1);
simulation_end_hour   = time_points(end);

% 定义坐标轴刻度 (假设是 6:00 到次日 6:00 的 24 小时制)
% 如果您的仿真时间段不同，可以手动调整这里
x_ticks = [6, 12, 18, 24, 30];
x_labels = {'06:00', '12:00', '18:00', '00:00 (次日)', '06:00 (次日)'};

%% 3. 执行绘图
figure('Name', '可靠调节域提取结果', 'Position', [100, 100, 1200, 800], 'Color', 'w');

% --- 子图 1: 空调 (AC) ---
subplot(2, 1, 1); 
hold on;

% 绘制所有场景的背景（半透明）
% 注意：如果场景数非常多(>1000)，绘制所有线条可能会卡顿，可以考虑只画前100个
plot_limit = min(num_scenarios, 1000); 
plot(time_points, Scenarios_AC_Up(:, 1:plot_limit), 'Color', [0.6, 0.8, 1, 0.15], 'HandleVisibility', 'off');
plot(time_points, Scenarios_AC_Down(:, 1:plot_limit), 'Color', [1, 0.6, 0.6, 0.15], 'HandleVisibility', 'off');

% 绘制可靠边界
p1 = plot(time_points, Reliable_AC_Up, 'b-', 'LineWidth', 2, 'DisplayName', '可靠上调边界');
p2 = plot(time_points, Reliable_AC_Down, 'r-', 'LineWidth', 2, 'DisplayName', '可靠下调边界');

% 辅助线与装饰
yline(0, 'k--'); 
title('空调 (AC) 聚合体可靠调节域', 'FontSize', 14); 
legend([p1, p2], 'Location', 'best', 'FontSize', 10); 
grid on;

% 设置坐标轴
xlim([simulation_start_hour, simulation_end_hour]);
set(gca, 'XTick', x_ticks, 'XTickLabel', x_labels, 'FontSize', 12);
xlabel('时间'); ylabel('功率 (kW)');


% --- 子图 2: 电动汽车 (EV) ---
subplot(2, 1, 2); 
hold on;

% 绘制所有场景的背景（半透明）
plot(time_points, Scenarios_EV_Up(:, 1:plot_limit), 'Color', [0.6, 0.8, 1, 0.15], 'HandleVisibility', 'off');
plot(time_points, Scenarios_EV_Down(:, 1:plot_limit), 'Color', [1, 0.6, 0.6, 0.15], 'HandleVisibility', 'off');

% 绘制可靠边界
p3 = plot(time_points, Reliable_EV_Up, 'b-', 'LineWidth', 2, 'DisplayName', '可靠上调边界');
p4 = plot(time_points, Reliable_EV_Down, 'r-', 'LineWidth', 2, 'DisplayName', '可靠下调边界');

% 辅助线与装饰
yline(0, 'k--'); 
title('电动汽车 (EV) 聚合体可靠调节域', 'FontSize', 14); 
legend([p3, p4], 'Location', 'best', 'FontSize', 10); 
grid on;

% 设置坐标轴
xlim([simulation_start_hour, simulation_end_hour]);
set(gca, 'XTick', x_ticks, 'XTickLabel', x_labels, 'FontSize', 12);
xlabel('时间'); ylabel('功率 (kW)');

fprintf('绘图完成。\n');