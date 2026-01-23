%% plot_reliable_regulation_domain.m
% 功能：读取可靠调节域数据文件，绘制无标题的高清图表并分别保存（中文文件名）
% 修正说明：
% 1. 字体大小与线宽已调整，以匹配 test_ieee30_func_8ampv_tly.m 的论文格式标准。
%    - 坐标轴刻度字号: 18
%    - 坐标轴标签字号: 20 (加粗)
%    - 图例字号: 18
%    - 坐标轴线宽: 1.5

clear; close all; clc;

%% 1. 加载数据
data_file = 'reliable_regulation_domain_soc_2bound.mat';

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

% 定义坐标轴刻度
x_ticks = [8, 14, 20, 26, 32];
x_labels = {'08:00', '14:00', '20:00', '02:00 (次日)', '08:00 (次日)'};

% 场景绘制限制（防止过卡）
plot_limit = min(num_scenarios, 1000); 

% 设置通用字体和字号 (已修改以匹配 test_ieee30_func_8ampv_tly.m)
default_font = 'Microsoft YaHei'; % 确保中文正常显示
axis_font_size = 18;   % 原为 14 -> 修改为 18
label_font_size = 20;  % 原为 16 -> 修改为 20

%% 3. 绘制并保存 AC (空调) 图像
fprintf('正在绘制 AC 调节域图像...\n');
fig_ac = figure('Name', 'AC Reliable Domain', 'Position', [100, 100, 1000, 600], 'Color', 'w');
hold on;

% 绘制所有场景的背景（半透明）
plot(time_points, 5*Scenarios_AC_Up(:, 1:plot_limit), 'Color', [0.6, 0.8, 1, 0.15], 'HandleVisibility', 'off');
plot(time_points, 5*Scenarios_AC_Down(:, 1:plot_limit), 'Color', [1, 0.6, 0.6, 0.15], 'HandleVisibility', 'off');

% 绘制可靠边界
p1 = plot(time_points, 5*Reliable_AC_Up, 'b-', 'LineWidth', 2.5, 'DisplayName', '可靠上调边界');
p2 = plot(time_points, 5*Reliable_AC_Down, 'r-', 'LineWidth', 2.5, 'DisplayName', '可靠下调边界');

% 辅助线与装饰
yline(0, 'k--', 'LineWidth', 1.2); 
% title('空调 (AC) 聚合体可靠调节域', 'FontSize', 14); % [已移除标题]

% 图例设置 (FontSize 修改为 18)
legend([p1, p2], 'Location', 'best', 'FontSize', 18, 'FontName', default_font); 
grid on;

% 设置坐标轴 (LineWidth 修改为 1.5, Label 增加 FontWeight bold)
xlim([simulation_start_hour, simulation_end_hour]);
set(gca, 'XTick', x_ticks, 'XTickLabel', x_labels, ...
    'FontSize', axis_font_size, 'FontName', default_font, 'LineWidth', 1.5);
xlabel('时间', 'FontSize', label_font_size, 'FontName', default_font); 
ylabel('功率 (kW)', 'FontSize', label_font_size, 'FontName', default_font);

% 保存为高DPI PNG (600 dpi) - 中文文件名
output_filename_ac = '空调可靠调节域.png';
print(fig_ac, output_filename_ac, '-dpng', '-r600');
fprintf('AC 图像已保存为: %s\n', output_filename_ac);


%% 4. 绘制并保存 EV (电动汽车) 图像
fprintf('正在绘制 EV 调节域图像...\n');
fig_ev = figure('Name', 'EV Reliable Domain', 'Position', [150, 150, 1000, 600], 'Color', 'w');
hold on;

% 绘制所有场景的背景（半透明）
plot(time_points, 10*Scenarios_EV_Up(:, 1:plot_limit), 'Color', [0.6, 0.8, 1, 0.15], 'HandleVisibility', 'off');
plot(time_points, 10*Scenarios_EV_Down(:, 1:plot_limit), 'Color', [1, 0.6, 0.6, 0.15], 'HandleVisibility', 'off');

% 绘制可靠边界
p3 = plot(time_points, 10*Reliable_EV_Up, 'b-', 'LineWidth', 2.5, 'DisplayName', '可靠上调边界');
p4 = plot(time_points, 10*Reliable_EV_Down, 'r-', 'LineWidth', 2.5, 'DisplayName', '可靠下调边界');

% 辅助线与装饰
yline(0, 'k--', 'LineWidth', 1.2); 
% title('电动汽车 (EV) 聚合体可靠调节域', 'FontSize', 14); % [已移除标题]

% 图例设置 (FontSize 修改为 18)
legend([p3, p4], 'Location', 'best', 'FontSize', 18, 'FontName', default_font); 
grid on;

% 设置坐标轴 (LineWidth 修改为 1.5, Label 增加 FontWeight bold)
xlim([simulation_start_hour, simulation_end_hour]);
set(gca, 'XTick', x_ticks, 'XTickLabel', x_labels, ...
    'FontSize', axis_font_size, 'FontName', default_font, 'LineWidth', 1.5);
xlabel('时间', 'FontSize', label_font_size, 'FontName', default_font); 
ylabel('功率 (kW)', 'FontSize', label_font_size, 'FontName', default_font);

% 保存为高DPI PNG (600 dpi) - 中文文件名
output_filename_ev = '电动汽车可靠调节域.png';
print(fig_ev, output_filename_ev, '-dpng', '-r600');
fprintf('EV 图像已保存为: %s\n', output_filename_ev);

fprintf('所有绘图任务完成。\n');