%% plot_reliable_100_senarios.m
% 功能：从现有数据中随机选取100个场景，处理时间轴平移（0-24h），
%       并绘制不同置信水平（100%, 95%, 80%）的可靠调节域概率边界。
%       [修改] 字体整体放大2号 (16->18, 20->22)。
%       [保留] 取消强制置零逻辑，EV数据保持缩小一半。

clear; close all; clc;

%% 1. 加载数据
data_file = 'reliable_regulation_domain_soc_residential.mat';

if ~exist(data_file, 'file')
    error(['数据文件 %s 不存在。\n' ...
           '请先运行 main_scenario_generation_soc_residential.m 生成数据，' ...
           '或确认文件路径正确。'], data_file);
end

fprintf('正在加载数据: %s ...\n', data_file);
load(data_file);

% 获取原始维度
[num_steps, num_scenarios_total] = size(Scenarios_AC_Up);

%% 2. 场景筛选 (随机100个)
target_scenarios = 100;
fprintf('正在筛选场景...\n');

if num_scenarios_total >= target_scenarios
    rng(42); % 固定随机种子，保证结果可重复
    selected_indices = randperm(num_scenarios_total, target_scenarios);
else
    selected_indices = 1:num_scenarios_total;
    fprintf('警告：可用场景不足100个，使用全部 %d 个场景。\n', num_scenarios_total);
    target_scenarios = num_scenarios_total;
end

% 提取选定场景数据
% [保留] 将EV数据缩小一半 (乘以 0.5)
S_AC_Up   = Scenarios_AC_Up(:, selected_indices);
S_AC_Down = Scenarios_AC_Down(:, selected_indices);
S_EV_Up   = Scenarios_EV_Up(:, selected_indices) * 0.8;   
S_EV_Down = Scenarios_EV_Down(:, selected_indices) * 0.8; 

fprintf('已选择 %d 个场景进行概率边界分析。\n', target_scenarios);

%% 3. 时间轴平移与数据重组 (08:00-32:00 -> 00:00-24:00)
fprintf('正在执行时间轴平移 (25-32h -> 0-8h)...\n');

% 找到分割点：24:00
idx_part1 = find(time_points >= 24); % 次日凌晨部分
idx_part2 = find(time_points < 24);  % 当日白天部分

% 构造新的 0-24h 时间轴
new_time_points = [time_points(idx_part1)-24, time_points(idx_part2)];

% 对数据矩阵按行（时间步）进行重排
process_shift = @(data) [data(idx_part1, :); data(idx_part2, :)];

AC_Up_Shifted   = process_shift(S_AC_Up);
AC_Down_Shifted = process_shift(S_AC_Down);
EV_Up_Shifted   = process_shift(S_EV_Up);
EV_Down_Shifted = process_shift(S_EV_Down);

%% 4. 计算概率边界 (Reliable Boundaries)
fprintf('正在计算 100%%, 95%%, 80%% 置信边界...\n');

% 辅助函数：计算边界
calc_bounds_up = @(data, n) [min(data, [], 2), ...                  % 100%
                             data(:, round(0.05 * n)), ...          % 95%
                             data(:, round(0.20 * n))];             % 80%

calc_bounds_down = @(data, n) [max(data, [], 2), ...                % 100%
                               data(:, round(0.95 * n)), ...        % 95%
                               data(:, round(0.80 * n))];           % 80%

% 先对每一行进行排序
AC_Up_Sorted   = sort(AC_Up_Shifted, 2);
AC_Down_Sorted = sort(AC_Down_Shifted, 2);
EV_Up_Sorted   = sort(EV_Up_Shifted, 2);
EV_Down_Sorted = sort(EV_Down_Shifted, 2);

% 提取边界 [100%, 95%, 80%]
Bounds_AC_Up   = calc_bounds_up(AC_Up_Sorted, target_scenarios);
Bounds_AC_Down = calc_bounds_down(AC_Down_Sorted, target_scenarios);
Bounds_EV_Up   = calc_bounds_up(EV_Up_Sorted, target_scenarios);
Bounds_EV_Down = calc_bounds_down(EV_Down_Sorted, target_scenarios);

%% 5. 绘图参数设置 (字体放大)
% [修改] 字体大小增加
x_ticks_new = [0, 6, 12, 18, 24];
x_labels_new = {'00:00', '06:00', '12:00', '18:00', '24:00'};
default_font = 'Microsoft YaHei';
axis_font_size = 18;  % 坐标轴刻度字体 (原16 -> 18)
label_font_size = 22; % 坐标轴标签字体 (原20 -> 22)

%% 6. 绘制 AC 调节域
fprintf('正在绘制 AC 概率边界图像...\n');
fig_ac = figure('Name', 'AC Reliable Probabilistic Domain', 'Position', [100, 100, 1000, 600], 'Color', 'w');
hold on;

% 1. 绘制背景场景
plot(new_time_points, AC_Up_Shifted, 'Color', [0.6, 0.8, 1, 0.08], 'HandleVisibility', 'off');
plot(new_time_points, AC_Down_Shifted, 'Color', [1, 0.6, 0.6, 0.08], 'HandleVisibility', 'off');

% 2. 绘制边界
% 上调
p80_u = plot(new_time_points, Bounds_AC_Up(:,3), 'b:', 'LineWidth', 1.5, 'DisplayName', '80% 置信边界');
p95_u = plot(new_time_points, Bounds_AC_Up(:,2), 'b--', 'LineWidth', 2.0, 'DisplayName', '95% 置信边界');
p100_u = plot(new_time_points, Bounds_AC_Up(:,1), 'b-', 'LineWidth', 2.5, 'DisplayName', '100% 可靠边界');

% 下调
plot(new_time_points, Bounds_AC_Down(:,3), 'r:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(new_time_points, Bounds_AC_Down(:,2), 'r--', 'LineWidth', 2.0, 'HandleVisibility', 'off');
plot(new_time_points, Bounds_AC_Down(:,1), 'r-', 'LineWidth', 2.5, 'HandleVisibility', 'off');

% 0轴线
yline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

% 装饰
grid on;
xlim([0, 24]);
% 图例设置 (字体放大)
legend([p100_u, p95_u, p80_u], 'Location', 'best', 'FontSize', 18, 'FontName', default_font);
set(gca, 'XTick', x_ticks_new, 'XTickLabel', x_labels_new, ...
    'FontSize', axis_font_size, 'FontName', default_font, 'LineWidth', 1.2);
xlabel('时间', 'FontSize', label_font_size, 'FontName', default_font); 
ylabel('AC 功率 (kW)', 'FontSize', label_font_size, 'FontName', default_font);

% 保存
save_name_ac = 'AC_可靠调节域_概率边界_100scenarios.png';
print(fig_ac, save_name_ac, '-dpng', '-r600');
fprintf('AC 图像已保存: %s\n', save_name_ac);

%% 7. 绘制 EV 调节域
fprintf('正在绘制 EV 概率边界图像...\n');
fig_ev = figure('Name', 'EV Reliable Probabilistic Domain', 'Position', [150, 150, 1000, 600], 'Color', 'w');
hold on;

% 1. 绘制背景场景
plot(new_time_points, EV_Up_Shifted, 'Color', [0.6, 0.8, 1, 0.08], 'HandleVisibility', 'off');
plot(new_time_points, EV_Down_Shifted, 'Color', [1, 0.6, 0.6, 0.08], 'HandleVisibility', 'off');

% 2. 绘制边界
% 上调
p80_ev = plot(new_time_points, Bounds_EV_Up(:,3), 'b:', 'LineWidth', 1.5, 'DisplayName', '80% 置信边界');
p95_ev = plot(new_time_points, Bounds_EV_Up(:,2), 'b--', 'LineWidth', 2.0, 'DisplayName', '95% 置信边界');
p100_ev = plot(new_time_points, Bounds_EV_Up(:,1), 'b-', 'LineWidth', 2.5, 'DisplayName', '100% 可靠边界');

% 下调
plot(new_time_points, Bounds_EV_Down(:,3), 'r:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(new_time_points, Bounds_EV_Down(:,2), 'r--', 'LineWidth', 2.0, 'HandleVisibility', 'off');
plot(new_time_points, Bounds_EV_Down(:,1), 'r-', 'LineWidth', 2.5, 'HandleVisibility', 'off');

% 0轴线
yline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

% 装饰
grid on;
xlim([0, 24]);
% 图例设置 (字体放大)
legend([p100_ev, p95_ev, p80_ev], 'Location', 'best', 'FontSize', 18, 'FontName', default_font);
set(gca, 'XTick', x_ticks_new, 'XTickLabel', x_labels_new, ...
    'FontSize', axis_font_size, 'FontName', default_font, 'LineWidth', 1.2);
xlabel('时间', 'FontSize', label_font_size, 'FontName', default_font); 
ylabel('EV 功率 (kW)', 'FontSize', label_font_size, 'FontName', default_font);

% 保存
save_name_ev = 'EV_可靠调节域_概率边界_100scenarios.png';
print(fig_ev, save_name_ev, '-dpng', '-r600');
fprintf('EV 图像已保存: %s\n', save_name_ev);

fprintf('所有程序执行完毕。\n');