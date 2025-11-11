%% 绘制由 ac_ev_simulation_block_abselute_hour.m 生成的 EV 潜力数据
%
% 功能:
% 1. 加载指定的 .mat 文件 (例如 'chunk_results_abs_hour/results_chunk_1.mat')。
% 2. 提取 EV_Up 和 EV_Down 聚合数据。
% 3. 创建一个从 "早上6点到第二天6点" (6:00 到 30:00) 的时间轴。
% 4. 绘制 EV_Up 和 EV_Down 在两个不同子图上的曲线。
% 5. 格式化 X 轴以显示 "6:00", "12:00", "0:00 (次日)" 等标签。

clear; close all; clc;

%% 1. --- 参数配置 ---

% !!! 请修改这里，指向您生成的 .mat 文件 !!!
mat_file_to_load = 'chunk_results_abs_hour/results_chunk_1.mat';

% --- 时间轴参数 (必须与 ac_ev_simulation_block_abselute_hour.m 脚本中的定义一致) ---
simulation_start_hour = 6;  % 仿真开始于 6:00
simulation_end_hour   = 30; % 仿真结束于次日 6:00 (24 + 6)
dt = 0.05;                  % 时间分辨率（小时）

%% 2. --- 加载数据 ---
fprintf('正在加载数据文件: %s\n', mat_file_to_load);
if ~exist(mat_file_to_load, 'file')
    error('错误: 未找到数据文件 "%s"。\n请先运行 ac_ev_simulation_block_abselute_hour.m 生成该文件。', mat_file_to_load);
end

try
    loaded_data = load(mat_file_to_load);
    if ~isfield(loaded_data, 'results')
        error('错误: 在 "%s" 中未找到 "results" 结构体。', mat_file_to_load);
    end
    results = loaded_data.results;
catch ME
    error('加载 .mat 文件时出错: %s', ME.message);
end

% --- 提取 EV 数据 ---
if ~isfield(results, 'EV_Up') || ~isfield(results, 'EV_Down')
    warning('警告: "results" 结构体中缺少 "EV_Up" 或 "EV_Down" 字段。');
    return;
end
ev_up_data = results.EV_Up;
ev_down_data = results.EV_Down;

fprintf('数据加载成功。\n');

%% 3. --- 创建和验证时间轴 ---
% (此逻辑基于 ac_ev_simulation_block_abselute_hour.m)
time_axis_absolute = (simulation_start_hour:dt:simulation_end_hour)'; % 创建绝对时间轴

% 确保数据长度与时间轴匹配
if length(time_axis_absolute) ~= length(ev_up_data)
    warning('时间轴长度 (%d) 与数据长度 (%d) 不匹配！将截断为最短长度。', ...
            length(time_axis_absolute), length(ev_up_data));
    min_len = min(length(time_axis_absolute), length(ev_up_data));
    time_axis_absolute = time_axis_absolute(1:min_len);
    ev_up_data = ev_up_data(1:min_len);
    ev_down_data = ev_down_data(1:min_len);
end

%% 4. --- 绘制图形 ---
fprintf('正在生成图表...\n');

% --- 图 1: EV 上调潜力 ---
fig1 = figure('Name', 'EV 上调潜力 (EV_Up)', 'Position', [100, 400, 1000, 450]);
ax1 = axes(fig1);
plot(ax1, time_axis_absolute, ev_up_data, 'b-', 'LineWidth', 1.5, 'DisplayName', 'EV_Up');
hold(ax1, 'on'); grid on;
ylabel(ax1, '上调潜力 (kW)', 'FontSize', 12);
legend(ax1, 'Location', 'best');
title(ax1, 'EV 集群上调潜力 (EV_Up)', 'FontSize', 14);

% --- 图 2: EV 下调潜力 ---
fig2 = figure('Name', 'EV 下调潜力 (EV_Down)', 'Position', [100, 100, 1000, 450]);
ax2 = axes(fig2);
plot(ax2, time_axis_absolute, ev_down_data, 'r-', 'LineWidth', 1.5, 'DisplayName', 'EV_Down');
hold(ax2, 'on'); grid on;
ylabel(ax2, '下调潜力 (kW)', 'FontSize', 12);
legend(ax2, 'Location', 'best');
title(ax2, 'EV 集群下调潜力 (EV_Down)', 'FontSize', 14);

%% 5. --- 统一格式化 X 轴 ---
% (此 X 轴格式化逻辑参考了 plot_pch.m 和 plot_diffr.m)
all_axes = [ax1; ax2]; % 获取两个图的坐标轴句柄

for ax = all_axes'
    % --- 设置 X 轴范围 ---
    xlim(ax, [simulation_start_hour, simulation_end_hour]);
    
    % --- 生成刻度位置 ---
    xticks_to_set = simulation_start_hour:6:simulation_end_hour; % 每6小时一个刻度
    
    % --- 生成刻度标签 ---
    xtick_labels_to_set = {};
    for tick_hour = xticks_to_set
        hour_of_day_display = mod(tick_hour, 24); % 计算24小时制的小时
        
        % 检查是否是第二天
        if tick_hour >= 24 && hour_of_day_display < simulation_start_hour
            xtick_labels_to_set{end+1} = sprintf('%.0f:00 (次日)', hour_of_day_display);
        else
            xtick_labels_to_set{end+1} = sprintf('%.0f:00', hour_of_day_display);
        end
    end
    
    % --- 应用刻度和标签 ---
    set(ax, 'XTick', xticks_to_set, 'XTickLabel', xtick_labels_to_set, 'FontSize', 11);
    
    % --- 设置 X 轴总标签 ---
    xlabel_text = sprintf('时间 (Day 1 %.0f:00 至 Day 2 %.0f:00)', mod(simulation_start_hour, 24), mod(simulation_end_hour, 24));
    xlabel(ax, xlabel_text, 'FontSize', 12);
end

fprintf('绘图完成。\n');