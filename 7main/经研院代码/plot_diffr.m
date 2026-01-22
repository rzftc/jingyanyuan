%% 绘制特定dt(调节时长)下、不同激励电价的集群总调节能力
% 上调和下调潜力分别在不同的图上显示，下调保留原始符号
% 横坐标表示实际时间：从第一天早上6点到第二天早上6点
% 数据从 main_diff_delt_48.m 等脚本运行后保存的 .mat 文件中加载
% -------------------------------------------------------------------------

clear; close all; clc; % 清理工作区和图形

% --- 1. 参数和文件设置 ---
mat_filename = 'dt_5m.mat'; % 确保文件名与实际保存的文件一致
p_indices_for_plot = [3, 8, 13, 18, 23]; % 选择要绘制曲线的激励电价的索引

% --- 2. 加载数据 ---
if exist(mat_filename, 'file')
    fprintf('加载文件: %s\n', mat_filename);
    loaded_data = load(mat_filename); 
else
    error('错误: 数据文件 %s 未找到。', mat_filename);
end

if ~isfield(loaded_data, 'results_3D') || ~isstruct(loaded_data.results_3D)
    error('错误: 文件 %s 中未找到 results_3D 结构体。', mat_filename);
end
results_3D_current_file = loaded_data.results_3D;

% 获取仿真参数
simulation_start_offset_hours = 6.0; % 默认值
if isfield(results_3D_current_file, 'simulation_start_hour_actual')
    simulation_start_offset_hours = results_3D_current_file.simulation_start_hour_actual;
end

time_points_absolute = [];
if isfield(results_3D_current_file, 'time_points_actual')
    time_points_absolute = results_3D_current_file.time_points_actual(:);
else
    error('错误: 文件 %s 缺少时间轴信息。', mat_filename);
end

p_incentive_range = linspace(0, 50, 25); % 默认值
if isfield(results_3D_current_file, 'p_incentive_range_actual')
    p_incentive_range = results_3D_current_file.p_incentive_range_actual;
end
if max(p_indices_for_plot) > length(p_incentive_range) || min(p_indices_for_plot) < 1
    error('错误: p_indices_for_plot 中的索引超出了 p_incentive_range 的范围。');
end

% --- 3. 设置绘图参数 ---
sim_start_absolute_hour = time_points_absolute(1);
sim_end_absolute_hour = time_points_absolute(end);
colors_up = lines(length(p_indices_for_plot));
colors_down = jet(length(p_indices_for_plot));

% --- 4. 绘制上调潜力曲线 ---
fig_up = figure('Position', [200 200 1000 600]);
hold on; grid on;
legend_entries_up = {};
for k = 1:length(p_indices_for_plot)
    p_idx = p_indices_for_plot(k);
    current_price = p_incentive_range(p_idx);
    cluster_up_potential = results_3D_current_file.EV_Up(:, p_idx);
    
    plot(time_points_absolute, cluster_up_potential, 'LineWidth', 1.8, 'Color', colors_up(k,:));
    legend_entries_up{end+1} = sprintf('激励 %.1f 分/kW', current_price);
end
hold off;

ax_up = gca;
ax_up.FontSize = 18;
ylabel('集群上调潜力 (kW)', 'FontSize', 18);
if ~isempty(legend_entries_up)
    legend(legend_entries_up, 'Location', 'best', 'FontSize', 14);
end

% --- 5. 绘制下调潜力曲线 ---
fig_down = figure('Position', [200 850 1000 600]); 
hold on; grid on;
legend_entries_down = {};
for k = 1:length(p_indices_for_plot)
    p_idx = p_indices_for_plot(k);
    current_price = p_incentive_range(p_idx);
    cluster_down_potential = results_3D_current_file.EV_Down(:, p_idx);
    
    plot(time_points_absolute, cluster_down_potential, 'LineWidth', 1.8, 'Color', colors_down(k,:));
    legend_entries_down{end+1} = sprintf('激励 %.1f 分/kW', current_price);
end
hold off;

ax_down = gca;
ax_down.FontSize = 18;
ylabel('集群下调潜力 (kW)', 'FontSize', 18); 
if ~isempty(legend_entries_down)
    legend(legend_entries_down, 'Location', 'best', 'FontSize', 14);
end

% --- 6. 统一格式化坐标轴 ---
figs = [fig_up, fig_down];
for i = 1:length(figs)
    ax = findobj(figs(i), 'type', 'axes');
    xlim(ax, [sim_start_absolute_hour, sim_end_absolute_hour]);
    
    xticks_to_set = sim_start_absolute_hour:6:sim_end_absolute_hour;
    xtick_labels_to_set = {};
    for tick_hour = xticks_to_set
        hour_of_day_display = mod(tick_hour, 24);
        if tick_hour >= 24 && hour_of_day_display < sim_start_absolute_hour
            xtick_labels_to_set{end+1} = sprintf('%.0f:00 (次日)', hour_of_day_display);
        else
            xtick_labels_to_set{end+1} = sprintf('%.0f:00', hour_of_day_display);
        end
    end
    set(ax, 'XTick', xticks_to_set, 'XTickLabel', xtick_labels_to_set);
    
    xlabel_start_hour = mod(sim_start_absolute_hour, 24);
    xlabel_text = sprintf('时间 (Day 1 %.0f:00 至 Day 2 %.0f:00)', xlabel_start_hour, xlabel_start_hour);
    xlabel(ax, xlabel_text, 'FontSize', 18);
end

% 动态调整Y轴范围
ylim_dynamic_plot(findobj(fig_up, 'type', 'axes'), {results_3D_current_file.EV_Up(:, p_indices_for_plot)}, true);
ylim_dynamic_plot(findobj(fig_down, 'type', 'axes'), {results_3D_current_file.EV_Down(:, p_indices_for_plot)}, false);

% --- 7. 保存图像 ---
file_dt_minutes_str = '';
if isfield(results_3D_current_file, 'dt_actual')
    file_dt_minutes_str = sprintf('_dt%dmin', round(results_3D_current_file.dt_actual*60));
end
plot_filename_up = sprintf('Cluster_Up_Capacity_DiffPrice%s.png', file_dt_minutes_str);
print(fig_up, plot_filename_up, '-dpng', '-r400');
fprintf('上调潜力图已保存为 %s\n', plot_filename_up);

plot_filename_down = sprintf('Cluster_Down_Capacity_DiffPrice%s.png', file_dt_minutes_str);
print(fig_down, plot_filename_down, '-dpng', '-r400');
fprintf('下调潜力图已保存为 %s\n', plot_filename_down);

fprintf('所有绘图完成。\n');

% --- 辅助函数 ---
function ylim_dynamic_plot(ax, data_cell_array, force_positive_range_start_from_zero)
    all_data_vector = [];
    for i_cell = 1:length(data_cell_array)
        if ~isempty(data_cell_array{i_cell})
            all_data_vector = [all_data_vector; data_cell_array{i_cell}(:)];
        end
    end
    
    if isempty(all_data_vector) || all(isnan(all_data_vector)), return; end
    
    min_val = min(all_data_vector);
    max_val = max(all_data_vector);
    
    if abs(max_val - min_val) < 1e-9, padding = 0.1; else, padding = (max_val - min_val) * 0.1; end
    
    ylim_lower = min_val - padding;
    ylim_upper = max_val + padding;

    if force_positive_range_start_from_zero && min_val >= -1e-9
        ylim_lower = 0;
    end
    
    ylim(ax, [ylim_lower, ylim_upper]);
end