%% 绘制不同dt(调节时长)下、特定激励电价的集群总调节能力
% 上调和下调潜力分别在不同的图上显示
% 横坐标表示实际时间：从第一天早上6点到第二天早上6点
% 数据从 main_diff_delt_48.m 运行后保存的 .mat 文件中加载
% -------------------------------------------------------------------------

clear; close all; clc; % 清理工作区和图形

% --- 1. 参数和文件设置 ---
% !!! 请根据您实际的文件名和对应的dt值进行修改 !!!
dt_values_for_plot_hours = [5/60, 15/60, 1]; % 您运行 main_diff_delt_48.m 时使用的 dt 值 (小时)
results_files = {
    'dt_5m.mat',   
    'dt_15m.mat',    
    'dt_60m.mat'
};
p_idx_for_plot = 13; % 选择要绘制的激励电价的索引

% --- 2. 加载和处理数据 ---
all_time_points_absolute = cell(1, length(results_files));
all_cluster_up_potentials = cell(1, length(results_files));
all_cluster_down_potentials = cell(1, length(results_files));
actual_dt_minutes_loaded = zeros(1, length(results_files));
valid_files_count = 0;
unified_simulation_start_offset = NaN; 

for k = 1:length(results_files)
    mat_filename = results_files{k};
    if exist(mat_filename, 'file')
        fprintf('加载文件: %s\n', mat_filename);
        loaded_data = load(mat_filename); 

        if ~isfield(loaded_data, 'results_3D') || ~isstruct(loaded_data.results_3D)
            warning('警告: 文件 %s 中未找到 results_3D 结构体，跳过。', mat_filename);
            continue;
        end
        results_3D_current_file = loaded_data.results_3D;
        
        % 提取仿真开始时间
        sim_start_hour_this_file = 6.0; % 默认值
        if isfield(results_3D_current_file, 'simulation_start_hour_actual')
            sim_start_hour_this_file = results_3D_current_file.simulation_start_hour_actual;
        end
        if isnan(unified_simulation_start_offset)
            unified_simulation_start_offset = sim_start_hour_this_file;
        end

        % 提取绝对时间轴
        time_points_abs_from_file = [];
        if isfield(results_3D_current_file, 'time_points_actual')
            time_points_abs_from_file = results_3D_current_file.time_points_actual(:);
        else
            warning('警告: 文件 %s 缺少时间轴信息，跳过。', mat_filename);
            continue;
        end

        valid_files_count = valid_files_count + 1;
        all_time_points_absolute{valid_files_count} = time_points_abs_from_file;
        all_cluster_up_potentials{valid_files_count} = results_3D_current_file.EV_Up(:, p_idx_for_plot);
        all_cluster_down_potentials{valid_files_count} = results_3D_current_file.EV_Down(:, p_idx_for_plot);
        actual_dt_minutes_loaded(valid_files_count) = round(dt_values_for_plot_hours(k) * 60);
    else
        warning('警告: 文件 %s 未找到，跳过。', mat_filename);
    end
end

if valid_files_count == 0
    error('错误: 未能成功加载任何有效数据文件，无法绘图。');
end

% --- 3. 设置绘图参数 ---
colors = lines(valid_files_count);
sim_start_absolute_hour = unified_simulation_start_offset;
sim_end_absolute_hour = sim_start_absolute_hour + 24;

% --- 4. 绘制上调潜力曲线 ---
fig_up = figure('Position', [200 200 1000 600]);
hold on; grid on;
legend_entries_up = {};
for k_plot = 1:valid_files_count
    plot(all_time_points_absolute{k_plot}, all_cluster_up_potentials{k_plot}, ...
         'LineWidth', 1.8, 'Color', colors(k_plot,:));
    legend_entries_up{end+1} = sprintf('%d 分钟', actual_dt_minutes_loaded(k_plot));
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
for k_plot = 1:valid_files_count
    plot(all_time_points_absolute{k_plot}, all_cluster_down_potentials{k_plot}, ... 
         'LineWidth', 1.8, 'Color', colors(k_plot,:));
    legend_entries_down{end+1} = sprintf('%d 分钟', actual_dt_minutes_loaded(k_plot));
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
ylim_dynamic_plot(findobj(fig_up, 'type', 'axes'), all_cluster_up_potentials, true); 
ylim_dynamic_plot(findobj(fig_down, 'type', 'axes'), all_cluster_down_potentials, false); 

% --- 7. 保存图像 ---
print(fig_up, 'Cluster_Up_Capacity_VarDT.png', '-dpng', '-r400');
fprintf('上调潜力图已保存为 Cluster_Up_Capacity_VarDT.png\n');
print(fig_down, 'Cluster_Down_Capacity_VarDT.png', '-dpng', '-r400');
fprintf('下调潜力图已保存为 Cluster_Down_Capacity_VarDT.png\n');

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