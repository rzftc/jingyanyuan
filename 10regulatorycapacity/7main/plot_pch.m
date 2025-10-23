%% 绘制特定激励下，EV集群总充电功率与入网车辆数的时变曲线
% 使用双Y轴展示两种数据
% 横坐标表示实际时间：从第一天早上6点到第二天早上6点
% 数据从 main_diff_delt_48.m 等脚本运行后保存的 .mat 文件中加载
% -------------------------------------------------------------------------

clear; close all; clc; % 清理工作区和图形

% --- 1. 参数和文件设置 ---
% !!! 请确保下面的文件名与您实际保存的对应特定dt的文件名完全一致 !!!
mat_filename = 'dt_5m.mat'; % 例如: 'results_3D_dt_5min.mat'

% 【修改】选择一个要绘制曲线的激励电价的索引
% 默认的激励范围是 linspace(0, 50, 25)，我们选择中间的第13个作为代表
p_idx_for_plot = 13; 
fprintf('将为激励电价索引 %d 绘制集群总充电功率与入网车辆数曲线。\n', p_idx_for_plot);

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

% 检查所需字段是否存在
if ~isfield(results_3D_current_file, 'EV_Total_Charge_Power')
    error('错误: 文件 %s 的 results_3D 结构体中缺少 EV_Total_Charge_Power 字段。', mat_filename);
end
if ~isfield(results_3D_current_file, 'Online_EV_Count')
    error('错误: 文件 %s 的 results_3D 结构体中缺少 Online_EV_Count 字段。', mat_filename);
end

% 获取仿真参数
% 默认仿真起始于第一天的早上6点
simulation_start_offset_hours = 6.0; 
if isfield(results_3D_current_file, 'simulation_start_hour_actual') && ~isempty(results_3D_current_file.simulation_start_hour_actual)
    simulation_start_offset_hours = results_3D_current_file.simulation_start_hour_actual;
end

% 获取仿真的相对时间点 (从0开始)
time_points_relative = []; 
if isfield(results_3D_current_file, 'time_points_actual') && ~isempty(results_3D_current_file.time_points_actual)
     time_points_relative = results_3D_current_file.time_points_actual(:) - results_3D_current_file.time_points_actual(1);
elseif isfield(loaded_data, 'time_points') && ~isempty(loaded_data.time_points) 
     time_points_relative = loaded_data.time_points(:) - loaded_data.time_points(1);
else 
    error('错误: 文件 %s 中无法确定时间轴信息。', mat_filename);
end

p_incentive_range = linspace(0, 50, 25); % 默认值
if isfield(results_3D_current_file, 'p_incentive_range_actual') && ~isempty(results_3D_current_file.p_incentive_range_actual)
    p_incentive_range = results_3D_current_file.p_incentive_range_actual;
end
if p_idx_for_plot > length(p_incentive_range) || p_idx_for_plot < 1
    error('错误: p_idx_for_plot (%d) 超出了 p_incentive_range (长度 %d) 的范围。', p_idx_for_plot, length(p_incentive_range));
end

% 提取所需的数据序列
cluster_total_charge_power = results_3D_current_file.EV_Total_Charge_Power(:, p_idx_for_plot);
online_ev_count = results_3D_current_file.Online_EV_Count(:, p_idx_for_plot);

% --- 3. 绘制双Y轴图 ---
fig_combined = figure('Position', [200 200 1000 600]);
ax = gca; % 获取当前坐标轴句柄

% 生成绝对时间轴用于绘图
time_axis_to_plot = time_points_relative + simulation_start_offset_hours;

% 左Y轴 - 充电功率
yyaxis left
plot(time_axis_to_plot, cluster_total_charge_power, 'b-', 'LineWidth', 2.0);
ylabel('集群总充电功率 (kW)', 'FontSize', 18, 'Color', 'b');
ax.YColor = 'b';
grid on;

% 右Y轴 - 入网车辆数
yyaxis right
plot(time_axis_to_plot, online_ev_count, 'r--', 'LineWidth', 2.0);
ylabel('入网电动汽车数量', 'FontSize', 18, 'Color', 'r');
ax.YColor = 'r';

% 图形整体设置
hold off;
legend({'总充电功率', '入网EV数量'}, 'Location', 'best', 'FontSize', 16);
set(ax, 'FontSize', 18);

% --- 4. 设置横坐标（X轴）格式 ---
sim_start_absolute_hour = time_axis_to_plot(1);
sim_end_absolute_hour = time_axis_to_plot(end);

% 设置X轴范围为数据的完整时间范围
xlim([sim_start_absolute_hour, sim_end_absolute_hour]);

% 生成并设置X轴刻度和标签
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

% 设置X轴标题
xlabel_start_hour = mod(sim_start_absolute_hour, 24);
xlabel_end_hour = mod(sim_end_absolute_hour, 24);
if abs(xlabel_end_hour - xlabel_start_hour) < 1e-6
    xlabel_text = sprintf('时间 (Day 1 %.0f:00 至 Day 2 %.0f:00)', xlabel_start_hour, xlabel_start_hour);
else
    xlabel_text = sprintf('时间 (从 %.0f:00 开始)', xlabel_start_hour);
end
xlabel(xlabel_text, 'FontSize', 18);

% --- 5. 文件名和保存 ---
file_dt_minutes_str = '';
if isfield(results_3D_current_file, 'dt_actual') && ~isempty(results_3D_current_file.dt_actual)
    file_dt_minutes_str = sprintf('_dt%dmin', round(results_3D_current_file.dt_actual*60));
end
plot_filename_combined = sprintf('Cluster_Power_vs_OnlineCount_P%d%s_48hr.png', p_idx_for_plot, file_dt_minutes_str);
print(fig_combined, plot_filename_combined, '-dpng', '-r400');
fprintf('集群功率与在线车辆数组合图已保存为 %s\n', plot_filename_combined);

fprintf('绘图完成。\n');