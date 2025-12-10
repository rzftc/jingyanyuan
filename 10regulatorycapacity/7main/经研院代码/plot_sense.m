%% 绘制特定dt下，不同激励电价对应的最大集群调节能力
% 上调和下调潜力分别在不同的图上显示
% 横坐标为激励电价
% 数据从单个 main_diff_delt_48.m 运行结果 .mat 文件中加载
% -------------------------------------------------------------------------

clear; close all; clc; % 清理工作区和图形

% --- 1. 参数和文件设置 ---
% !!! 请确保下面的文件名与您实际保存的对应特定dt的文件名完全一致 !!!
mat_filename = 'dt_5_60.mat'; % 例如: 'results_3D_dt_3min_sim_starts_at_6.mat'
                             % 这是您要分析的、包含特定dt下所有激励电价结果的文件

p_incentive_range_default = linspace(0, 50, 25); % 默认的激励电价范围

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

if ~isfield(results_3D_current_file, 'EV_Up') || ~isfield(results_3D_current_file, 'EV_Down')
    error('错误: 文件 %s 中缺少 EV_Up 或 EV_Down 字段。', mat_filename);
end

% 获取激励电价范围
p_incentive_range = p_incentive_range_default; % 默认值
if isfield(results_3D_current_file, 'p_incentive_range_actual') && ~isempty(results_3D_current_file.p_incentive_range_actual)
    p_incentive_range = results_3D_current_file.p_incentive_range_actual;
elseif isfield(loaded_data, 'p_incentive_range') && ~isempty(loaded_data.p_incentive_range) % 兼容旧的保存方式
    p_incentive_range = loaded_data.p_incentive_range;
end
if size(results_3D_current_file.EV_Up, 2) ~= length(p_incentive_range)
    warning('警告: results_3D.EV_Up 的列数与 p_incentive_range 长度不匹配，请检查数据一致性。将尝试使用 p_incentive_range 的长度。');
    if size(results_3D_current_file.EV_Up, 2) < length(p_incentive_range)
        p_incentive_range = p_incentive_range(1:size(results_3D_current_file.EV_Up, 2));
    end
end


% --- 3. 计算最大调节能力 ---
% results_3D.EV_Up 和 EV_Down 的维度是 [length(time_points), length(p_incentive_range)]
cluster_up_all_prices_time_series = results_3D_current_file.EV_Up;
cluster_down_all_prices_time_series = results_3D_current_file.EV_Down;

% 计算每个激励电价下的最大上调潜力
max_up_capacity = zeros(1, length(p_incentive_range));
for p_idx = 1:length(p_incentive_range)
    if p_idx <= size(cluster_up_all_prices_time_series, 2)
        max_up_capacity(p_idx) = max(cluster_up_all_prices_time_series(:, p_idx));
    else
        max_up_capacity(p_idx) = NaN; % 如果数据列不足
    end
end

% 计算每个激励电价下的最大下调潜力（绝对值）
% EV_Down 通常为负值或0，所以取其绝对值的最大值，或者 -min(EV_Down)
max_down_capacity_magnitude = zeros(1, length(p_incentive_range));
for p_idx = 1:length(p_incentive_range)
    if p_idx <= size(cluster_down_all_prices_time_series, 2)
        max_down_capacity_magnitude(p_idx) = max(abs(cluster_down_all_prices_time_series(:, p_idx)));
        % 或者: max_down_capacity_magnitude(p_idx) = -min(cluster_down_all_prices_time_series(:, p_idx)); 
        % 两者在EV_Down <= 0 时等价
    else
        max_down_capacity_magnitude(p_idx) = NaN;
    end
end

% 获取当前 dt 用于文件名 (如果存在)
file_dt_minutes_str = '';
if isfield(results_3D_current_file, 'dt_actual') && ~isempty(results_3D_current_file.dt_actual)
    file_dt_minutes_str = sprintf('_dt%dmin', round(results_3D_current_file.dt_actual*60));
end

% --- 4. 绘制最大上调潜力曲线 ---
fig_up = figure('Position', [200 200 1000 600]);
plot(p_incentive_range, max_up_capacity, ...
     'LineWidth', 1.8, ...
     'Marker', 'o', ...
     'Color', [0 0.4470 0.7410]); % 蓝色
grid on;
xlabel('激励电价 (分/kW)', 'FontSize', 18);
ylabel('最大集群上调潜力 (kW)', 'FontSize', 18);
ax_up = gca; ax_up.FontSize = 18;
% Y轴从0开始，并根据数据动态调整上限
ylim_min_up = 0;
ylim_max_up = max(max_up_capacity(:)) * 1.1;
if ylim_max_up < 1, ylim_max_up = 1; end % 避免上限过小
if isnan(ylim_max_up), ylim_max_up = 1; end
ylim([ylim_min_up, ylim_max_up]);

plot_filename_up = sprintf('Max_Cluster_Up_Capacity_vs_Price%s.png', file_dt_minutes_str);
print(fig_up, plot_filename_up, '-dpng', '-r400');
fprintf('最大上调潜力图已保存为 %s\n', plot_filename_up);

% --- 5. 绘制最大下调潜力曲线 ---
fig_down = figure('Position', [200 850 1000 600]); 
plot(p_incentive_range, max_down_capacity_magnitude, ... 
     'LineWidth', 1.8, ...
     'Marker', 's', ...
     'Color', [0.8500 0.3250 0.0980]); % 橙色
grid on;
xlabel('激励电价 (分/kW)', 'FontSize', 18);
ylabel('最大集群下调潜力 (kW)', 'FontSize', 18); 
ax_down = gca; ax_down.FontSize = 18;
% Y轴从0开始，并根据数据动态调整上限
ylim_min_down = 0;
ylim_max_down = max(max_down_capacity_magnitude(:)) * 1.1;
if ylim_max_down < 1, ylim_max_down = 1; end % 避免上限过小
if isnan(ylim_max_down), ylim_max_down = 1; end
ylim([ylim_min_down, ylim_max_down]);

plot_filename_down = sprintf('Max_Cluster_Down_Capacity_vs_Price%s.png', file_dt_minutes_str);
print(fig_down, plot_filename_down, '-dpng', '-r400');
fprintf('最大下调潜力图已保存为 %s\n', plot_filename_down);

fprintf('所有绘图完成。\n');