%% plot_single_ev_details_v2.m
%
% 功能:
% 1. 加载由 (修正后的) ac_ev_simulation_block_abselute_hour.m 
%    生成的 .mat 结果文件 (必须包含 E_current_EV, P_current_EV, SOC_EV)。
% 2. 加载原始的 EV Excel 输入文件 (获取 E_in, C_EV 等元数据)。
% 3. 允许用户在下方 "用户配置参数" 部分指定要绘制的 EV 索引。
% 4. 绘制该 EV 的三个独立曲线：
%    - 图 1: 充电功率 (P_current_EV)
%    - 图 2: 电池能量 (E_current_EV)
%    - 图 3: 虚拟 SOC (SOC_EV)
% 5. 格式化X轴以显示 "6:00" 到 "第二天 6:00" 的绝对时间。

clear; close all; clc;

%% 1. --- 用户配置参数 ---

% !!! (用户必须修改) 指向您 *新生成* 的 .mat 结果文件
mat_file_to_load = 'chunk_results_abs_hour/results_chunk_1.mat'; 

% !!! (用户必须修改) 指向用于生成该结果的 EV Excel 文件
ev_excel_file = '2EV_residential.xlsx'; 

% !!! (用户必须修改) 选择要绘制的 EV 索引
ev_index_to_plot = 5; 

% --- (可选) 自动检测字段名 (如果您的保存脚本不同) ---
field_name_e_current = 'E_current_EV';
field_name_p_current = 'P_current_EV'; % [!!! 新增 !!!]
field_name_soc_virtual = 'SOC_EV';

%% 2. --- 加载数据 ---

% --- 加载 .mat 结果文件 ---
fprintf('正在加载 .mat 结果文件: %s\n', mat_file_to_load);
if ~exist(mat_file_to_load, 'file')
    error('错误: 未找到 .mat 文件 "%s"。', mat_file_to_load);
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

% --- 加载 Excel 输入文件 ---
fprintf('正在加载 Excel 输入文件: %s\n', ev_excel_file);
if ~exist(ev_excel_file, 'file')
    error('错误: 未找到 Excel 文件 "%s"。', ev_excel_file);
end
try
    EVs_all = initializeEVsFromExcel(ev_excel_file);
catch ME
    error('加载或初始化 Excel 文件时出错 (请确保 initializeEVsFromExcel.m 在路径中): %s', ME.message);
end

%% 3. --- 提取和验证数据 ---

% --- 提取时间轴 ---
% (ac_ev_simulation_block_abselute_hour.m V2 已保存这些)
if isfield(results, 'time_points_absolute')
    time_axis_absolute = results.time_points_absolute(:);
else
    error('错误: 在 "results" 结构体中未找到 "time_points_absolute"。请确保您使用的是 V2 版本的仿真脚本。');
end
if isfield(results, 'dt')
    dt = results.dt;
else
     error('错误: 在 "results" 结构体中未找到 "dt"。');
end
fprintf('已加载时间轴 (%.0f:00 到 %.0f:00, dt=%.4f)\n', ...
    time_axis_absolute(1), time_axis_absolute(end), dt);

% --- 提取 EV 特定参数 (从Excel加载的完整列表中) ---
num_ev_in_excel = length(EVs_all);

chunkIndex_str = regexp(mat_file_to_load, 'chunk_(\d+)', 'tokens');
chunkSize = 10000; % 假设分块大小
absolute_ev_index = ev_index_to_plot; 

if ~isempty(chunkIndex_str)
    chunkNum = str2double(chunkIndex_str{1}{1});
    absolute_ev_index = (chunkNum - 1) * chunkSize + ev_index_to_plot;
    fprintf('检测到文件为 分块 %d。绘图索引 %d 对应 Excel 绝对索引 %d。\n', ...
        chunkNum, ev_index_to_plot, absolute_ev_index);
else
    fprintf('未检测到分块编号。假设 ev_index_to_plot (%d) 是 Excel 中的绝对索引。\n', ev_index_to_plot);
end

if absolute_ev_index < 1 || absolute_ev_index > num_ev_in_excel
    error('错误: 计算出的绝对索引 (%d) 超出了 Excel 文件中的 EV 数量 (1 到 %d)。', ...
        absolute_ev_index, num_ev_in_excel);
end

% 从Excel获取该EV的元数据
E_in_ev = EVs_all(absolute_ev_index).E_in;
C_EV_ev = EVs_all(absolute_ev_index).C_EV;

% --- 自动检测并提取 E_current, P_current, 和 SOC_EV 字段 ---
if ~isfield(results, field_name_e_current)
    error('错误: 在 "results" 结构体中未找到 "%s" 字段。', field_name_e_current);
end
if ~isfield(results, field_name_p_current)
    error('错误: 在 "results" 结构体中未找到 "%s" 字段。', field_name_p_current);
end
if ~isfield(results, field_name_soc_virtual)
    error('错误: 在 "results" 结构体中未找到 "%s" 字段。', field_name_soc_virtual);
end
fprintf('使用字段: "%s", "%s", "%s"\n', ...
    field_name_e_current, field_name_p_current, field_name_soc_virtual);

% --- 验证索引是否在 .mat 结果范围内 ---
num_ev_in_mat = size(results.(field_name_e_current), 1);
if ev_index_to_plot > num_ev_in_mat
    error('错误: ev_index_to_plot (%d) 超出了 .mat 文件中的 EV 数量 (%d)。', ...
        ev_index_to_plot, num_ev_in_mat);
end

% --- 提取历史数据 (使用相对索引 ev_index_to_plot) ---
E_current_history = results.(field_name_e_current)(ev_index_to_plot, :);
P_current_history = results.(field_name_p_current)(ev_index_to_plot, :); % [!!! 修改 !!!]
SOC_virtual_history = results.(field_name_soc_virtual)(ev_index_to_plot, :);


%% 4. --- 绘图 ---
% (注意: 第4节 "派生充电功率" 已被移除，因为我们现在直接加载 P_current_history)
fprintf('正在生成图表...\n');

% --- 图 1: 充电功率 ---
fig1 = figure('Name', sprintf('EV %d (Abs Index %d) - 充电功率', ev_index_to_plot, absolute_ev_index), 'Position', [100, 600, 1000, 450]);
ax1 = axes(fig1);
% [!!! 修改 !!!] 直接绘制 P_current_history
plot(ax1, time_axis_absolute, P_current_history, 'r-', 'LineWidth', 1.5);
hold(ax1, 'on'); grid on;
ylabel(ax1, '充电功率 (kW)', 'FontSize', 12);
title(ax1, sprintf('EV %d (Excel行 %d) - 充电功率 (P_{current})', ev_index_to_plot, absolute_ev_index+1), 'FontSize', 14);
yline(ax1, 0, 'k:');

% --- 图 2: 电池能量 ---
fig2 = figure('Name', sprintf('EV %d (Abs Index %d) - 电池能量', ev_index_to_plot, absolute_ev_index), 'Position', [100, 350, 1000, 450]);
ax2 = axes(fig2);
plot(ax2, time_axis_absolute, E_current_history, 'b-', 'LineWidth', 1.5, 'DisplayName', 'E_{current} (kWh)');
hold(ax2, 'on'); grid on;
yline(ax2, E_in_ev, 'k:', 'LineWidth', 1, 'DisplayName', sprintf('E_{in} (%.1f kWh)', E_in_ev));
yline(ax2, C_EV_ev, 'k--', 'LineWidth', 1, 'DisplayName', sprintf('C_{EV} (%.1f kWh)', C_EV_ev));
ylabel(ax2, '电池能量 (kWh)', 'FontSize', 12);
title(ax2, sprintf('EV %d (Excel行 %d) - 电池能量 (E_{current})', ev_index_to_plot, absolute_ev_index+1), 'FontSize', 14);
legend(ax2, 'Location', 'best');

% --- 图 3: 虚拟 SOC ---
fig3 = figure('Name', sprintf('EV %d (Abs Index %d) - 虚拟 SOC', ev_index_to_plot, absolute_ev_index), 'Position', [100, 100, 1000, 450]);
ax3 = axes(fig3);
plot(ax3, time_axis_absolute, SOC_virtual_history, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Virtual SOC');
hold(ax3, 'on'); grid on;
try
    SOC_physical_history = E_current_history / C_EV_ev;
    plot(ax3, time_axis_absolute, SOC_physical_history, 'k:', 'LineWidth', 1.0, 'DisplayName', 'Physical SOC (E_{current} / C_{EV})');
catch
end
ylabel(ax3, '虚拟 SOC', 'FontSize', 12);
title(ax3, sprintf('EV %d (Excel行 %d) - 虚拟 SOC (SOC_{EV})', ev_index_to_plot, absolute_ev_index+1), 'FontSize', 14);
legend(ax3, 'Location', 'best');
ylim(ax3, 'auto'); 

%% 5. --- 统一格式化 X 轴 ---
all_axes = [ax1; ax2; ax3];
sim_start_absolute_hour_plot = time_axis_absolute(1);
sim_end_absolute_hour_plot = time_axis_absolute(end);

for ax = all_axes'
    xlim(ax, [sim_start_absolute_hour_plot, sim_end_absolute_hour_plot]);
    
    xticks_to_set = sim_start_absolute_hour_plot:6:sim_end_absolute_hour_plot;
    
    xtick_labels_to_set = {};
    for tick_hour = xticks_to_set
        hour_of_day_display = mod(tick_hour, 24);
        
        if tick_hour >= 24 && hour_of_day_display < sim_start_absolute_hour_plot
            xtick_labels_to_set{end+1} = sprintf('%.0f:00 (次日)', hour_of_day_display);
        else
            xtick_labels_to_set{end+1} = sprintf('%.0f:00', hour_of_day_display);
        end
    end
    
    set(ax, 'XTick', xticks_to_set, 'XTickLabel', xtick_labels_to_set, 'FontSize', 11);
    
    xlabel_start_hour_plot = mod(sim_start_absolute_hour_plot, 24);
    xlabel_text = sprintf('时间 (Day 1 %.0f:00 至 Day 2 %.0f:00)', xlabel_start_hour_plot, xlabel_start_hour_plot);
    xlabel(ax, xlabel_text, 'FontSize', 12);
end

fprintf('绘图完成。\n');