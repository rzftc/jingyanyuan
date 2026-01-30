% plot_vpp_analysis_Optimal_Style_V2.m
%
% 脚本功能：
% 1. 加载 2030, 2035, 2040 年的VPP优化数据。
% 2. 绘制2030年详细图 (任务一)。
% 3. 绘制2030-2040年趋势对比图 (任务二)。
% 4. 【优化风格 V2】：
%    - 字体: 'Arial' (此处代码实际指定为 'Yahei' 以支持中文)
%    - 字体大小: 增大 (14pt 标签, 12pt 刻度)
%    - 【修改】颜色: 优化后的数据使用亮色，参考线使用黑色。
%    - 【修改】线条:
%        - 优化调度 (最重要): 粗实线 (2.0pt)
%        - 电网需求 (次重要): 黑色虚线 (1.2pt)
%        - 系统潜力 (参考): 黑色点线 (1.0pt)
%    - 【修改】背景: 移除灰色面积图 (area)。
%    - 坐标轴: 移除顶部和右侧 (Box='off')，刻度朝外 (TickDir='out')
%    - 图例: 无边框，放置在左上角
%    - 网格: 仅 Y 轴, 浅色
% 5. 【保留】：横坐标维持 "Day 1 6:00" 到 "Day 2 6:00"
% 6. 【新增】：所有绘图数据除以10，保存为PNG格式。

clear; close all; clc;

%% 1. 全局配置与数据加载
fprintf('开始加载绘图数据...\n');

% --- 文件配置 ---
data_files = {
    'multi_obj_plot_data.mat',  % 2030年
    'multi_obj_plot_data2.mat', % 2035年
    'multi_obj_plot_data3.mat'  % 2040年
};
years = [2030, 2035, 2040];

% --- 【保留】时间轴转换参数 ---
% 基于 ac_ev_simulation_block.m (T=24, dt=0.05)
dt = 0.05; % 5 分钟 (0.0833) 小时步长
simulation_start_hour = 6.0; % 假设时间步 1 对应 6:00
fprintf('应用时间轴转换：dt=%.4f 小时, 仿真开始=%.1f:00\n', dt, simulation_start_hour);

% --- 【优化】绘图风格配置 ---
font_name = 'Yahei'; 
dpi = 300; 
axis_line_width = 1.0; % 坐标轴线宽 (pt)
grid_color = [0.15 0.15 0.15]; 
grid_alpha = 0.1; 
axis_color = [0.2 0.2 0.2]; % 坐标轴颜色 (深灰)

% 检查字体是否可用
if ~any(strcmpi(listfonts, font_name))
    warning('字体 "%s" 未在您的系统中找到。MATLAB将回退到默认字体。', font_name);
end
fprintf('绘图将使用字体: %s\n', font_name);


% --- 【优化 V2】色盲友好配色方案 ---
% 任务一 (2030年详细图)
color_optimized_up = [0.0, 0.447, 0.741];   % 蓝色
color_optimized_down = [0.850, 0.325, 0.098]; % 橙色
color_demand_and_potential = [0, 0, 0];      % 黑色

% 任务二 (跨年份对比 - 保留)
color_comparison = [
    0, 114, 178;   % 2030年 - 蓝色
    0, 158, 115;   % 2035年 - 绿色
    213, 94, 0     % 2040年 - 橙红色
]/255;

% --- 加载数据 (保留) ---
plotDataAllYears = cell(length(data_files), 1);
data_loaded = false;
for i = 1:length(data_files)
    if exist(data_files{i}, 'file')
        fprintf('  正在加载: %s\n', data_files{i});
        loaded_data = load(data_files{i});
        if isfield(loaded_data, 'plotData')
            plotDataAllYears{i} = loaded_data.plotData;
            data_loaded = true;
        else
            warning('文件 %s 中未找到 "plotData" 结构体。', data_files{i});
        end
    else
        warning('数据文件 %s 未找到。', data_files{i});
    end
end

if ~data_loaded
    error('未能成功加载任何数据。请检查 .mat 文件是否存在。');
end
fprintf('数据加载完毕。\n\n');


%% 2. 绘制任务一：2030年详细调节能力图
fprintf('开始绘制任务一：2030年详细图 (优化风格 V2)...\n');

% --- 2.1 绘制上调能力图 (2030年) ---
try
    fig1 = figure('Color', 'white', 'Position', [100, 100, 900, 450]); 
    ax1 = axes(fig1);
    data2030 = plotDataAllYears{1}.Up; 
    
    % 【保留】计算绝对时间轴
    time_axis_steps = data2030.time_axis;
    time_axis_absolute = (time_axis_steps - 1) * dt + simulation_start_hour;
    sim_end_absolute_hour = time_axis_absolute(end);
    
    hold(ax1, 'on');
    
    % 【V2 优化】图层 1: 系统总潜力 (使用黑色细点线) - 【修改：除以10】
    plot(ax1, time_axis_absolute, data2030.P_potential_total / 10, ...
        'LineWidth', 1.0, ...
        'Color', color_demand_and_potential, ...
        'LineStyle', ':', ... % Dotted line
        'DisplayName', '系统总潜力');
    
    % 【V2 优化】图层 2: 电网需求 (使用黑色中等虚线) - 【修改：除以10】
    plot(ax1, time_axis_absolute, data2030.P_req / 10, ...
        'LineWidth', 1.2, ...
        'Color', color_demand_and_potential, ...
        'LineStyle', '--', ... % Dashed line
        'DisplayName', '电网需求');
        
    % 【V2 优化】图层 3: 优化调度出力 (使用蓝色粗实线，放在最上层) - 【修改：除以10】
    plot(ax1, time_axis_absolute, data2030.p_opt_t / 10, ...
        'LineWidth', 2.0, ...
        'Color', color_optimized_up, ...
        'LineStyle', '-', ...
        'DisplayName', '优化调度出力');
    
    hold(ax1, 'off');
    
    % --- 【优化】设置风格 ---
    xlabel(ax1, '时间', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', font_name); 
    ylabel(ax1, '上调功率 (kW)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', font_name);
    legend(ax1, 'Location', 'northwest', 'Box', 'off', 'FontSize', 12, 'FontName', font_name);
    set(ax1, ...
        'FontName', font_name, ...
        'FontSize', 12, ...
        'Box', 'off', ...
        'TickDir', 'out', ... % 刻度朝外
        'LineWidth', axis_line_width, ...
        'YGrid', 'on', ...
        'XGrid', 'off', ...
        'GridColor', grid_color, ...
        'GridAlpha', grid_alpha, ...
        'XColor', axis_color, ... 
        'YColor', axis_color);
        
    xlim(ax1, [simulation_start_hour, sim_end_absolute_hour]); 
    
    % 【保留】设置X轴时间刻度标签 (6:00, 12:00, 18:00, 0:00(次日), 6:00(次日))
    xticks_to_set = simulation_start_hour:6:sim_end_absolute_hour;
    xtick_labels_to_set = {};
    for tick_hour = xticks_to_set
        hour_of_day_display = mod(tick_hour, 24);
        if tick_hour >= 24 && hour_of_day_display < simulation_start_hour
            xtick_labels_to_set{end+1} = sprintf('%.0f:00 (次日)', hour_of_day_display);
        else
            xtick_labels_to_set{end+1} = sprintf('%.0f:00', hour_of_day_display);
        end
    end
    set(ax1, 'XTick', xticks_to_set, 'XTickLabel', xtick_labels_to_set);
    
    % 保存图像 - 【修改：PNG格式】
    output_filename1 = 'VPP_Up_Regulation_2030_OptimalStyle_V2.png';
    print(fig1, output_filename1, '-dpng', sprintf('-r%d', dpi));
    fprintf('  上调能力图已保存为: %s\n', output_filename1);
catch ME
    fprintf('*** 绘制2030年上调图时出错: %s ***\n', ME.message);
end

% --- 2.2 绘制下调能力图 (2030年) ---
try
    fig2 = figure('Color', 'white', 'Position', [150, 150, 900, 450]);
    ax2 = axes(fig2);
    data2030 = plotDataAllYears{1}.Down; 

    % 【保留】计算绝对时间轴
    time_axis_steps = data2030.time_axis;
    time_axis_absolute = (time_axis_steps - 1) * dt + simulation_start_hour;
    sim_end_absolute_hour = time_axis_absolute(end);
    
    hold(ax2, 'on');
    
    % 【V2 优化】图层 1: 系统总潜力 (使用黑色细点线) - 【修改：除以10】
    plot(ax2, time_axis_absolute, data2030.P_potential_total / 10, ...
        'LineWidth', 1.0, ...
        'Color', color_demand_and_potential, ...
        'LineStyle', ':', ... % Dotted line
        'DisplayName', '系统总潜力');
    
    % 【V2 优化】图层 2: 电网需求 (使用黑色中等虚线) - 【修改：除以10】
    plot(ax2, time_axis_absolute, data2030.P_req / 10, ...
        'LineWidth', 1.2, ...
        'Color', color_demand_and_potential, ...
        'LineStyle', '--', ... % Dashed line
        'DisplayName', '电网需求');
        
    % 【V2 优化】图层 3: 优化调度出力 (使用橙色粗实线) - 【修改：除以10】
    plot(ax2, time_axis_absolute, data2030.p_opt_t / 10, ...
        'LineWidth', 2.0, ...
        'Color', color_optimized_down, ...
        'LineStyle', '-', ...
        'DisplayName', '优化调度出力');
        
    hold(ax2, 'off');
    
    % --- 【优化】设置风格 ---
    xlabel(ax2, '时间', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', font_name); 
    ylabel(ax2, '下调功率 (kW)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', font_name);
    legend(ax2, 'Location', 'northwest', 'Box', 'off', 'FontSize', 12, 'FontName', font_name);
    set(ax2, ...
        'FontName', font_name, ...
        'FontSize', 12, ...
        'Box', 'off', ...
        'TickDir', 'out', ... % 刻度朝外
        'LineWidth', axis_line_width, ...
        'YGrid', 'on', ...
        'XGrid', 'off', ...
        'GridColor', grid_color, ...
        'GridAlpha', grid_alpha, ...
        'XColor', axis_color, ...
        'YColor', axis_color);
        
    xlim(ax2, [simulation_start_hour, sim_end_absolute_hour]); 

    % 【保留】设置X轴时间刻度标签
    xticks_to_set = simulation_start_hour:6:sim_end_absolute_hour;
    xtick_labels_to_set = {};
    for tick_hour = xticks_to_set
        hour_of_day_display = mod(tick_hour, 24);
        if tick_hour >= 24 && hour_of_day_display < simulation_start_hour
            xtick_labels_to_set{end+1} = sprintf('%.0f:00 (次日)', hour_of_day_display);
        else
            xtick_labels_to_set{end+1} = sprintf('%.0f:00', hour_of_day_display);
        end
    end
    set(ax2, 'XTick', xticks_to_set, 'XTickLabel', xtick_labels_to_set);

    % 保存图像 - 【修改：PNG格式】
    output_filename2 = 'VPP_Down_Regulation_2030_OptimalStyle_V2.png';
    print(fig2, output_filename2, '-dpng', sprintf('-r%d', dpi));
    fprintf('  下调能力图已保存为: %s\n\n', output_filename2);
catch ME
     fprintf('*** 绘制2030年下调图时出错: %s ***\n\n', ME.message);
end


%% 3. 绘制任务二：2030-2040年调节能力对比图
fprintf('开始绘制任务二：多年对比图 (优化风格 V2)...\n');

% --- 3.1 绘制上调能力对比图 ---
try
    fig3 = figure('Color', 'white', 'Position', [200, 200, 900, 450]);
    ax3 = axes(fig3);
    hold(ax3, 'on');
    
    legend_entries_up = {};
    time_axis_common_abs = []; 
    
    for i = 1:length(years)
        data = plotDataAllYears{i}.Up;
        
        % 【保留】转换为绝对时间
        time_axis_absolute = (data.time_axis - 1) * dt + simulation_start_hour;
        if isempty(time_axis_common_abs)
            time_axis_common_abs = time_axis_absolute;
            sim_end_absolute_hour = time_axis_common_abs(end); 
        end
        
        % 【优化】线条加粗 - 【修改：除以10】
        plot(ax3, time_axis_absolute, data.p_opt_t / 10, 'LineWidth', 1.8, 'Color', color_comparison(i,:));
        legend_entries_up{end+1} = sprintf('%d年', years(i));
    end
    hold(ax3, 'off');
    
    % --- 【优化】设置风格 ---
    xlabel(ax3, '时间', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', font_name); 
    ylabel(ax3, '优化后上调功率 (kW)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', font_name);
    legend(ax3, legend_entries_up, 'Location', 'northwest', 'Box', 'off', 'FontSize', 12, 'FontName', font_name);
    set(ax3, ...
        'FontName', font_name, ...
        'FontSize', 12, ...
        'Box', 'off', ...
        'TickDir', 'out', ...
        'LineWidth', axis_line_width, ...
        'YGrid', 'on', ...
        'XGrid', 'off', ...
        'GridColor', grid_color, ...
        'GridAlpha', grid_alpha, ...
        'XColor', axis_color, ...
        'YColor', axis_color);
        
    if ~isempty(time_axis_common_abs)
        xlim(ax3, [simulation_start_hour, sim_end_absolute_hour]); 
    end
    
    % 【保留】设置X轴时间刻度标签
    xticks_to_set = simulation_start_hour:6:sim_end_absolute_hour;
    xtick_labels_to_set = {};
    for tick_hour = xticks_to_set
        hour_of_day_display = mod(tick_hour, 24);
        if tick_hour >= 24 && hour_of_day_display < simulation_start_hour
            xtick_labels_to_set{end+1} = sprintf('%.0f:00 (次日)', hour_of_day_display);
        else
            xtick_labels_to_set{end+1} = sprintf('%.0f:00', hour_of_day_display);
        end
    end
    set(ax3, 'XTick', xticks_to_set, 'XTickLabel', xtick_labels_to_set);

    % 保存图像 - 【修改：PNG格式】
    output_filename3 = 'VPP_Up_Comparison_2030-2040_OptimalStyle_V2.png';
    print(fig3, output_filename3, '-dpng', sprintf('-r%d', dpi));
    fprintf('  上调能力对比图已保存为: %s\n', output_filename3);
catch ME
    fprintf('*** 绘制多年上调对比图时出错: %s ***\n', ME.message);
end

% --- 3.2 绘制下调能力对比图 ---
try
    fig4 = figure('Color', 'white', 'Position', [250, 250, 900, 450]);
    ax4 = axes(fig4);
    hold(ax4, 'on');
    
    legend_entries_down = {};
    time_axis_common_abs = []; % 重置
    
    for i = 1:length(years)
        data = plotDataAllYears{i}.Down;
          
        % 【保留】转换为绝对时间
        time_axis_absolute = (data.time_axis - 1) * dt + simulation_start_hour;
        if isempty(time_axis_common_abs)
            time_axis_common_abs = time_axis_absolute;
            sim_end_absolute_hour = time_axis_common_abs(end); 
        end
          
        % 【优化】线条加粗 - 【修改：除以10】
        plot(ax4, time_axis_absolute, data.p_opt_t / 10, 'LineWidth', 1.8, 'Color', color_comparison(i,:));
        legend_entries_down{end+1} = sprintf('%d年', years(i));
    end
    hold(ax4, 'off');
    
    % --- 【优化】设置风格 ---
    xlabel(ax4, '时间', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', font_name); 
    ylabel(ax4, '优化后下调功率 (kW)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', font_name);
    legend(ax4, legend_entries_down, 'Location', 'northwest', 'Box', 'off', 'FontSize', 12, 'FontName', font_name);
    set(ax4, ...
        'FontName', font_name, ...
        'FontSize', 12, ...
        'Box', 'off', ...
        'TickDir', 'out', ...
        'LineWidth', axis_line_width, ...
        'YGrid', 'on', ...
        'XGrid', 'off', ...
        'GridColor', grid_color, ...
        'GridAlpha', grid_alpha, ...
        'XColor', axis_color, ...
        'YColor', axis_color);
        
    if ~isempty(time_axis_common_abs)
        xlim(ax4, [simulation_start_hour, sim_end_absolute_hour]); 
    end
    
    % 【保留】设置X轴时间刻度标签
    xticks_to_set = simulation_start_hour:6:sim_end_absolute_hour;
    xtick_labels_to_set = {};
    for tick_hour = xticks_to_set
        hour_of_day_display = mod(tick_hour, 24);
        if tick_hour >= 24 && hour_of_day_display < simulation_start_hour
            xtick_labels_to_set{end+1} = sprintf('%.0f:00 (次日)', hour_of_day_display);
        else
            xtick_labels_to_set{end+1} = sprintf('%.0f:00', hour_of_day_display);
        end
    end
    set(ax4, 'XTick', xticks_to_set, 'XTickLabel', xtick_labels_to_set);
    
    % 保存图像 - 【修改：PNG格式】
    output_filename4 = 'VPP_Down_Comparison_2030-2040_OptimalStyle_V2.png';
    print(fig4, output_filename4, '-dpng', sprintf('-r%d', dpi));
    fprintf('  下调能力对比图已保存为: %s\n', output_filename4);
catch ME
    fprintf('*** 绘制多年下调对比图时出错: %s ***\n', ME.message);
end

fprintf('\n所有绘图任务完成。\n');