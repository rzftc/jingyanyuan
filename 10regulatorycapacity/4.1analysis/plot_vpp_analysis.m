% plot_vpp_analysis_final.m
%
% 脚本功能：
% 1. 加载 2030, 2035, 2040 年的VPP优化数据。
% 2. 绘制2030年详细图 (任务一)。
% 3. 绘制2030-2040年趋势对比图 (任务二)。
% 4. 严格遵循学术期刊风格指南：
%    - 字体: 'Microsoft YaHei' (中文设为 "yahei")
%    - 线条: 0.8 pt
%    - 刻度: 向内 (TickDir = 'in')
%    - 网格: 仅 Y 轴, 浅色, 10% 透明
%    - 标题: 无
%    - 图例: 无边框 (Box = 'off')
%    - 保存: 300 DPI PNG

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

% --- 绘图风格配置 ---
font_name = 'Microsoft YaHei'; % '微软雅黑'
dpi = 300; % 保存图像的分辨率
line_width = 0.8; % 曲线线条粗细 (pt)
axis_line_width = 0.8; % 坐标轴粗细 (pt)
grid_color = [0.15 0.15 0.15]; % 网格线颜色
grid_alpha = 0.1; % 网格透明度 (10%)

% 检查字体是否可用
if ~any(strcmpi(listfonts, font_name))
    warning('字体 "%s" (YaHei) 未在您的系统中找到。MATLAB将尝试使用此设置，但如果失败，将回退到系统默认字体。中文可能无法正常显示。', font_name);
end
fprintf('绘图将强制尝试使用字体: %s\n', font_name);


% --- 色盲友好配色方案 ---
% 任务一 (2030年详细图)
color_demand = [0 0 0];       % 黑色 (电网需求)
color_potential = [0 0 0];    % 黑色 (系统总潜力)
color_optimized_up = [0, 114, 178]/255;   % 蓝色
color_optimized_down = [213, 94, 0]/255; % 橙红色

% 任务二 (跨年份对比)
color_comparison = [
    0, 114, 178;   % 2030年 - 蓝色
    0, 158, 115;   % 2035年 - 绿色
    213, 94, 0     % 2040年 - 橙红色
]/255;

% --- 加载数据 ---
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
fprintf('开始绘制任务一：2030年详细图...\n');

% --- 2.1 绘制上调能力图 (2030年) ---
try
    fig1 = figure('Color', 'white', 'Position', [100, 100, 900, 450]); % 宽幅
    ax1 = axes(fig1);
    data2030 = plotDataAllYears{1}.Up; % 获取2030年上调数据
    
    hold(ax1, 'on');
    plot(ax1, data2030.time_axis, data2030.P_potential_total, ...
        'LineWidth', line_width*0.8, 'Color', color_potential, 'LineStyle', ':', 'DisplayName', '系统总潜力');
    plot(ax1, data2030.time_axis, data2030.P_req, ...
        'LineWidth', line_width, 'Color', color_demand, 'LineStyle', '--', 'DisplayName', '电网需求');
    plot(ax1, data2030.time_axis, data2030.p_opt_t, ...
        'LineWidth', line_width, 'Color', color_optimized_up, 'LineStyle', '-', 'DisplayName', '优化调度出力');
    hold(ax1, 'off');
    
    % --- 设置风格 (无标题) ---
    xlabel(ax1, '时间步', 'FontSize', 12, 'FontName', font_name);
    ylabel(ax1, '上调功率 (kW)', 'FontSize', 12, 'FontName', font_name);
    legend(ax1, 'Location', 'northeast', 'Box', 'off', 'FontSize', 10, 'FontName', font_name);
    set(ax1, ...
        'FontName', font_name, ...
        'FontSize', 10, ...
        'Box', 'off', ...
        'TickDir', 'in', ...
        'LineWidth', axis_line_width, ...
        'YGrid', 'on', ...
        'XGrid', 'off', ...
        'GridColor', grid_color, ...
        'GridAlpha', grid_alpha);
    xlim(ax1, [data2030.time_axis(1), data2030.time_axis(end)]);
    
    % 保存图像
    output_filename1 = 'VPP_Up_Regulation_2030_Nature.png';
    print(fig1, output_filename1, '-dpng', sprintf('-r%d', dpi));
    fprintf('  上调能力图已保存为: %s\n', output_filename1);
catch ME
    fprintf('*** 绘制2030年上调图时出错: %s ***\n', ME.message);
end

% --- 2.2 绘制下调能力图 (2030年) ---
try
    fig2 = figure('Color', 'white', 'Position', [150, 150, 900, 450]);
    ax2 = axes(fig2);
    data2030 = plotDataAllYears{1}.Down; % 获取2030年下调数据
    
    hold(ax2, 'on');
    plot(ax2, data2030.time_axis, data2030.P_potential_total, ...
        'LineWidth', line_width*0.8, 'Color', color_potential, 'LineStyle', ':', 'DisplayName', '系统总潜力');
    plot(ax2, data2030.time_axis, data2030.P_req, ...
        'LineWidth', line_width, 'Color', color_demand, 'LineStyle', '--', 'DisplayName', '电网需求');
    plot(ax2, data2030.time_axis, data2030.p_opt_t, ...
        'LineWidth', line_width, 'Color', color_optimized_down, 'LineStyle', '-', 'DisplayName', '优化调度出力');
    hold(ax2, 'off');
    
    % --- 设置风格 (无标题) ---
    xlabel(ax2, '时间步', 'FontSize', 12, 'FontName', font_name);
    ylabel(ax2, '下调功率 (kW)', 'FontSize', 12, 'FontName', font_name);
    legend(ax2, 'Location', 'northeast', 'Box', 'off', 'FontSize', 10, 'FontName', font_name);
    set(ax2, ...
        'FontName', font_name, ...
        'FontSize', 10, ...
        'Box', 'off', ...
        'TickDir', 'in', ...
        'LineWidth', axis_line_width, ...
        'YGrid', 'on', ...
        'XGrid', 'off', ...
        'GridColor', grid_color, ...
        'GridAlpha', grid_alpha);
    xlim(ax2, [data2030.time_axis(1), data2030.time_axis(end)]);

    % 保存图像
    output_filename2 = 'VPP_Down_Regulation_2030_Nature.png';
    print(fig2, output_filename2, '-dpng', sprintf('-r%d', dpi));
    fprintf('  下调能力图已保存为: %s\n\n', output_filename2);
catch ME
     fprintf('*** 绘制2030年下调图时出错: %s ***\n\n', ME.message);
end


%% 3. 绘制任务二：2030-2040年调节能力对比图
fprintf('开始绘制任务二：多年对比图...\n');

% --- 3.1 绘制上调能力对比图 ---
try
    fig3 = figure('Color', 'white', 'Position', [200, 200, 900, 450]);
    ax3 = axes(fig3);
    hold(ax3, 'on');
    
    legend_entries_up = {};
    time_axis_common = []; 
    
    for i = 1:length(years)
        data = plotDataAllYears{i}.Up;
        if isempty(time_axis_common)
            time_axis_common = data.time_axis;
        end
        plot(ax3, data.time_axis, data.p_opt_t, 'LineWidth', line_width, 'Color', color_comparison(i,:));
        legend_entries_up{end+1} = sprintf('%d年', years(i));
    end
    hold(ax3, 'off');
    
    % --- 设置风格 (无标题) ---
    xlabel(ax3, '时间步', 'FontSize', 12, 'FontName', font_name);
    ylabel(ax3, '优化后上调功率 (kW)', 'FontSize', 12, 'FontName', font_name);
    legend(ax3, legend_entries_up, 'Location', 'northeast', 'Box', 'off', 'FontSize', 10, 'FontName', font_name);
    set(ax3, ...
        'FontName', font_name, ...
        'FontSize', 10, ...
        'Box', 'off', ...
        'TickDir', 'in', ...
        'LineWidth', axis_line_width, ...
        'YGrid', 'on', ...
        'XGrid', 'off', ...
        'GridColor', grid_color, ...
        'GridAlpha', grid_alpha);
    if ~isempty(time_axis_common)
        xlim(ax3, [time_axis_common(1), time_axis_common(end)]);
    end
    
    % 保存图像
    output_filename3 = 'VPP_Up_Comparison_2030-2040_Nature.png';
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
    time_axis_common = []; % 重置
    
    for i = 1:length(years)
        data = plotDataAllYears{i}.Down;
         if isempty(time_axis_common)
            time_axis_common = data.time_axis;
         end
        plot(ax4, data.time_axis, data.p_opt_t, 'LineWidth', line_width, 'Color', color_comparison(i,:));
        legend_entries_down{end+1} = sprintf('%d年', years(i));
    end
    hold(ax4, 'off');
    
    % --- 设置风格 (无标题) ---
    xlabel(ax4, '时间步', 'FontSize', 12, 'FontName', font_name);
    ylabel(ax4, '优化后下调功率 (kW)', 'FontSize', 12, 'FontName', font_name);
    legend(ax4, legend_entries_down, 'Location', 'northeast', 'Box', 'off', 'FontSize', 10, 'FontName', font_name);
    set(ax4, ...
        'FontName', font_name, ...
        'FontSize', 10, ...
        'Box', 'off', ...
        'TickDir', 'in', ...
        'LineWidth', axis_line_width, ...
        'YGrid', 'on', ...
        'XGrid', 'off', ...
        'GridColor', grid_color, ...
        'GridAlpha', grid_alpha);
    if ~isempty(time_axis_common)
        xlim(ax4, [time_axis_common(1), time_axis_common(end)]);
    end
    
    % 保存图像
    output_filename4 = 'VPP_Down_Comparison_2030-2040_Nature.png';
    print(fig4, output_filename4, '-dpng', sprintf('-r%d', dpi));
    fprintf('  下调能力对比图已保存为: %s\n', output_filename4);
catch ME
    fprintf('*** 绘制多年下调对比图时出错: %s ***\n', ME.message);
end

fprintf('\n所有绘图任务完成。\n');