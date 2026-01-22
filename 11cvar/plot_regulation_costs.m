%% plot_regulation_costs.m
% 功能：绘制 AC、EV、火电 (Gen) 的分时调节成本电价曲线
% 输出：高 DPI PNG 图片，无标题，中文标注
% 时间轴：08:00 - 次日 08:00
clear; close all; clc;
%% 1. 初始化时间轴参数
simulation_start_hour = 8;
simulation_end_hour   = 32; % 次日 8 点
dt = 15/60;
t_axis = simulation_start_hour : dt : simulation_end_hour;
T_steps = length(t_axis);
Hour_Vector = mod(t_axis, 24)'; 
%% 2. 复现分时成本逻辑
% --- 火电 (Gen) ---
C1_Gen_Vec = 800 * ones(T_steps, 1); 
idx_gen_peak   = (Hour_Vector >= 10 & Hour_Vector < 12) | ...
                 (Hour_Vector >= 18 & Hour_Vector < 21);
idx_gen_valley = (Hour_Vector >= 23) | (Hour_Vector < 5);
C1_Gen_Vec(idx_gen_peak)   = 1000;
C1_Gen_Vec(idx_gen_valley) = 600;
% --- 电动汽车 (EV) ---
C1_EV_Vec = 500 * ones(T_steps, 1); 
idx_ev_urgent = (Hour_Vector >= 17 & Hour_Vector < 21);
idx_ev_low    = (Hour_Vector >= 23) | (Hour_Vector < 7);
C1_EV_Vec(idx_ev_urgent) = 800;
C1_EV_Vec(idx_ev_low)    = 300;
% --- 空调 (AC) ---
C1_AC_Vec = 400 * ones(T_steps, 1); 
idx_ac_hot  = (Hour_Vector >= 13 & Hour_Vector < 16);
idx_ac_cool = (Hour_Vector >= 22) | (Hour_Vector < 8);
C1_AC_Vec(idx_ac_hot)  = 600;
C1_AC_Vec(idx_ac_cool) = 300;
%% 3. 绘图设置
font_name = 'Microsoft YaHei'; 
% 修改 1: 适当增大字号
font_size_axis = 16; 
font_size_label = 18;
font_size_legend = 16;
line_width = 2.5;
figure('Color', 'w', 'Position', [100, 100, 1000, 600]);
hold on;
p1 = plot(t_axis, C1_Gen_Vec, 'r-', 'LineWidth', line_width, 'DisplayName', '火电机组 (Gen)');
p2 = plot(t_axis, C1_EV_Vec,  'g--', 'LineWidth', line_width, 'DisplayName', '电动汽车 (EV)');
p3 = plot(t_axis, C1_AC_Vec,  'b-.', 'LineWidth', line_width, 'DisplayName', '空调集群 (AC)');
% 修改 2: 去掉网格线
grid off;
box on;
ax = gca;
ax.LineWidth = 1.5;
ax.FontName = font_name;
ax.FontSize = font_size_axis;
% --- 坐标轴设置 ---
xlim([simulation_start_hour, simulation_end_hour]);
x_ticks = [8, 12, 16, 20, 24, 28, 32];
x_labels = {'08:00', '12:00', '16:00', '20:00', '00:00 (次日)', '04:00 (次日)', '08:00 (次日)'};
set(gca, 'XTick', x_ticks, 'XTickLabel', x_labels);
ylim([0, 1200]); 
ylabel('调节成本 (元/MW)', 'FontName', font_name, 'FontSize', font_size_label);
xlabel('时间', 'FontName', font_name, 'FontSize', font_size_label);
% --- 图例设置 ---
legend([p1, p2, p3], 'Location', 'north', 'Orientation', 'horizontal', ...
       'FontSize', font_size_legend, 'FontName', font_name, 'Box', 'off');
%% 4. 保存为高DPI图片
output_filename = '分时调节成本电价.png';
% print 函数保存为 600 DPI 的 PNG
print(gcf, output_filename, '-dpng', '-r600');