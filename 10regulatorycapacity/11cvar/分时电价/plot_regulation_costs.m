%% plot_regulation_costs.m
% 功能：绘制 AC、EV、火电 (Gen) 的分时调节成本电价曲线
% 逻辑来源：test_ieee30_func_8ampv_tly.m
% 输出：高 DPI PNG 图片，无标题，中文标注
% 时间轴：08:00 - 次日 08:00

clear; close all; clc;

%% 1. 初始化时间轴参数
% 仿真时间设置：从第8小时开始，持续24小时
simulation_start_hour = 8;
simulation_end_hour   = 32; % 次日 8 点
dt = 15/60;                 % 时间步长 15分钟 (0.25小时)

% 生成时间序列
t_axis = simulation_start_hour : dt : simulation_end_hour;
T_steps = length(t_axis);

% 生成用于判断的小时向量 (0-23循环)
Hour_Vector = mod(t_axis, 24)'; 

%% 2. 复现分时成本逻辑 (基于 test_ieee30_func_8ampv_tly.m)

% --- 火电 (Gen) ---
% 基准: 800
% 高峰 (10-12, 18-21): 1000
% 低谷 (23-5): 600
C1_Gen_Vec = 800 * ones(T_steps, 1); 
idx_gen_peak   = (Hour_Vector >= 10 & Hour_Vector < 12) | ...
                 (Hour_Vector >= 18 & Hour_Vector < 21);
idx_gen_valley = (Hour_Vector >= 23) | (Hour_Vector < 5);
C1_Gen_Vec(idx_gen_peak)   = 1000;
C1_Gen_Vec(idx_gen_valley) = 600;

% --- 电动汽车 (EV) ---
% 基准: 500
% 紧迫 (17-21): 800
% 低谷 (23-7): 300
C1_EV_Vec = 500 * ones(T_steps, 1); 
idx_ev_urgent = (Hour_Vector >= 17 & Hour_Vector < 21);
idx_ev_low    = (Hour_Vector >= 23) | (Hour_Vector < 7);
C1_EV_Vec(idx_ev_urgent) = 800;
C1_EV_Vec(idx_ev_low)    = 300;

% --- 空调 (AC) ---
% 基准: 400
% 高温 (13-16): 600
% 低谷 (22-8): 300
C1_AC_Vec = 400 * ones(T_steps, 1); 
idx_ac_hot  = (Hour_Vector >= 13 & Hour_Vector < 16);
idx_ac_cool = (Hour_Vector >= 22) | (Hour_Vector < 8);
C1_AC_Vec(idx_ac_hot)  = 600;
C1_AC_Vec(idx_ac_cool) = 300;

%% 3. 绘图设置
% 设置字体以支持中文
font_name = 'Microsoft YaHei'; 
font_size_axis = 14;
font_size_label = 16;
font_size_legend = 14;
line_width = 2.5;

figure('Color', 'w', 'Position', [100, 100, 1000, 600]);
hold on;

% 绘制阶梯图或折线图 (此处使用 plot 模拟连续变化，符合原逻辑中的向量定义)
% 为了视觉区分，使用不同的线型和颜色
p1 = plot(t_axis, C1_Gen_Vec, 'r-', 'LineWidth', line_width, 'DisplayName', '火电机组 (Gen)');
p2 = plot(t_axis, C1_EV_Vec,  'g--', 'LineWidth', line_width, 'DisplayName', '电动汽车 (EV)');
p3 = plot(t_axis, C1_AC_Vec,  'b-.', 'LineWidth', line_width, 'DisplayName', '空调集群 (AC)');

% 网格与边框
grid on;
box on;
ax = gca;
ax.LineWidth = 1.5;
ax.FontName = font_name;
ax.FontSize = font_size_axis;

% --- 坐标轴设置 ---
% X轴范围
xlim([simulation_start_hour, simulation_end_hour]);

% 自定义 X 轴刻度标签 (08:00 到 次日 08:00)
x_ticks = [8, 12, 16, 20, 24, 28, 32];
x_labels = {'08:00', '12:00', '16:00', '20:00', '00:00 (次日)', '04:00 (次日)', '08:00 (次日)'};
set(gca, 'XTick', x_ticks, 'XTickLabel', x_labels);

% Y轴范围微调 (留出一点上下边距)
ylim([0, 1200]); 
ylabel('调节成本 (元/MW)', 'FontName', font_name, 'FontSize', font_size_label);
xlabel('时间', 'FontName', font_name, 'FontSize', font_size_label);

% --- 图例设置 ---
legend([p1, p2, p3], 'Location', 'north', 'Orientation', 'horizontal', ...
       'FontSize', font_size_legend, 'FontName', font_name, 'Box', 'off');

%% 4. 保存为高DPI图片
output_filename = '分时调节成本电价.png';
fprintf('正在保存图片到: %s ...\n', output_filename);

% 使用 print 函数保存为 600 DPI 的 PNG
print(gcf, output_filename, '-dpng', '-r600');

fprintf('完成。\n');