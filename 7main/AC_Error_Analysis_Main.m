%% AC_Error_Analysis_Main.m
% =========================================================================
% 主程序：空调聚合体功率误差分析
% 功能：
%   1. 读取仿真结果 (AC_Stateful_Simulation_Results_5min_pi.mat)
%   2. 提取“单体累加实际功率”与“聚合模型推算功率”
%   3. 计算各项误差指标 (RMSE, MAE, Max Error)
%   4. 绘制对比曲线和误差分布图
% =========================================================================

clear; close all; clc;

%% 1. 设置与加载
result_file = 'AC_Stateful_Simulation_Results_5min_pi.mat';
fprintf('正在加载仿真结果文件: %s ...\n', result_file);

if ~exist(result_file, 'file')
    error('未找到结果文件 %s。\n请先运行主仿真程序 AC_main_1_inc_pi.m 生成数据。', result_file);
end

data = load(result_file);
if ~isfield(data, 'results')
    error('数据文件结构不正确，缺失 "results" 字段。');
end
res = data.results;

%% 2. 提取关键数据
% 检查是否存在必要的字段
if ~isfield(res, 'Agg_Total_Power') || ~isfield(res, 'Agg_Model_Total_Power')
    error('数据中缺失功率数据 (Agg_Total_Power 或 Agg_Model_Total_Power)。');
end

% 提取数据向量
P_physics = res.Agg_Total_Power;       % 单体累加功率 (物理真值)
P_model   = res.Agg_Model_Total_Power; % 聚合模型功率 (理论推算值)
Time      = res.time_points;           % 时间轴

% 确保数据维度一致
P_physics = P_physics(:);
P_model = P_model(:);
Time = Time(:);

if length(P_physics) ~= length(P_model)
    warning('物理功率与模型功率数据长度不一致，将截断至最短长度进行比较。');
    min_len = min(length(P_physics), length(P_model));
    P_physics = P_physics(1:min_len);
    P_model = P_model(1:min_len);
    Time = Time(1:min_len);
end

%% 3. 计算误差指标
fprintf('\n正在计算误差指标...\n');

% 3.1 绝对误差向量 (Difference)
Error_Abs_Vec = P_physics - P_model;  % 物理值 - 模型值

% 3.2 相对误差 (避免分母为0)
% 设定一个微小的阈值防止除零
valid_idx = abs(P_physics) > 1e-3; 
Error_Rel_Vec = zeros(size(P_physics));
Error_Rel_Vec(valid_idx) = abs(Error_Abs_Vec(valid_idx)) ./ abs(P_physics(valid_idx)) * 100; % 百分比

% 3.3 统计指标
MAE  = mean(abs(Error_Abs_Vec));              % 平均绝对误差
RMSE = sqrt(mean(Error_Abs_Vec.^2));          % 均方根误差
[Max_Err_Val, Max_Err_Idx] = max(abs(Error_Abs_Vec)); % 最大绝对误差
Max_Err_Time = Time(Max_Err_Idx);             % 最大误差发生时刻

% 3.4 打印报告
fprintf('======================================================\n');
fprintf('             空调聚合模型误差分析报告             \n');
fprintf('======================================================\n');
fprintf('1. 数据点总数       : %d\n', length(P_physics));
fprintf('2. 平均绝对误差(MAE): %.4f kW\n', MAE);
fprintf('3. 均方根误差 (RMSE): %.4f kW\n', RMSE);
fprintf('4. 最大绝对误差     : %.4f kW (发生于 %.2f 小时)\n', Max_Err_Val, Max_Err_Time);
fprintf('5. 平均相对误差     : %.2f %%\n', mean(Error_Rel_Vec));
fprintf('------------------------------------------------------\n');
fprintf('物理功率范围        : [%.2f, %.2f] kW\n', min(P_physics), max(P_physics));
fprintf('模型功率范围        : [%.2f, %.2f] kW\n', min(P_model), max(P_model));
fprintf('======================================================\n');

%% 4. 绘图分析

% 设置字体
set(0, 'DefaultAxesFontName', 'Microsoft YaHei'); 
set(0, 'DefaultTextFontName', 'Microsoft YaHei');

figure('Name', '聚合模型误差分析', 'Position', [100, 100, 1000, 700], 'Color', 'w');

% --- 子图 1: 功率对比 ---
subplot(2, 1, 1);
hold on;
plot(Time, P_physics, 'r-', 'LineWidth', 1.5, 'DisplayName', '单体累加功率 (物理真值)');
plot(Time, P_model, 'g--', 'LineWidth', 1.5, 'DisplayName', '聚合模型功率 (理论推算)');
% 标记最大误差处
x_max = Max_Err_Time;
y_phy = P_physics(Max_Err_Idx);
y_mod = P_model(Max_Err_Idx);
plot([x_max, x_max], [y_phy, y_mod], 'k-', 'LineWidth', 1); % 误差连线
plot(x_max, y_phy, 'ro', 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
plot(x_max, y_mod, 'go', 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
text(x_max + 0.5, max(y_phy, y_mod), sprintf('Max Err: %.2f kW', Max_Err_Val), 'FontSize', 10);

hold off;
title('聚合体运行功率对比: 物理仿真 vs 数学模型');
ylabel('功率 (kW)');
legend('Location', 'best');
grid on;
xlim([0, 24]);

% --- 子图 2: 误差曲线 ---
subplot(2, 1, 2);
plot(Time, Error_Abs_Vec, 'k-', 'LineWidth', 1.2);
hold on;
yline(0, 'r--', 'LineWidth', 1);
yline(MAE, 'b:', 'LineWidth', 1.5, 'DisplayName', 'MAE');
yline(-MAE, 'b:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
hold off;
title('绝对误差曲线 (物理值 - 模型值)');
xlabel('时间 (小时)');
ylabel('功率误差 (kW)');
legend({'误差曲线', '零线', sprintf('MAE=%.2f', MAE)}, 'Location', 'best');
grid on;
xlim([0, 24]);

% 保存图片
save_img_name = 'AC_Model_Error_Analysis.png';
print(gcf, save_img_name, '-dpng', '-r300');
fprintf('分析图表已保存至: %s\n', save_img_name);

%% 5. (可选) 保存误差数据
save('AC_Error_Analysis_Data.mat', 'Error_Abs_Vec', 'MAE', 'RMSE', 'Max_Err_Val');
fprintf('误差数据已保存至: AC_Error_Analysis_Data.mat\n');