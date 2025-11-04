%% 主程序：使用遗传算法优化并检验互补性指标
clc; close all; clear;

%% 步骤1：定义参数
load('results_single.mat');
T = length(results.AC_Down);              % 时间段数
n_AC_max = 100;                           % 空调的最大数量
n_EV_max = 100;                           % 电动汽车的最大数量
dt = 0.05; 
steps_per_hour = 1 / dt;
% 原始调节能力（未归一化）
AC_Down_raw = abs(results.AC_Down);  
AC_Up_raw = results.AC_Up;         
EV_Down_raw = abs(results.EV_Down);  
EV_Up_raw = results.EV_Up;         

% 归一化处理（用于优化）
AC_Down = AC_Down_raw / n_AC_max;  
AC_Up = AC_Up_raw / n_AC_max;      
EV_Down = EV_Down_raw / n_EV_max;  
EV_Up = EV_Up_raw / n_EV_max; 

AC_Up_ind = abs(results.AC_Up_Individual);
AC_Down_ind = abs(results.AC_Down_Individual);
EV_Up_ind = abs(results.EV_Up_Individual);
EV_Down_ind = abs(results.EV_Down_Individual);
% 关键参数调整
theta_1 = 0.2;                            % 空调互补性阈值
theta_2 = 0.2;                           % 电动汽车互补性阈值
theta_3 = -0.3;                           % 相关系数阈值
lambda_1 = 0.5;                           % SDCI惩罚系数
lambda_2 = 1;                           % ρ惩罚系数

%% 步骤2：设置遗传算法参数
nvars = 2 * T;                            % 决策变量数
lb = zeros(1, nvars);                     % 下界
ub = [n_AC_max * ones(1, T), n_EV_max * ones(1, T)];  % 上界
intcon = 1:2*T;                           % 整数变量位置

options = optimoptions('ga', ...
    'PopulationSize', 50, ...
    'MaxGenerations', 500, ...
    'Display', 'iter', ...
    'PlotFcn', @gaplotbestf, ...
    'ConstraintTolerance', 1e-6);

%% 步骤3：运行优化算法
% [x_opt_up, fval_up] = ga(@(x) objectiveFunction_up_aver(x, n_AC_max, n_EV_max, AC_Up, EV_Up, T, ...
%                         theta_1, theta_3, lambda_1, lambda_2), ...
%                         nvars, [], [], [], [], lb, ub, [], intcon, options);
% 
% [x_opt_down, fval_down] = ga(@(x) objectiveFunction_down_aver(x, n_AC_max, n_EV_max, AC_Down, EV_Down, T, ...
%                         theta_2, theta_3, lambda_1, lambda_2), ...
%                         nvars, [], [], [], [], lb, ub, [], intcon, options);
% [x_opt_up, fval_up] = ga(@(x) objectiveFunction_up(x, n_AC_max, n_EV_max, AC_Up_ind, EV_Up_ind, T, ...
%                         theta_1, theta_3, lambda_1, lambda_2), ...
%                         nvars, [], [], [], [], lb, ub, [], intcon, options);
% 
% [x_opt_down, fval_down] = ga(@(x) objectiveFunction_down(x, n_AC_max, n_EV_max, AC_Down_ind, EV_Down_ind, T, ...
%                         theta_2, theta_3, lambda_1, lambda_2), ...
%                         nvars, [], [], [], [], lb, ub, [], intcon, options);
[x_opt_up, fval_up] = ga(@(x) objectiveFunction_up_real_time(x, n_AC_max, n_EV_max, AC_Up_ind, EV_Up_ind, T-1, steps_per_hour, theta_1, theta_3, lambda_1, lambda_2), ...
                        nvars, [], [], [], [], lb, ub, [], intcon, options);

[x_opt_down, fval_down] = ga(@(x) objectiveFunction_real_time(x, n_AC_max, n_EV_max, AC_Down_ind, EV_Down_ind, T-1, steps_per_hour, theta_2, theta_3, lambda_1, lambda_2), ...
                        nvars, [], [], [], [], lb, ub, [], intcon, options);

%% 步骤4：计算调节能力
n_AC_up = x_opt_up(1:T);
n_EV_up = x_opt_up(T+1:end);
n_AC_down = x_opt_down(1:T);
n_EV_down = x_opt_down(T+1:end);

AC_Up_opt = AC_Up_raw .* n_AC_up' / n_AC_max;
EV_Up_opt = EV_Up_raw .* n_EV_up' / n_EV_max;
AC_Down_opt = AC_Down_raw .* n_AC_down' / n_AC_max;
EV_Down_opt = EV_Down_raw .* n_EV_down' / n_EV_max;

Total_Up_raw = AC_Up_raw + EV_Up_raw;
Total_Down_raw = AC_Down_raw + EV_Down_raw;
Total_Up_opt = AC_Up_opt + EV_Up_opt;
Total_Down_opt = AC_Down_opt + EV_Down_opt;

%% 步骤5：检验互补性指标（确保输出为标量）
n_AC_all = n_AC_max * ones(T,1);
n_EV_all = n_EV_max * ones(T,1);

% 上调指标
SDCI_up_raw = ensureScalar(calculateSDCI(n_AC_all, n_EV_all, AC_Up, EV_Up));
rho_up_raw = ensureScalar(calculateSpearmanRho(n_AC_all, AC_Up, n_EV_all, EV_Up));
SDCI_up_opt = ensureScalar(calculateSDCI(n_AC_up, n_EV_up, AC_Up, EV_Up));
rho_up_opt = ensureScalar(calculateSpearmanRho(n_AC_up, AC_Up, n_EV_up, EV_Up));

% 下调指标（强制负相关）
SDCI_down_raw = ensureScalar(calculateSDCI(n_AC_all, n_EV_all, AC_Down, EV_Down));
rho_down_raw = -abs(ensureScalar(calculateSpearmanRho(n_AC_all, AC_Down, n_EV_all, EV_Down)));
SDCI_down_opt = ensureScalar(calculateSDCI(n_AC_down, n_EV_down, AC_Down, EV_Down));
rho_down_opt = -abs(ensureScalar(calculateSpearmanRho(n_AC_down, AC_Down, n_EV_down, EV_Down)));

%% 步骤6：绘制综合对比图
time_points = 1:T;
figure('Position', [100, 100, 1200, 900])

% 子图1：上调能力对比
subplot(2,2,1)
plot(time_points, Total_Up_raw, 'k-', 'LineWidth', 2); hold on;
plot(time_points, Total_Up_opt, 'r--', 'LineWidth', 2);
title('总上调能力对比');
xlabel('时间点'); ylabel('功率 (kW)');
legend('原始总上调', '优化后总上调', 'Location', 'best');
grid on;

% 子图2：下调能力对比
subplot(2,2,2)
plot(time_points, Total_Down_raw, 'k-', 'LineWidth', 2); hold on;
plot(time_points, Total_Down_opt, 'b--', 'LineWidth', 2);
title('总下调能力对比');
xlabel('时间点'); ylabel('功率 (kW)');
legend('原始总下调', '优化后总下调', 'Location', 'best');
grid on;

% 子图3：上调指标对比（修正维度问题）
subplot(2,2,3)
bar_data = [SDCI_up_raw, SDCI_up_opt; rho_up_raw, rho_up_opt];
bar(bar_data');
set(gca, 'XTickLabel', {'SDCI⁺', 'ρ'});
legend('优化前', '优化后', 'Location', 'northoutside');
title('上调互补性指标对比');
ylabel('指标值');
grid on;

% 子图4：下调指标对比（修正维度问题）
subplot(2,2,4)
bar_data = [SDCI_down_raw, SDCI_down_opt; rho_down_raw, rho_down_opt];
bar(bar_data');
set(gca, 'XTickLabel', {'SDCI⁻', 'ρ'});
legend('优化前', '优化后', 'Location', 'northoutside');
title('下调互补性指标对比（ρ强制负相关）');
ylabel('指标值');
grid on;

%% 步骤7：显示优化结果
disp('=== 优化结果汇总 ===');
disp(['上调节最优目标值: ', num2str(fval_up)]);
disp(['下调节最优目标值: ', num2str(fval_down)]);
disp(['空调平均参与数量(上/下): ', num2str(mean(n_AC_up)), ' / ', num2str(mean(n_AC_down))]);
disp(['EV平均参与数量(上/下): ', num2str(mean(n_EV_up)), ' / ', num2str(mean(n_EV_down))]);

disp('=== 互补性指标检验 ===');
disp('【上调】');
disp(['原始 SDCI⁺: ', num2str(SDCI_up_raw), ' | 优化后: ', num2str(SDCI_up_opt)]);
disp(['原始 ρ: ', num2str(rho_up_raw), ' | 优化后: ', num2str(rho_up_opt)]);
disp('【下调】');
disp(['原始 SDCI⁻: ', num2str(SDCI_down_raw), ' | 优化后: ', num2str(SDCI_down_opt)]);
disp(['原始 ρ: ', num2str(rho_down_raw), ' | 优化后: ', num2str(rho_down_opt), ' (强制负相关)']);

%% 辅助函数：确保输出为标量
function val = ensureScalar(inputVal)
    if isscalar(inputVal)
        val = inputVal;
    else
        val = mean(inputVal(:)); % 对非标量取均值
    end
end
