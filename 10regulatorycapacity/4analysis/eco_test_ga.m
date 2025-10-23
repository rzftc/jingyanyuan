%% 测试脚本 (已修改为调用 solve_total_time_dispatch_ga - 精简功能版)
clc; close all; clear;

% 假设 solve_total_time_dispatch_ga.m, calculateSDCI.m, calculateSpearmanRho.m 在MATLAB路径中

%% 1. 加载数据与核心参数提取
load('results_single.mat'); % 假设此mat文件包含所有必要的个体数据和可能的电网需求

% 直接从 results 结构体中提取个体数据 (假设存在且格式正确)
P_ac_up_individual = results.AC_Up_Individual;
P_ac_down_individual = abs(results.AC_Down_Individual); % 确保为正
P_ev_up_individual = results.EV_Up_Individual;
P_ev_down_individual = abs(results.EV_Down_Individual); % 确保为正

num_ac_total = size(P_ac_up_individual, 1);
num_ev_total = size(P_ev_up_individual, 1);
T = size(P_ac_up_individual, 2); % 假设所有个体数据的时间维度一致

% 定义/加载电网调节需求
% 如果 results_single.mat 中没有提供，则生成示例需求
if isfield(results, 'P_grid_up_regulation_demand') && ~isempty(results.P_grid_up_regulation_demand) && length(results.P_grid_up_regulation_demand) == T
    P_grid_up_demand = results.P_grid_up_regulation_demand(:);
else
    avg_raw_total_up_power_per_step = 0;
    if T > 0
        avg_raw_total_up_power_per_step = (sum(sum(P_ac_up_individual,1)) + sum(sum(P_ev_up_individual,1))) / T;
    end
    if isnan(avg_raw_total_up_power_per_step) || avg_raw_total_up_power_per_step == 0; avg_raw_total_up_power_per_step = 100; end % 默认基础值
    demand_factor_up_low = 0.2; demand_factor_up_high = 0.5;
    P_grid_up_demand = avg_raw_total_up_power_per_step * (demand_factor_up_low + (demand_factor_up_high-demand_factor_up_low) * rand(T,1));
    % results.P_grid_up_regulation_demand = P_grid_up_demand; % 可选：如果希望在当前工作区更新results结构
end

if isfield(results, 'P_grid_down_regulation_demand') && ~isempty(results.P_grid_down_regulation_demand) && length(results.P_grid_down_regulation_demand) == T
    P_grid_down_demand = results.P_grid_down_regulation_demand(:);
else
    avg_raw_total_down_power_per_step = 0;
    if T > 0
        avg_raw_total_down_power_per_step = (sum(sum(P_ac_down_individual,1)) + sum(sum(P_ev_down_individual,1))) / T;
    end
    if isnan(avg_raw_total_down_power_per_step) || avg_raw_total_down_power_per_step == 0; avg_raw_total_down_power_per_step = 80; end % 默认基础值
    demand_factor_down_low = 0.15; demand_factor_down_high = 0.4;
    P_grid_down_demand = avg_raw_total_down_power_per_step * (demand_factor_down_low + (demand_factor_down_high-demand_factor_down_low) * sin(linspace(0,2*pi,T)'*0.5 + pi/4).^2 );
    % results.P_grid_down_regulation_demand = P_grid_down_demand; % 可选
end

eps_val = 1e-6;

%% 2. 定义成本参数
C_ac_up = 0.05 * ones(num_ac_total, 1); 
C_ev_up = 0.04 * ones(num_ev_total, 1);
C_ac_down = 0.03 * ones(num_ac_total, 1);
C_ev_down = 0.02 * ones(num_ev_total, 1);

%% 3. 计算优化前的"原始"聚合能力和指标
AC_Up_Aggregated_Raw = zeros(T,1); EV_Up_Aggregated_Raw = zeros(T,1);
n_ac_raw_up_t = zeros(T,1); avg_P_ac_raw_up_t = zeros(T,1);
n_ev_raw_up_t = zeros(T,1); avg_P_ev_raw_up_t = zeros(T,1);

if num_ac_total > 0
    AC_Up_Aggregated_Raw = sum(P_ac_up_individual, 1)';
    for t_idx = 1:T
        n_ac_raw_up_t(t_idx) = num_ac_total;
        avg_P_ac_raw_up_t(t_idx) = sum(P_ac_up_individual(:, t_idx)) / (num_ac_total + eps_val);
    end
end
if num_ev_total > 0
    EV_Up_Aggregated_Raw = sum(P_ev_up_individual, 1)';
    for t_idx = 1:T
        n_ev_raw_up_t(t_idx) = num_ev_total;
        avg_P_ev_raw_up_t(t_idx) = sum(P_ev_up_individual(:, t_idx)) / (num_ev_total + eps_val);
    end
end
Total_Up_Aggregated_Raw = AC_Up_Aggregated_Raw + EV_Up_Aggregated_Raw;

if sum(n_ac_raw_up_t) == 0 && sum(n_ev_raw_up_t) == 0
    SDCI_up_raw = 0; rho_up_raw = 0;
else
    SDCI_up_raw = calculateSDCI(n_ac_raw_up_t, n_ev_raw_up_t, avg_P_ac_raw_up_t, avg_P_ev_raw_up_t);
    rho_up_raw = calculateSpearmanRho(n_ac_raw_up_t, avg_P_ac_raw_up_t, n_ev_raw_up_t, avg_P_ev_raw_up_t);
end

AC_Down_Aggregated_Raw = zeros(T,1); EV_Down_Aggregated_Raw = zeros(T,1);
n_ac_raw_down_t = zeros(T,1); avg_P_ac_raw_down_t = zeros(T,1);
n_ev_raw_down_t = zeros(T,1); avg_P_ev_raw_down_t = zeros(T,1);

if num_ac_total > 0
    AC_Down_Aggregated_Raw = sum(P_ac_down_individual, 1)';
    for t_idx = 1:T
        n_ac_raw_down_t(t_idx) = num_ac_total;
        avg_P_ac_raw_down_t(t_idx) = sum(P_ac_down_individual(:, t_idx)) / (num_ac_total + eps_val);
    end
end
if num_ev_total > 0
    EV_Down_Aggregated_Raw = sum(P_ev_down_individual, 1)';
    for t_idx = 1:T
        n_ev_raw_down_t(t_idx) = num_ev_total;
        avg_P_ev_raw_down_t(t_idx) = sum(P_ev_down_individual(:, t_idx)) / (num_ev_total + eps_val);
    end
end
Total_Down_Aggregated_Raw = AC_Down_Aggregated_Raw + EV_Down_Aggregated_Raw;

if sum(n_ac_raw_down_t) == 0 && sum(n_ev_raw_down_t) == 0
    SDCI_down_raw = 0; rho_down_raw = 0;
else
    SDCI_down_raw = calculateSDCI(n_ac_raw_down_t, n_ev_raw_down_t, avg_P_ac_raw_down_t, avg_P_ev_raw_down_t);
    rho_down_raw = calculateSpearmanRho(n_ac_raw_down_t, avg_P_ac_raw_down_t, n_ev_raw_down_t, avg_P_ev_raw_down_t);
end

%% 4. 设置遗传算法参数
pop_size_calc = max(20, (num_ac_total + num_ev_total));
gen_calc = max(10, T);
N_vars_for_ga_opts = (num_ac_total + num_ev_total) * T; % 用于确定合理的种群大小和代数

ga_opts = optimoptions('ga', ...
    'PopulationType', 'bitstring', ...
    'PopulationSize', min(200, max(50, round(N_vars_for_ga_opts / 1000 * 50))), ... % 调整与变量数的比例
    'MaxGenerations', min(500, max(100, round(N_vars_for_ga_opts / 1000 * 100))),... % 调整与变量数的比例
    'EliteCount', ceil(0.05 * min(200, max(50, round(N_vars_for_ga_opts / 1000 * 50)))), ...
    'CrossoverFraction', 0.8, ...
    'MutationFcn', {@mutationuniform, 0.01}, ...
    'Display', 'iter', ... % 减少迭代过程输出
    'UseParallel', true,...
    'PlotFcn', @gaplotbestf); 
% 如果 N_vars_for_ga_opts 为0 (例如T=0或无设备)，则上述计算种群大小等可能需要进一步保护，但调用GA前会检查

%% 5. 执行优化
% --- 上调优化 ---
U_ac_up_optimal = zeros(num_ac_total, T); 
U_ev_up_optimal = zeros(num_ev_total, T);
cost_up_optimal = 0; 
% exitflag_up = -1; % 如果不使用exitflag，可以不初始化

if (num_ac_total > 0 && ~isempty(P_ac_up_individual)) || (num_ev_total > 0 && ~isempty(P_ev_up_individual))
    C_ac_up_eff = C_ac_up; if num_ac_total == 0; C_ac_up_eff = []; end
    C_ev_up_eff = C_ev_up; if num_ev_total == 0; C_ev_up_eff = []; end
    
    [U_ac_up_optimal, U_ev_up_optimal, cost_up_optimal, ~, ~] = ... % 忽略 exitflag, output, population, scores
        solve_total_time_dispatch_ga(num_ac_total, num_ev_total, ...
                                     P_ac_up_individual, P_ev_up_individual, ...
                                     C_ac_up_eff, C_ev_up_eff, P_grid_up_demand, T, ...
                                     SDCI_up_raw, rho_up_raw, ga_opts);
end

% --- 下调优化 ---
U_ac_down_optimal = zeros(num_ac_total, T);
U_ev_down_optimal = zeros(num_ev_total, T);
cost_down_optimal = 0;
% exitflag_down = -1;

if (num_ac_total > 0 && ~isempty(P_ac_down_individual)) || (num_ev_total > 0 && ~isempty(P_ev_down_individual))
    C_ac_down_eff = C_ac_down; if num_ac_total == 0; C_ac_down_eff = []; end
    C_ev_down_eff = C_ev_down; if num_ev_total == 0; C_ev_down_eff = []; end

    [U_ac_down_optimal, U_ev_down_optimal, cost_down_optimal, ~, ~] = ...
        solve_total_time_dispatch_ga(num_ac_total, num_ev_total, ...
                                     P_ac_down_individual, P_ev_down_individual, ...
                                     C_ac_down_eff, C_ev_down_eff, P_grid_down_demand, T, ...
                                     SDCI_down_raw, rho_down_raw, ga_opts); 
end

%% 6. 计算优化后的聚合能力和指标
% --- 上调 Optimized ---
AC_Up_Optimal_Agg = zeros(T,1); EV_Up_Optimal_Agg = zeros(T,1);
n_ac_opt_up_t = zeros(T,1); avg_P_ac_opt_up_t = zeros(T,1);
n_ev_opt_up_t = zeros(T,1); avg_P_ev_opt_up_t = zeros(T,1);

if num_ac_total > 0 && size(U_ac_up_optimal,1) == num_ac_total && size(U_ac_up_optimal,2) == T
    for t_idx = 1:T
        AC_Up_Optimal_Agg(t_idx) = sum(P_ac_up_individual(:, t_idx) .* U_ac_up_optimal(:, t_idx));
        n_ac_opt_up_t(t_idx) = sum(U_ac_up_optimal(:, t_idx));
        avg_P_ac_opt_up_t(t_idx) = AC_Up_Optimal_Agg(t_idx) / (n_ac_opt_up_t(t_idx) + eps_val);
    end
end
if num_ev_total > 0 && size(U_ev_up_optimal,1) == num_ev_total && size(U_ev_up_optimal,2) == T
    for t_idx = 1:T
        EV_Up_Optimal_Agg(t_idx) = sum(P_ev_up_individual(:, t_idx) .* U_ev_up_optimal(:, t_idx));
        n_ev_opt_up_t(t_idx) = sum(U_ev_up_optimal(:, t_idx));
        avg_P_ev_opt_up_t(t_idx) = EV_Up_Optimal_Agg(t_idx) / (n_ev_opt_up_t(t_idx) + eps_val);
    end
end
Total_Up_Optimal_Agg = AC_Up_Optimal_Agg + EV_Up_Optimal_Agg;

if sum(n_ac_opt_up_t) == 0 && sum(n_ev_opt_up_t) == 0
    SDCI_up_opt = 0; rho_up_opt = 0;
else
    SDCI_up_opt = calculateSDCI(n_ac_opt_up_t, n_ev_opt_up_t, avg_P_ac_opt_up_t, avg_P_ev_opt_up_t);
    rho_up_opt = calculateSpearmanRho(n_ac_opt_up_t, avg_P_ac_opt_up_t, n_ev_opt_up_t, avg_P_ev_opt_up_t);
end

% --- 下调 Optimized ---
AC_Down_Optimal_Agg = zeros(T,1); EV_Down_Optimal_Agg = zeros(T,1);
n_ac_opt_down_t = zeros(T,1); avg_P_ac_opt_down_t = zeros(T,1);
n_ev_opt_down_t = zeros(T,1); avg_P_ev_opt_down_t = zeros(T,1);

if num_ac_total > 0 && size(U_ac_down_optimal,1) == num_ac_total && size(U_ac_down_optimal,2) == T
    for t_idx = 1:T
        AC_Down_Optimal_Agg(t_idx) = sum(P_ac_down_individual(:, t_idx) .* U_ac_down_optimal(:, t_idx));
        n_ac_opt_down_t(t_idx) = sum(U_ac_down_optimal(:, t_idx));
        avg_P_ac_opt_down_t(t_idx) = AC_Down_Optimal_Agg(t_idx) / (n_ac_opt_down_t(t_idx) + eps_val);
    end
end
if num_ev_total > 0 && size(U_ev_down_optimal,1) == num_ev_total && size(U_ev_down_optimal,2) == T
    for t_idx = 1:T
        EV_Down_Optimal_Agg(t_idx) = sum(P_ev_down_individual(:, t_idx) .* U_ev_down_optimal(:, t_idx));
        n_ev_opt_down_t(t_idx) = sum(U_ev_down_optimal(:, t_idx));
        avg_P_ev_opt_down_t(t_idx) = EV_Down_Optimal_Agg(t_idx) / (n_ev_opt_down_t(t_idx) + eps_val);
    end
end
Total_Down_Optimal_Agg = AC_Down_Optimal_Agg + EV_Down_Optimal_Agg;

if sum(n_ac_opt_down_t) == 0 && sum(n_ev_opt_down_t) == 0
    SDCI_down_opt = 0; rho_down_opt = 0;
else
    SDCI_down_opt = calculateSDCI(n_ac_opt_down_t, n_ev_opt_down_t, avg_P_ac_opt_down_t, avg_P_ev_opt_down_t);
    rho_down_opt = calculateSpearmanRho(n_ac_opt_down_t, avg_P_ac_opt_down_t, n_ev_opt_down_t, avg_P_ev_opt_down_t);
end

%% 7. 结果可视化
time_axis = (1:T)'; 

figure('Position', [100, 100, 1200, 800]);
sgtitle('VPP Dispatch Optimization Results (GA with SDCI/Rho Constraints)', 'FontSize', 16);

subplot(2,2,1);
plot(time_axis, Total_Up_Aggregated_Raw, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Raw Aggregated Up-Power'); hold on;
plot(time_axis, Total_Up_Optimal_Agg, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Optimized Up-Power (GA)');
plot(time_axis, P_grid_up_demand, 'k:', 'LineWidth', 1, 'DisplayName', 'Up-Regulation Demand');
hold off;
title('Up-Regulation Power Comparison'); xlabel('Time Step'); ylabel('Power (kW)');
legend('show', 'Location', 'best'); grid on;

subplot(2,2,2);
plot(time_axis, Total_Down_Aggregated_Raw, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Raw Aggregated Down-Power'); hold on;
plot(time_axis, Total_Down_Optimal_Agg, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Optimized Down-Power (GA)');
plot(time_axis, P_grid_down_demand, 'k:', 'LineWidth', 1, 'DisplayName', 'Down-Regulation Demand');
hold off;
title('Down-Regulation Power Comparison'); xlabel('Time Step'); ylabel('Power (kW)');
legend('show', 'Location', 'best'); grid on;

subplot(2,2,3);
indicator_names_up = {'SDCI⁺', 'Spearman ρ⁺'};
values_raw_up = [SDCI_up_raw; rho_up_raw];
values_opt_up = [SDCI_up_opt; rho_up_opt];
bar_data_up = [values_raw_up, values_opt_up];
b_up = bar(bar_data_up);
set(gca, 'XTickLabel', indicator_names_up);
ylabel('Indicator Value'); ylim([-1.1, 1.1]); 
legend([b_up(1) b_up(2)], {'Raw Aggregated', 'Optimized (GA)'}, 'Location', 'northoutside', 'Orientation','horizontal');
title('Up-Regulation Complementarity & Correlation'); grid on;
for k_bar = 1:size(bar_data_up,1) 
    for j_bar = 1:size(bar_data_up,2)
        text(b_up(j_bar).XData(k_bar) + b_up(j_bar).XOffset, bar_data_up(k_bar,j_bar), sprintf('%.3f', bar_data_up(k_bar,j_bar)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize',8, 'Color','k');
    end
end

subplot(2,2,4);
indicator_names_down = {'SDCI⁻', 'Spearman ρ⁻'};
values_raw_down = [SDCI_down_raw; rho_down_raw];
values_opt_down = [SDCI_down_opt; rho_down_opt];
bar_data_down = [values_raw_down, values_opt_down];
b_down = bar(bar_data_down);
set(gca, 'XTickLabel', indicator_names_down);
ylabel('Indicator Value'); ylim([-1.1, 1.1]);
legend([b_down(1) b_down(2)], {'Raw Aggregated', 'Optimized (GA)'}, 'Location', 'northoutside', 'Orientation','horizontal');
title('Down-Regulation Complementarity & Correlation'); grid on;
for k_bar = 1:size(bar_data_down,1)
    for j_bar = 1:size(bar_data_down,2)
        text(b_down(j_bar).XData(k_bar) + b_down(j_bar).XOffset, bar_data_down(k_bar,j_bar), sprintf('%.3f', bar_data_down(k_bar,j_bar)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize',8, 'Color','k');
    end
end

%% 8. 命令行输出汇总
disp(' ');
disp('=== 互补性与相关性指标对比 (GA Optimization) ===');
disp('【上调】');
fprintf('优化前 (Raw Aggregated): SDCI⁺ = %.4f, ρ⁺ = %.4f\n', SDCI_up_raw, rho_up_raw);
fprintf('优化后 (GA Optimized):  SDCI⁺ = %.4f, ρ⁺ = %.4f, 优化成本 = %.2f\n', SDCI_up_opt, rho_up_opt, cost_up_optimal);
disp(' ');
disp('【下调】');
fprintf('优化前 (Raw Aggregated): SDCI⁻ = %.4f, ρ⁻ = %.4f\n', SDCI_down_raw, rho_down_raw);
fprintf('优化后 (GA Optimized):  SDCI⁻ = %.4f, ρ⁻ = %.4f, 优化成本 = %.2f\n', SDCI_down_opt, rho_down_opt, cost_down_optimal);
disp(' ');