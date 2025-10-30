% eco_test_ga_ptdf.m
clc; close all; clear;

% --- 1. 加载数据与核心参数提取 (同 eco_test_ga.m) ---
load('results_single.mat');
P_ac_up_individual = results.AC_Up_Individual; P_ac_down_individual = abs(results.AC_Down_Individual);
P_ev_up_individual = results.EV_Up_Individual; P_ev_down_individual = abs(results.EV_Down_Individual);
num_ac_total = size(P_ac_up_individual, 1); num_ev_total = size(P_ev_up_individual, 1);
T = size(P_ac_up_individual, 2);
% (加载或生成 P_grid_up_demand, P_grid_down_demand 的逻辑同 eco_test_ga.m)
avg_raw_total_up_power = (sum(sum(P_ac_up_individual)) + sum(sum(P_ev_up_individual))) / (T+eps);
P_grid_up_demand = avg_raw_total_up_power * (0.2 + (0.5-0.2) * rand(T,1));
avg_raw_total_down_power = (sum(sum(P_ac_down_individual)) + sum(sum(P_ev_down_individual))) / (T+eps);
P_grid_down_demand = avg_raw_total_down_power * (0.15 + (0.4-0.15) * rand(T,1));
eps_val = 1e-6;

% --- 2. 定义成本参数 (同 eco_test_ga.m) ---
C_ac_up = 0.05 * ones(num_ac_total, 1); C_ev_up = 0.04 * ones(num_ev_total, 1);
C_ac_down = 0.03 * ones(num_ac_total, 1); C_ev_down = 0.02 * ones(num_ev_total, 1);

% --- 3. 计算原始指标 (同 eco_test_ga.m) ---
% (计算 SDCI_up_raw, rho_up_raw, SDCI_down_raw, rho_down_raw 的逻辑同 eco_test_ga.m)
AC_Up_Aggregated_Raw = sum(P_ac_up_individual, 1)'; EV_Up_Aggregated_Raw = sum(P_ev_up_individual, 1)';
n_ac_raw_up_t = ones(T,1)*num_ac_total; avg_P_ac_raw_up_t = AC_Up_Aggregated_Raw./(n_ac_raw_up_t+eps);
n_ev_raw_up_t = ones(T,1)*num_ev_total; avg_P_ev_raw_up_t = EV_Up_Aggregated_Raw./(n_ev_raw_up_t+eps);
SDCI_up_raw = calculateSDCI(n_ac_raw_up_t, n_ev_raw_up_t, avg_P_ac_raw_up_t, avg_P_ev_raw_up_t);
rho_up_raw = calculateSpearmanRho(n_ac_raw_up_t, avg_P_ac_raw_up_t, n_ev_raw_up_t, avg_P_ev_raw_up_t);
% (下调同理)
AC_Down_Aggregated_Raw=sum(P_ac_down_individual, 1)'; EV_Down_Aggregated_Raw=sum(P_ev_down_individual, 1)';
n_ac_raw_down_t=ones(T,1)*num_ac_total; avg_P_ac_raw_down_t=AC_Down_Aggregated_Raw./(n_ac_raw_down_t+eps);
n_ev_raw_down_t=ones(T,1)*num_ev_total; avg_P_ev_raw_down_t=EV_Down_Aggregated_Raw./(n_ev_raw_down_t+eps);
SDCI_down_raw = calculateSDCI(n_ac_raw_down_t, n_ev_raw_down_t, avg_P_ac_raw_down_t, avg_P_ev_raw_down_t);
rho_down_raw = calculateSpearmanRho(n_ac_raw_down_t, avg_P_ac_raw_down_t, n_ev_raw_down_t, avg_P_ev_raw_down_t);


% --- *** 新增：4. 定义网络拓扑参数 (示例) *** ---
N_bus = 5;  % 假设有5个节点
N_line = 4; % 假设有4条线路

% 随机生成设备位置 (确保 Location 向量长度与设备数一致)
if num_ac_total > 0; Location_AC = randi([1, N_bus], num_ac_total, 1); else; Location_AC = []; end
if num_ev_total > 0; Location_EV = randi([1, N_bus], num_ev_total, 1); else; Location_EV = []; end

% 生成示例 PTDF 矩阵 (实际应用中需要根据电网拓扑计算)
PTDF_matrix = rand(N_line, N_bus) * 0.2 - 0.1; % 示例值

% 生成示例线路基础潮流 (例如，随时间波动)
P_Line_Base = rand(N_line, T) * 50 - 25; % 示例值

% 定义线路容量
P_Line_Max = rand(N_line, 1) * 50 + 50; % 示例值 (50-100 kW)

% --- 5. 设置遗传算法参数 (同 eco_test_ga.m) ---
ga_opts = optimoptions('ga', 'PopulationType', 'bitstring', 'Display', 'iter', 'UseParallel', true, 'PlotFcn', @gaplotbestf);
% (可以根据 eco_test_ga.m 中的逻辑动态调整 PopulationSize, MaxGenerations 等)
N_vars_for_ga_opts = (num_ac_total + num_ev_total) * T;
ga_opts.PopulationSize = min(200, max(50, round(N_vars_for_ga_opts / 1000 * 50)));
ga_opts.MaxGenerations = min(500, max(100, round(N_vars_for_ga_opts / 1000 * 100)));
ga_opts.EliteCount = ceil(0.05 * ga_opts.PopulationSize);

% --- 6. 执行优化 (调用新函数) ---
U_ac_up_optimal = zeros(num_ac_total, T); U_ev_up_optimal = zeros(num_ev_total, T); cost_up_optimal = 0;
if num_ac_total > 0 || num_ev_total > 0
    fprintf('\n--- 开始上调优化 (带PTDF约束) ---\n');
    [U_ac_up_optimal, U_ev_up_optimal, cost_up_optimal, ~, ~, ~, ~] = ...
        solve_total_time_dispatch_ga_ptdf(num_ac_total, num_ev_total, ...
                                     P_ac_up_individual, P_ev_up_individual, ...
                                     C_ac_up, C_ev_up, P_grid_up_demand, T, ...
                                     SDCI_up_raw, rho_up_raw, ga_opts, ...
                                     Location_AC, Location_EV, PTDF_matrix, ... % 传递新参数
                                     P_Line_Base, P_Line_Max, N_bus, N_line);   % 传递新参数
end

U_ac_down_optimal = zeros(num_ac_total, T); U_ev_down_optimal = zeros(num_ev_total, T); cost_down_optimal = 0;
if num_ac_total > 0 || num_ev_total > 0
     fprintf('\n--- 开始下调优化 (带PTDF约束) ---\n');
     [U_ac_down_optimal, U_ev_down_optimal, cost_down_optimal, ~, ~, ~, ~] = ...
        solve_total_time_dispatch_ga_ptdf(num_ac_total, num_ev_total, ...
                                     P_ac_down_individual, P_ev_down_individual, ...
                                     C_ac_down, C_ev_down, P_grid_down_demand, T, ...
                                     SDCI_down_raw, rho_down_raw, ga_opts, ...
                                     Location_AC, Location_EV, PTDF_matrix, ... % 传递新参数
                                     P_Line_Base, P_Line_Max, N_bus, N_line);   % 传递新参数
end

% --- 7. 计算优化后指标 (同 eco_test_ga.m) ---
% (计算 SDCI_up_opt, rho_up_opt, SDCI_down_opt, rho_down_opt 的逻辑同 eco_test_ga.m)
AC_Up_Optimal_Agg = zeros(T,1); EV_Up_Optimal_Agg = zeros(T,1);
n_ac_opt_up_t = zeros(T,1); avg_P_ac_opt_up_t = zeros(T,1);
n_ev_opt_up_t = zeros(T,1); avg_P_ev_opt_up_t = zeros(T,1);
if num_ac_total > 0; for t=1:T; AC_Up_Optimal_Agg(t)=sum(P_ac_up_individual(:,t).*U_ac_up_optimal(:,t)); n_ac_opt_up_t(t)=sum(U_ac_up_optimal(:,t)); avg_P_ac_opt_up_t(t)=AC_Up_Optimal_Agg(t)/(n_ac_opt_up_t(t)+eps); end; end
if num_ev_total > 0; for t=1:T; EV_Up_Optimal_Agg(t)=sum(P_ev_up_individual(:,t).*U_ev_up_optimal(:,t)); n_ev_opt_up_t(t)=sum(U_ev_up_optimal(:,t)); avg_P_ev_opt_up_t(t)=EV_Up_Optimal_Agg(t)/(n_ev_opt_up_t(t)+eps); end; end
SDCI_up_opt = calculateSDCI(n_ac_opt_up_t, n_ev_opt_up_t, avg_P_ac_opt_up_t, avg_P_ev_opt_up_t);
rho_up_opt = calculateSpearmanRho(n_ac_opt_up_t, avg_P_ac_opt_up_t, n_ev_opt_up_t, avg_P_ev_opt_up_t);
% (下调同理)
AC_Down_Optimal_Agg = zeros(T,1); EV_Down_Optimal_Agg = zeros(T,1);
n_ac_opt_down_t = zeros(T,1); avg_P_ac_opt_down_t = zeros(T,1);
n_ev_opt_down_t = zeros(T,1); avg_P_ev_opt_down_t = zeros(T,1);
if num_ac_total > 0; for t=1:T; AC_Down_Optimal_Agg(t)=sum(P_ac_down_individual(:,t).*U_ac_down_optimal(:,t)); n_ac_opt_down_t(t)=sum(U_ac_down_optimal(:,t)); avg_P_ac_opt_down_t(t)=AC_Down_Optimal_Agg(t)/(n_ac_opt_down_t(t)+eps); end; end
if num_ev_total > 0; for t=1:T; EV_Down_Optimal_Agg(t)=sum(P_ev_down_individual(:,t).*U_ev_down_optimal(:,t)); n_ev_opt_down_t(t)=sum(U_ev_down_optimal(:,t)); avg_P_ev_opt_down_t(t)=EV_Down_Optimal_Agg(t)/(n_ev_opt_down_t(t)+eps); end; end
SDCI_down_opt = calculateSDCI(n_ac_opt_down_t, n_ev_opt_down_t, avg_P_ac_opt_down_t, avg_P_ev_opt_down_t);
rho_down_opt = calculateSpearmanRho(n_ac_opt_down_t, avg_P_ac_opt_down_t, n_ev_opt_down_t, avg_P_ev_opt_down_t);


% --- 8. 结果可视化与输出 (同 eco_test_ga.m) ---
time_axis = (1:T)'; Total_Up_Optimal_Agg = AC_Up_Optimal_Agg + EV_Up_Optimal_Agg; Total_Down_Optimal_Agg = AC_Down_Optimal_Agg + EV_Down_Optimal_Agg;
Total_Up_Aggregated_Raw = AC_Up_Aggregated_Raw + EV_Up_Aggregated_Raw; Total_Down_Aggregated_Raw = AC_Down_Aggregated_Raw + EV_Down_Aggregated_Raw;
figure('Position', [100, 100, 1200, 800]); sgtitle('VPP Dispatch Optimization (GA with PTDF)', 'FontSize', 16);
subplot(2,2,1); plot(time_axis, Total_Up_Aggregated_Raw, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Raw Aggregated'); hold on; plot(time_axis, Total_Up_Optimal_Agg, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Optimized (GA+PTDF)'); plot(time_axis, P_grid_up_demand, 'k:', 'LineWidth', 1, 'DisplayName', 'Demand'); hold off; title('Up-Regulation Power'); legend('show'); grid on;
subplot(2,2,2); plot(time_axis, Total_Down_Aggregated_Raw, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Raw Aggregated'); hold on; plot(time_axis, Total_Down_Optimal_Agg, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Optimized (GA+PTDF)'); plot(time_axis, P_grid_down_demand, 'k:', 'LineWidth', 1, 'DisplayName', 'Demand'); hold off; title('Down-Regulation Power'); legend('show'); grid on;
subplot(2,2,3); bar_data_up = [SDCI_up_raw, SDCI_up_opt; rho_up_raw, rho_up_opt]; b_up = bar(bar_data_up); set(gca, 'XTickLabel', {'SDCI⁺', 'Spearman ρ⁺'}); legend([b_up(1) b_up(2)], {'Raw', 'Optimized (GA+PTDF)'}); title('Up-Regulation Metrics'); grid on; ylim([-1.1 1.1]); % add_bar_labels(b_up, bar_data_up);
subplot(2,2,4); bar_data_down = [SDCI_down_raw, SDCI_down_opt; rho_down_raw, rho_down_opt]; b_down = bar(bar_data_down); set(gca, 'XTickLabel', {'SDCI⁻', 'Spearman ρ⁻'}); legend([b_down(1) b_down(2)], {'Raw', 'Optimized (GA+PTDF)'}); title('Down-Regulation Metrics'); grid on; ylim([-1.1 1.1]); % add_bar_labels(b_down, bar_data_down);

disp(' '); disp('=== 指标对比 (GA Optimization with PTDF) ===');
disp('【上调】'); fprintf('优化前: SDCI⁺ = %.4f, ρ⁺ = %.4f\n', SDCI_up_raw, rho_up_raw); fprintf('优化后: SDCI⁺ = %.4f, ρ⁺ = %.4f, 成本 = %.2f\n', SDCI_up_opt, rho_up_opt, cost_up_optimal);
disp(' '); disp('【下调】'); fprintf('优化前: SDCI⁻ = %.4f, ρ⁻ = %.4f\n', SDCI_down_raw, rho_down_raw); fprintf('优化后: SDCI⁻ = %.4f, ρ⁻ = %.4f, 成本 = %.2f\n', SDCI_down_opt, rho_down_opt, cost_down_optimal); disp(' ');