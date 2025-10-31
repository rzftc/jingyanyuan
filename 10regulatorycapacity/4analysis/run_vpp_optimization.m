% run_vpp_optimization.m
clear; clc; close all;

%% 1. 加载和准备数据 (已按您的要求修改)
fprintf('正在加载仿真数据...\n');

% --- 指定要加载的分块文件名 ---
chunk_file_to_load = 'chunk_results/results_chunk_1.mat'; % 您可以修改这个文件名

try
    data_struct = load(chunk_file_to_load);
    % 数据位于加载后的 .results 字段中
    if isfield(data_struct, 'results')
        sim_results = data_struct.results;
    else
        error('在 %s 文件中未找到 "results" 结构体。', chunk_file_to_load);
    end
catch ME
    error('加载数据文件 %s 失败。请确保文件路径正确。\n错误信息: %s', chunk_file_to_load, ME.message);
end

% --- 从 results 结构体中提取数据 ---
% 假设是上调场景, 使用 AC_Up_Individual 和 EV_Up_Individual
P_ac_potential = [];
P_ev_potential = [];

if isfield(sim_results, 'AC_Up_Individual')
    P_ac_potential = sim_results.AC_Up_Individual;
    fprintf('成功提取空调(AC)数据。\n');
else
    warning('在加载的文件中未找到 AC_Up_Individual 字段。');
end

if isfield(sim_results, 'EV_Up_Individual')
    P_ev_potential = sim_results.EV_Up_Individual;
    fprintf('成功提取电动汽车(EV)数据。\n');
else
    warning('在加载的文件中未找到 EV_Up_Individual 字段。');
end

if isempty(P_ac_potential) && isempty(P_ev_potential)
    error('未能从文件中加载任何设备潜力数据，优化无法进行。');
end

% 获取维度信息
[num_ac, T_ac] = size(P_ac_potential);
[num_ev, T_ev] = size(P_ev_potential);
T = max(T_ac, T_ev); % 取两者中最大的时间步数

% 填充空数据以保持维度一致
if isempty(P_ac_potential)
    P_ac_potential = zeros(0, T);
    num_ac = 0;
end
if isempty(P_ev_potential)
    P_ev_potential = zeros(0, T);
    num_ev = 0;
end

fprintf('数据加载完成: %d 台空调, %d 台电动汽车, %d 个时间步。\n', num_ac, num_ev, T);

%% 2. 构建优化输入参数
% --- a. 设备参数 ---
ac_params.num = num_ac;
ac_params.P_potential = P_ac_potential;
ac_params.Cost = 0.05 * ones(num_ac, 1); % 示例成本: 0.05 元/kW
ac_params.Location = randi([1, 10], num_ac, 1); % 随机分配到10个节点中的一个

ev_params.num = num_ev;
ev_params.P_potential = P_ev_potential;
ev_params.Cost = 0.04 * ones(num_ev, 1); % 示例成本: 0.04 元/kW
ev_params.Location = randi([1, 10], num_ev, 1);

% --- b. 电网与需求参数 ---
grid_params.T = T;
% 示例需求：聚合潜力的20%到50%
total_potential = sum(P_ac_potential,1) + sum(P_ev_potential,1);
grid_params.P_req = (0.2 + 0.3*rand(T,1)) .* total_potential';

% 计算基准指标 (全员参与时)
n_ac_orig = num_ac * ones(T,1);
n_ev_orig = num_ev * ones(T,1);
avg_p_ac_orig = (sum(P_ac_potential,1)./(num_ac+eps))';
avg_p_ev_orig = (sum(P_ev_potential,1)./(num_ev+eps))';
grid_params.SDCI_orig = calculateSDCI(n_ac_orig, n_ev_orig, avg_p_ac_orig, avg_p_ev_orig);
grid_params.rho_orig = calculateSpearmanRho(n_ac_orig, avg_p_ac_orig, n_ev_orig, avg_p_ev_orig);
fprintf('基准指标: SDCI = %.4f, Rho = %.4f\n', grid_params.SDCI_orig, grid_params.rho_orig);


% --- c. 网络拓扑参数 (示例) ---
grid_params.N_bus = 10;
grid_params.N_line = 8;
grid_params.PTDF = rand(grid_params.N_line, grid_params.N_bus) * 0.2 - 0.1; % 示例 PTDF
grid_params.P_line_base = (rand(grid_params.N_line, T) - 0.5) * 50; % 示例基础潮流
grid_params.P_line_max = 80 * ones(grid_params.N_line, 1); % 线路潮流上限 80 kW

%% 3. 配置并运行遗传算法
% 针对大规模问题，建议增大种群和代数，并启用并行计算
ga_options = optimoptions('ga', ...
    'PopulationType', 'doubleVector',... % ***【修正】***: 必须为 'doubleVector' 才能配合 intcon 使用
    'PopulationSize', 200, ...         
    'MaxGenerations', 500, ...        
    'EliteCount', 10, ...
    'CrossoverFraction', 0.8, ...
    'MutationFcn', {@mutationadaptfeasible},... % ***【修正】***: 使用 'doubleVector' 对应的变异函数
    'Display', 'iter', ...
    'PlotFcn', @gaplotbestf, ...
    'UseParallel', true);             % 启用并行计算

% 调用核心求解器 (不变)
[U_ac_opt, U_ev_opt, cost_opt, exitflag, ga_output] = solve_vpp_optimization_main(ac_params, ev_params, grid_params, ga_options);

% (上一版中的备用方案已移除，因为此配置是GA处理整数约束的标准方式)

disp(ga_output.message);

%% 4. 结果分析与可视化
fprintf('优化完成！最低总成本: %.2f 元\n', cost_opt);

% --- a. 功率满足情况 ---
P_ac_dispatch_opt = sum(U_ac_opt .* P_ac_potential, 1);
P_ev_dispatch_opt = sum(U_ev_opt .* P_ev_potential, 1);
P_total_dispatch_opt = P_ac_dispatch_opt + P_ev_dispatch_opt;

figure('Position', [100, 100, 1200, 600]);
subplot(2,1,1)
hold on;
plot(1:T, grid_params.P_req, 'k--', 'LineWidth', 1.5, 'DisplayName', '需求功率');
plot(1:T, P_total_dispatch_opt, 'r-', 'LineWidth', 1.5, 'DisplayName', '调度功率');
area(1:T, P_ac_dispatch_opt, 'FaceColor', 'b', 'FaceAlpha', 0.3, 'DisplayName', '空调出力');
area(1:T, P_ev_dispatch_opt, 'FaceColor', 'g', 'FaceAlpha', 0.3, 'DisplayName', 'EV出力');
hold off;
title('优化调度功率 vs. 电网需求');
xlabel('时间步'); ylabel('功率 (kW)');
legend('Location', 'best'); grid on;

% --- b. 指标改善情况 ---
n_ac_opt_t = sum(U_ac_opt, 1)';
n_ev_opt_t = sum(U_ev_opt, 1)';
avg_p_ac_opt_t = (P_ac_dispatch_opt ./ (n_ac_opt_t' + eps))';
avg_p_ev_opt_t = (P_ev_dispatch_opt ./ (n_ev_opt_t' + eps))';

sdci_final = calculateSDCI(n_ac_opt_t, n_ev_opt_t, avg_p_ac_opt_t, avg_p_ev_opt_t);
rho_final = calculateSpearmanRho(n_ac_opt_t, avg_p_ac_opt_t, n_ev_opt_t, avg_p_ev_opt_t);

subplot(2,2,3);
bar_data = [grid_params.SDCI_orig, sdci_final; grid_params.rho_orig, rho_final];
bar(bar_data);
set(gca, 'XTickLabel', {'SDCI', 'Spearman Rho'});
ylabel('指标值'); ylim([-1.1, 1.1]);
title('互补性与相关性指标对比');
legend('优化前', '优化后'); grid on;

% --- c. 线路潮流检查 (抽样) ---
Delta_Pinj_opt = zeros(grid_params.N_bus, T);
for k = 1:grid_params.N_bus
    ac_at_bus_k = (ac_params.Location == k);
    ev_at_bus_k = (ev_params.Location == k);
    if any(ac_at_bus_k), Delta_Pinj_opt(k,:) = Delta_Pinj_opt(k,:) + sum(U_ac_opt(ac_at_bus_k,:) .* ac_params.P_potential(ac_at_bus_k,:), 1); end
    if any(ev_at_bus_k), Delta_Pinj_opt(k,:) = Delta_Pinj_opt(k,:) + sum(U_ev_opt(ev_at_bus_k,:) .* ev_params.P_potential(ev_at_bus_k,:), 1); end
end
Delta_P_line_opt = grid_params.PTDF * Delta_Pinj_opt;
P_line_final_opt = grid_params.P_line_base + Delta_P_line_opt;

max_loading_pct = max(abs(P_line_final_opt) ./ grid_params.P_line_max, [], 'all') * 100;
subplot(2,2,4)
plot(1:T, P_line_final_opt');
hold on;
yline(grid_params.P_line_max(1), 'r--', 'LineWidth', 1.5);
yline(-grid_params.P_line_max(1), 'r--', 'LineWidth', 1.5);
hold off;
title(sprintf('线路潮流 (最高负载率: %.1f%%)', max_loading_pct));
xlabel('时间步'); ylabel('潮流 (kW)');
grid on;