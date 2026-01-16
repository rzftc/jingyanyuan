%% 测试脚本 (已修改为调用 solve_total_time_dispatch_ga_pv - 含光伏PV版)
clc; close all; clear;

% 假设 solve_total_time_dispatch_ga_pv.m, objective_function_ga_pv.m, 
% calculateSDCI.m, calculateSpearmanRho.m 在MATLAB路径中

%% 1. 加载数据与核心参数提取
load('results_single.mat'); % 假设此mat文件包含所有必要的个体数据和可能的电网需求

% 直接从 results 结构体中提取个体数据
P_ac_up_individual = results.AC_Up_Individual;
P_ac_down_individual = abs(results.AC_Down_Individual); % 确保为正
P_ev_up_individual = results.EV_Up_Individual;
P_ev_down_individual = abs(results.EV_Down_Individual); % 确保为正

num_ac_total = size(P_ac_up_individual, 1);
num_ev_total = size(P_ev_up_individual, 1);
T = size(P_ac_up_individual, 2); 

% --- 新增：生成光伏 (PV) 调节能力数据 ---
% 依据：根据AC和EV的调节能力进行生成 (取均值后附加随机波动)
num_pv_total = round((num_ac_total + num_ev_total) / 2); % 假设PV数量为两者的均值
if num_pv_total == 0; num_pv_total = 10; end

avg_ac_up = mean(P_ac_up_individual(:));
avg_ev_up = mean(P_ev_up_individual(:));
base_pv_up = (avg_ac_up + avg_ev_up) / 2;
% 生成PV上调能力 (假设有一定的随机性)
P_pv_up_individual = base_pv_up * (0.8 + 0.4 * rand(num_pv_total, T));

avg_ac_down = mean(P_ac_down_individual(:));
avg_ev_down = mean(P_ev_down_individual(:));
base_pv_down = (avg_ac_down + avg_ev_down) / 2;
% 生成PV下调能力
P_pv_down_individual = base_pv_down * (0.8 + 0.4 * rand(num_pv_total, T));
% ---------------------------------------

% 定义/加载电网调节需求 (需包含PV的潜在贡献考虑)
if isfield(results, 'P_grid_up_regulation_demand') && ~isempty(results.P_grid_up_regulation_demand) && length(results.P_grid_up_regulation_demand) == T
    P_grid_up_demand = results.P_grid_up_regulation_demand(:);
else
    avg_raw_total_up_power_per_step = 0;
    if T > 0
        % 需求计算加入PV的基础能力
        avg_raw_total_up_power_per_step = (sum(sum(P_ac_up_individual,1)) + sum(sum(P_ev_up_individual,1)) + sum(sum(P_pv_up_individual,1))) / T;
    end
    if isnan(avg_raw_total_up_power_per_step) || avg_raw_total_up_power_per_step == 0; avg_raw_total_up_power_per_step = 100; end 
    demand_factor_up_low = 0.2; demand_factor_up_high = 0.5;
    P_grid_up_demand = avg_raw_total_up_power_per_step * (demand_factor_up_low + (demand_factor_up_high-demand_factor_up_low) * rand(T,1));
end

if isfield(results, 'P_grid_down_regulation_demand') && ~isempty(results.P_grid_down_regulation_demand) && length(results.P_grid_down_regulation_demand) == T
    P_grid_down_demand = results.P_grid_down_regulation_demand(:);
else
    avg_raw_total_down_power_per_step = 0;
    if T > 0
        % 需求计算加入PV的基础能力
        avg_raw_total_down_power_per_step = (sum(sum(P_ac_down_individual,1)) + sum(sum(P_ev_down_individual,1)) + sum(sum(P_pv_down_individual,1))) / T;
    end
    if isnan(avg_raw_total_down_power_per_step) || avg_raw_total_down_power_per_step == 0; avg_raw_total_down_power_per_step = 80; end 
    demand_factor_down_low = 0.15; demand_factor_down_high = 0.4;
    P_grid_down_demand = avg_raw_total_down_power_per_step * (demand_factor_down_low + (demand_factor_down_high-demand_factor_down_low) * sin(linspace(0,2*pi,T)'*0.5 + pi/4).^2 );
end

eps_val = 1e-6;

%% 2. 定义成本参数
C_ac_up = 0.05 * ones(num_ac_total, 1); 
C_ev_up = 0.04 * ones(num_ev_total, 1);
C_pv_up = 0.03 * ones(num_pv_total, 1); % 新增PV成本 (假设比AC/EV稍便宜)

C_ac_down = 0.03 * ones(num_ac_total, 1);
C_ev_down = 0.02 * ones(num_ev_total, 1);
C_pv_down = 0.015 * ones(num_pv_total, 1); % 新增PV成本

%% 3. 计算优化前的"原始"聚合能力和指标
AC_Up_Aggregated_Raw = zeros(T,1); EV_Up_Aggregated_Raw = zeros(T,1); PV_Up_Aggregated_Raw = zeros(T,1);
n_ac_raw_up_t = zeros(T,1); avg_P_ac_raw_up_t = zeros(T,1);
n_ev_raw_up_t = zeros(T,1); avg_P_ev_raw_up_t = zeros(T,1);
% PV原始统计暂不纳入SDCI计算（SDCI原定义为AC与EV互补），但计入总功率

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
if num_pv_total > 0
    PV_Up_Aggregated_Raw = sum(P_pv_up_individual, 1)';
end

Total_Up_Aggregated_Raw = AC_Up_Aggregated_Raw + EV_Up_Aggregated_Raw + PV_Up_Aggregated_Raw; % 加入PV

% 指标计算保持原样 (AC vs EV)
if sum(n_ac_raw_up_t) == 0 && sum(n_ev_raw_up_t) == 0
    SDCI_up_raw = 0; rho_up_raw = 0;
else
    SDCI_up_raw = calculateSDCI(n_ac_raw_up_t, n_ev_raw_up_t, avg_P_ac_raw_up_t, avg_P_ev_raw_up_t);
    rho_up_raw = calculateSpearmanRho(n_ac_raw_up_t, avg_P_ac_raw_up_t, n_ev_raw_up_t, avg_P_ev_raw_up_t);
end

AC_Down_Aggregated_Raw = zeros(T,1); EV_Down_Aggregated_Raw = zeros(T,1); PV_Down_Aggregated_Raw = zeros(T,1);
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
if num_pv_total > 0
    PV_Down_Aggregated_Raw = sum(P_pv_down_individual, 1)';
end

Total_Down_Aggregated_Raw = AC_Down_Aggregated_Raw + EV_Down_Aggregated_Raw + PV_Down_Aggregated_Raw; % 加入PV

if sum(n_ac_raw_down_t) == 0 && sum(n_ev_raw_down_t) == 0
    SDCI_down_raw = 0; rho_down_raw = 0;
else
    SDCI_down_raw = calculateSDCI(n_ac_raw_down_t, n_ev_raw_down_t, avg_P_ac_raw_down_t, avg_P_ev_raw_down_t);
    rho_down_raw = calculateSpearmanRho(n_ac_raw_down_t, avg_P_ac_raw_down_t, n_ev_raw_down_t, avg_P_ev_raw_down_t);
end

%% 4. 设置遗传算法参数
pop_size_calc = max(20, (num_ac_total + num_ev_total + num_pv_total));
gen_calc = max(10, T);
% 变量总数增加 PV 部分
N_vars_for_ga_opts = (num_ac_total + num_ev_total + num_pv_total) * T; 

ga_opts = optimoptions('ga', ...
    'PopulationType', 'bitstring', ...
    'PopulationSize', min(50, max(50, round(N_vars_for_ga_opts / 1000 * 50))), ... 
    'MaxGenerations', min(100, max(100, round(N_vars_for_ga_opts / 1000 * 100))),... 
    'EliteCount', ceil(0.05 * min(200, max(50, round(N_vars_for_ga_opts / 1000 * 50)))), ...
    'CrossoverFraction', 0.8, ...
    'MutationFcn', {@mutationuniform, 0.01}, ...
    'Display', 'iter', ... 
    'UseParallel', true,...
    'PlotFcn', @gaplotbestf); 

%% 5. 执行优化
% --- 上调优化 ---
U_ac_up_optimal = zeros(num_ac_total, T); 
U_ev_up_optimal = zeros(num_ev_total, T);
U_pv_up_optimal = zeros(num_pv_total, T); % 新增PV结果容器
cost_up_optimal = 0; 

if (num_ac_total > 0 && ~isempty(P_ac_up_individual)) || (num_ev_total > 0 && ~isempty(P_ev_up_individual)) || (num_pv_total > 0)
    C_ac_up_eff = C_ac_up; if num_ac_total == 0; C_ac_up_eff = []; end
    C_ev_up_eff = C_ev_up; if num_ev_total == 0; C_ev_up_eff = []; end
    C_pv_up_eff = C_pv_up; if num_pv_total == 0; C_pv_up_eff = []; end
    
    % 调用新命名的函数
    [U_ac_up_optimal, U_ev_up_optimal, U_pv_up_optimal, cost_up_optimal, ~, ~] = ... 
        solve_total_time_dispatch_ga_pv(num_ac_total, num_ev_total, num_pv_total, ...
                                     P_ac_up_individual, P_ev_up_individual, P_pv_up_individual, ...
                                     C_ac_up_eff, C_ev_up_eff, C_pv_up_eff, ...
                                     P_grid_up_demand, T, ...
                                     SDCI_up_raw, rho_up_raw, ga_opts);
end

% --- 下调优化 ---
U_ac_down_optimal = zeros(num_ac_total, T);
U_ev_down_optimal = zeros(num_ev_total, T);
U_pv_down_optimal = zeros(num_pv_total, T); % 新增PV结果容器
cost_down_optimal = 0;

if (num_ac_total > 0 && ~isempty(P_ac_down_individual)) || (num_ev_total > 0 && ~isempty(P_ev_down_individual)) || (num_pv_total > 0)
    C_ac_down_eff = C_ac_down; if num_ac_total == 0; C_ac_down_eff = []; end
    C_ev_down_eff = C_ev_down; if num_ev_total == 0; C_ev_down_eff = []; end
    C_pv_down_eff = C_pv_down; if num_pv_total == 0; C_pv_down_eff = []; end

    % 调用新命名的函数
    [U_ac_down_optimal, U_ev_down_optimal, U_pv_down_optimal, cost_down_optimal, ~, ~] = ...
        solve_total_time_dispatch_ga_pv(num_ac_total, num_ev_total, num_pv_total, ...
                                     P_ac_down_individual, P_ev_down_individual, P_pv_down_individual, ...
                                     C_ac_down_eff, C_ev_down_eff, C_pv_down_eff, ...
                                     P_grid_down_demand, T, ...
                                     SDCI_down_raw, rho_down_raw, ga_opts); 
end

%% 6. 计算优化后的聚合能力和指标
% --- 上调 Optimized ---
AC_Up_Optimal_Agg = zeros(T,1); EV_Up_Optimal_Agg = zeros(T,1); PV_Up_Optimal_Agg = zeros(T,1);
n_ac_opt_up_t = zeros(T,1); avg_P_ac_opt_up_t = zeros(T,1);
n_ev_opt_up_t = zeros(T,1); avg_P_ev_opt_up_t = zeros(T,1);
% PV Stats
n_pv_opt_up_t = zeros(T,1);

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
% PV Optimized Aggregation
if num_pv_total > 0 && size(U_pv_up_optimal,1) == num_pv_total && size(U_pv_up_optimal,2) == T
    for t_idx = 1:T
        PV_Up_Optimal_Agg(t_idx) = sum(P_pv_up_individual(:, t_idx) .* U_pv_up_optimal(:, t_idx));
        n_pv_opt_up_t(t_idx) = sum(U_pv_up_optimal(:, t_idx));
    end
end

Total_Up_Optimal_Agg = AC_Up_Optimal_Agg + EV_Up_Optimal_Agg + PV_Up_Optimal_Agg;

if sum(n_ac_opt_up_t) == 0 && sum(n_ev_opt_up_t) == 0
    SDCI_up_opt = 0; rho_up_opt = 0;
else
    SDCI_up_opt = calculateSDCI(n_ac_opt_up_t, n_ev_opt_up_t, avg_P_ac_opt_up_t, avg_P_ev_opt_up_t);
    rho_up_opt = calculateSpearmanRho(n_ac_opt_up_t, avg_P_ac_opt_up_t, n_ev_opt_up_t, avg_P_ev_opt_up_t);
end

% --- 下调 Optimized ---
AC_Down_Optimal_Agg = zeros(T,1); EV_Down_Optimal_Agg = zeros(T,1); PV_Down_Optimal_Agg = zeros(T,1);
n_ac_opt_down_t = zeros(T,1); avg_P_ac_opt_down_t = zeros(T,1);
n_ev_opt_down_t = zeros(T,1); avg_P_ev_opt_down_t = zeros(T,1);
n_pv_opt_down_t = zeros(T,1);

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
% PV Optimized Aggregation
if num_pv_total > 0 && size(U_pv_down_optimal,1) == num_pv_total && size(U_pv_down_optimal,2) == T
    for t_idx = 1:T
        PV_Down_Optimal_Agg(t_idx) = sum(P_pv_down_individual(:, t_idx) .* U_pv_down_optimal(:, t_idx));
        n_pv_opt_down_t(t_idx) = sum(U_pv_down_optimal(:, t_idx));
    end
end

Total_Down_Optimal_Agg = AC_Down_Optimal_Agg + EV_Down_Optimal_Agg + PV_Down_Optimal_Agg;

if sum(n_ac_opt_down_t) == 0 && sum(n_ev_opt_down_t) == 0
    SDCI_down_opt = 0; rho_down_opt = 0;
else
    SDCI_down_opt = calculateSDCI(n_ac_opt_down_t, n_ev_opt_down_t, avg_P_ac_opt_down_t, avg_P_ev_opt_down_t);
    rho_down_opt = calculateSpearmanRho(n_ac_opt_down_t, avg_P_ac_opt_down_t, n_ev_opt_down_t, avg_P_ev_opt_down_t);
end

%% 7. 结果可视化 (分开存储为4个高清中文图)
time_axis = (1:T)'; 
font_name = 'Microsoft YaHei'; % 设置支持中文的字体
font_size_axis = 16;
font_size_label = 18;
font_size_legend = 16;
line_width = 2.0;

% --- 图1：上调功率对比 ---
figure('Position', [100, 100, 1000, 700], 'Color', 'w');
plot(time_axis, Total_Up_Aggregated_Raw, 'b--', 'LineWidth', line_width, 'DisplayName', '原始聚合上调功率'); hold on;
plot(time_axis, Total_Up_Optimal_Agg, 'r-', 'LineWidth', line_width, 'DisplayName', '优化后上调功率');
plot(time_axis, PV_Up_Optimal_Agg, 'g-.', 'LineWidth', line_width, 'DisplayName', '光伏PV贡献'); % 新增PV曲线
plot(time_axis, P_grid_up_demand, 'k:', 'LineWidth', line_width, 'DisplayName', '上调需求');
hold off;
% 无标题
xlabel('时间步', 'FontSize', font_size_label, 'FontName', font_name); 
ylabel('功率 (kW)', 'FontSize', font_size_label, 'FontName', font_name);
legend('show', 'Location', 'best', 'FontSize', font_size_legend, 'FontName', font_name); 
grid on;
set(gca, 'FontSize', font_size_axis, 'FontName', font_name);
print('上调功率对比.png', '-dpng', '-r600');
fprintf('已保存: 上调功率对比.png\n');

% --- 图2：下调功率对比 ---
figure('Position', [150, 150, 1000, 700], 'Color', 'w');
plot(time_axis, Total_Down_Aggregated_Raw, 'b--', 'LineWidth', line_width, 'DisplayName', '原始聚合下调功率'); hold on;
plot(time_axis, Total_Down_Optimal_Agg, 'r-', 'LineWidth', line_width, 'DisplayName', '优化后下调功率');
plot(time_axis, PV_Down_Optimal_Agg, 'g-.', 'LineWidth', line_width, 'DisplayName', '光伏PV贡献'); % 新增PV曲线
plot(time_axis, P_grid_down_demand, 'k:', 'LineWidth', line_width, 'DisplayName', '下调需求');
hold off;
% 无标题
xlabel('时间步', 'FontSize', font_size_label, 'FontName', font_name); 
ylabel('功率 (kW)', 'FontSize', font_size_label, 'FontName', font_name);
legend('show', 'Location', 'best', 'FontSize', font_size_legend, 'FontName', font_name); 
grid on;
set(gca, 'FontSize', font_size_axis, 'FontName', font_name);
print('下调功率对比.png', '-dpng', '-r600');
fprintf('已保存: 下调功率对比.png\n');

% --- 图3：上调互补性与相关性 (保持 AC vs EV) ---
figure('Position', [200, 200, 1000, 700], 'Color', 'w');
indicator_names_up = {'SDCI⁺', 'Spearman ρ⁺'};
values_raw_up = [SDCI_up_raw; rho_up_raw];
values_opt_up = [SDCI_up_opt; rho_up_opt];
bar_data_up = [values_raw_up, values_opt_up];

b_up = bar(bar_data_up);
set(gca, 'XTickLabel', indicator_names_up, 'FontSize', font_size_axis, 'FontName', font_name);
ylabel('指标值 (AC vs EV)', 'FontSize', font_size_label, 'FontName', font_name); 
ylim([-1.2, 1.2]); 
legend([b_up(1) b_up(2)], {'原始聚合', '优化后'}, 'Location', 'north', 'Orientation','horizontal', 'FontSize', font_size_legend, 'FontName', font_name);
% 无标题
grid on;

% 添加数值标签
for k_bar = 1:size(bar_data_up,1) 
    for j_bar = 1:size(bar_data_up,2)
        text(b_up(j_bar).XData(k_bar) + b_up(j_bar).XOffset, bar_data_up(k_bar,j_bar), sprintf('%.3f', bar_data_up(k_bar,j_bar)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize', 12, 'FontName', font_name, 'Color','k');
    end
end
set(gca, 'FontSize', font_size_axis, 'FontName', font_name);
print('上调互补性与相关性.png', '-dpng', '-r600');
fprintf('已保存: 上调互补性与相关性.png\n');

% --- 图4：下调互补性与相关性 (保持 AC vs EV) ---
figure('Position', [250, 250, 1000, 700], 'Color', 'w');
indicator_names_down = {'SDCI⁻', 'Spearman ρ⁻'};
values_raw_down = [SDCI_down_raw; rho_down_raw];
values_opt_down = [SDCI_down_opt; rho_down_opt];
bar_data_down = [values_raw_down, values_opt_down];

b_down = bar(bar_data_down);
set(gca, 'XTickLabel', indicator_names_down, 'FontSize', font_size_axis, 'FontName', font_name);
ylabel('指标值 (AC vs EV)', 'FontSize', font_size_label, 'FontName', font_name); 
ylim([-1.2, 1.2]);
legend([b_down(1) b_down(2)], {'原始聚合', '优化后'}, 'Location', 'north', 'Orientation','horizontal', 'FontSize', font_size_legend, 'FontName', font_name);
% 无标题
grid on;

% 添加数值标签
for k_bar = 1:size(bar_data_down,1)
    for j_bar = 1:size(bar_data_down,2)
        text(b_down(j_bar).XData(k_bar) + b_down(j_bar).XOffset, bar_data_down(k_bar,j_bar), sprintf('%.3f', bar_data_down(k_bar,j_bar)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize', 12, 'FontName', font_name, 'Color','k');
    end
end
set(gca, 'FontSize', font_size_axis, 'FontName', font_name);
print('下调互补性与相关性.png', '-dpng', '-r600');
fprintf('已保存: 下调互补性与相关性.png\n');

%% 8. 命令行输出汇总
disp(' ');
disp('=== 互补性与相关性指标对比 (GA 优化结果 - 含光伏) ===');
disp('【上调】');
fprintf('优化前 (原始聚合): SDCI⁺ = %.4f, ρ⁺ = %.4f\n', SDCI_up_raw, rho_up_raw);
fprintf('优化后 (GA 优化):  SDCI⁺ = %.4f, ρ⁺ = %.4f, 优化成本 = %.2f\n', SDCI_up_opt, rho_up_opt, cost_up_optimal);
disp(' ');
disp('【下调】');
fprintf('优化前 (原始聚合): SDCI⁻ = %.4f, ρ⁻ = %.4f\n', SDCI_down_raw, rho_down_raw);
fprintf('优化后 (GA 优化):  SDCI⁻ = %.4f, ρ⁻ = %.4f, 优化成本 = %.2f\n', SDCI_down_opt, rho_down_opt, cost_down_optimal);
disp(' ');