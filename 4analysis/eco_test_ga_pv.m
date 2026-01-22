%% 测试脚本：调用 solve_total_time_dispatch_ga_pv (含光伏特性)
clc; close all; clear;

%% 1. 数据加载与参数初始化
load('results_single.mat'); 

P_ac_up_individual = results.AC_Up_Individual;
P_ac_down_individual = abs(results.AC_Down_Individual); 
P_ev_up_individual = results.EV_Up_Individual;
P_ev_down_individual = abs(results.EV_Down_Individual); 

num_ac_total = size(P_ac_up_individual, 1);
num_ev_total = size(P_ev_up_individual, 1);
T = size(P_ac_up_individual, 2); 

% --- 生成光伏(PV)模拟数据 (基于昼夜特性) ---
num_pv_total = round((num_ac_total + num_ev_total) / 2); 
if num_pv_total == 0; num_pv_total = 10; end

% 计算功率基准 (与AC/EV同量级)
avg_ac_mag = mean(mean(P_ac_up_individual(:)) + mean(P_ac_down_individual(:))) / 2;
avg_ev_mag = mean(mean(P_ev_up_individual(:)) + mean(P_ev_down_individual(:))) / 2;
base_pv_peak_mag = (avg_ac_mag + avg_ev_mag) / 2; 

% 生成光照曲线 (6:00-18:00 有光照)
t_hours = linspace(0, 24, T); 
solar_profile = max(0, sin((t_hours - 6) * pi / 12)); 
solar_noise = 0.9 + 0.2 * rand(1, T); 
solar_profile = solar_profile .* solar_noise;

% 生成个体PV调节能力
P_pv_up_individual = zeros(num_pv_total, T);
P_pv_down_individual = zeros(num_pv_total, T);

for i = 1:num_pv_total
    ind_scale = 0.8 + 0.4 * rand();
    % 下调能力: 对应当前最大出力
    P_pv_down_individual(i, :) = base_pv_peak_mag * ind_scale * solar_profile;
    % 上调能力: 假设预留40%备用
    reserve_ratio = 0.4; 
    P_pv_up_individual(i, :) = base_pv_peak_mag * ind_scale * solar_profile * reserve_ratio;
end

% --- 生成/加载电网调节需求 ---
if isfield(results, 'P_grid_up_regulation_demand') && ~isempty(results.P_grid_up_regulation_demand)
    P_grid_up_demand = results.P_grid_up_regulation_demand(:);
else
    avg_power = (sum(sum(P_ac_up_individual,1)) + sum(sum(P_ev_up_individual,1)) + sum(sum(P_pv_up_individual,1))) / T;
    if avg_power == 0; avg_power = 100; end
    P_grid_up_demand = avg_power * (0.2 + 0.3 * rand(T,1));
end

if isfield(results, 'P_grid_down_regulation_demand') && ~isempty(results.P_grid_down_regulation_demand)
    P_grid_down_demand = results.P_grid_down_regulation_demand(:);
else
    avg_power = (sum(sum(P_ac_down_individual,1)) + sum(sum(P_ev_down_individual,1)) + sum(sum(P_pv_down_individual,1))) / T;
    if avg_power == 0; avg_power = 80; end
    P_grid_down_demand = avg_power * (0.15 + 0.25 * sin(linspace(0,2*pi,T)'*0.5 + pi/4).^2 );
end

eps_val = 1e-6;

%% 2. 成本参数设置
C_ac_up = 0.05 * ones(num_ac_total, 1); 
C_ev_up = 0.04 * ones(num_ev_total, 1);
C_pv_up = 0.03 * ones(num_pv_total, 1); 

C_ac_down = 0.03 * ones(num_ac_total, 1);
C_ev_down = 0.02 * ones(num_ev_total, 1);
C_pv_down = 0.015 * ones(num_pv_total, 1); 

%% 3. 计算原始聚合指标 (无优化)
% 上调
AC_Up_Raw = sum(P_ac_up_individual, 1)';
EV_Up_Raw = sum(P_ev_up_individual, 1)';
PV_Up_Raw = sum(P_pv_up_individual, 1)';
Total_Up_Raw = AC_Up_Raw + EV_Up_Raw + PV_Up_Raw;

n_ac = repmat(num_ac_total, T, 1); avg_ac_up = AC_Up_Raw ./ (n_ac + eps_val);
n_ev = repmat(num_ev_total, T, 1); avg_ev_up = EV_Up_Raw ./ (n_ev + eps_val);

if sum(n_ac) == 0 && sum(n_ev) == 0
    SDCI_up_raw = 0; rho_up_raw = 0;
else
    SDCI_up_raw = calculateSDCI(n_ac, n_ev, avg_ac_up, avg_ev_up);
    rho_up_raw = calculateSpearmanRho(n_ac, avg_ac_up, n_ev, avg_ev_up);
end

% 下调
AC_Down_Raw = sum(P_ac_down_individual, 1)';
EV_Down_Raw = sum(P_ev_down_individual, 1)';
PV_Down_Raw = sum(P_pv_down_individual, 1)';
Total_Down_Raw = AC_Down_Raw + EV_Down_Raw + PV_Down_Raw;

avg_ac_down = AC_Down_Raw ./ (n_ac + eps_val);
avg_ev_down = EV_Down_Raw ./ (n_ev + eps_val);

if sum(n_ac) == 0 && sum(n_ev) == 0
    SDCI_down_raw = 0; rho_down_raw = 0;
else
    SDCI_down_raw = calculateSDCI(n_ac, n_ev, avg_ac_down, avg_ev_down);
    rho_down_raw = calculateSpearmanRho(n_ac, avg_ac_down, n_ev, avg_ev_down);
end

%% 4. GA 参数配置
N_vars = (num_ac_total + num_ev_total + num_pv_total) * T; 
ga_opts = optimoptions('ga', ...
    'PopulationType', 'bitstring', ...
    'PopulationSize', min(50, max(50, round(N_vars / 1000 * 50))), ... 
    'MaxGenerations', min(100, max(100, round(N_vars / 1000 * 100))),... 
    'EliteCount', ceil(0.05 * min(200, max(50, round(N_vars / 1000 * 50)))), ...
    'CrossoverFraction', 0.8, ...
    'MutationFcn', {@mutationuniform, 0.01}, ...
    'Display', 'iter', ... 
    'UseParallel', true,...
    'PlotFcn', @gaplotbestf); 

%% 5. 执行优化
% --- 上调优化 ---
U_ac_up_opt = zeros(num_ac_total, T); 
U_ev_up_opt = zeros(num_ev_total, T);
U_pv_up_opt = zeros(num_pv_total, T);
cost_up_opt = 0; 

if (num_ac_total > 0) || (num_ev_total > 0) || (num_pv_total > 0)
    [U_ac_up_opt, U_ev_up_opt, U_pv_up_opt, cost_up_opt, ~, ~] = ... 
        solve_total_time_dispatch_ga_pv(num_ac_total, num_ev_total, num_pv_total, ...
                                     P_ac_up_individual, P_ev_up_individual, P_pv_up_individual, ...
                                     C_ac_up, C_ev_up, C_pv_up, ...
                                     P_grid_up_demand, T, ...
                                     SDCI_up_raw, rho_up_raw, ga_opts);
end

% --- 下调优化 ---
U_ac_down_opt = zeros(num_ac_total, T);
U_ev_down_opt = zeros(num_ev_total, T);
U_pv_down_opt = zeros(num_pv_total, T);
cost_down_opt = 0;

if (num_ac_total > 0) || (num_ev_total > 0) || (num_pv_total > 0)
    [U_ac_down_opt, U_ev_down_opt, U_pv_down_opt, cost_down_opt, ~, ~] = ...
        solve_total_time_dispatch_ga_pv(num_ac_total, num_ev_total, num_pv_total, ...
                                     P_ac_down_individual, P_ev_down_individual, P_pv_down_individual, ...
                                     C_ac_down, C_ev_down, C_pv_down, ...
                                     P_grid_down_demand, T, ...
                                     SDCI_down_raw, rho_down_raw, ga_opts); 
end

%% 6. 计算优化后聚合结果
% 上调
AC_Up_Opt = sum(P_ac_up_individual .* U_ac_up_opt, 1)';
EV_Up_Opt = sum(P_ev_up_individual .* U_ev_up_opt, 1)';
PV_Up_Opt = sum(P_pv_up_individual .* U_pv_up_opt, 1)';
Total_Up_Opt = AC_Up_Opt + EV_Up_Opt + PV_Up_Opt;

n_ac_opt = sum(U_ac_up_opt, 1)'; avg_ac_opt = AC_Up_Opt ./ (n_ac_opt + eps_val);
n_ev_opt = sum(U_ev_up_opt, 1)'; avg_ev_opt = EV_Up_Opt ./ (n_ev_opt + eps_val);

if sum(n_ac_opt) == 0 && sum(n_ev_opt) == 0
    SDCI_up_opt = 0; rho_up_opt = 0;
else
    SDCI_up_opt = calculateSDCI(n_ac_opt, n_ev_opt, avg_ac_opt, avg_ev_opt);
    rho_up_opt = calculateSpearmanRho(n_ac_opt, avg_ac_opt, n_ev_opt, avg_ev_opt);
end

% 下调
AC_Down_Opt = sum(P_ac_down_individual .* U_ac_down_opt, 1)';
EV_Down_Opt = sum(P_ev_down_individual .* U_ev_down_opt, 1)';
PV_Down_Opt = sum(P_pv_down_individual .* U_pv_down_opt, 1)';
Total_Down_Opt = AC_Down_Opt + EV_Down_Opt + PV_Down_Opt;

n_ac_opt_d = sum(U_ac_down_opt, 1)'; avg_ac_opt_d = AC_Down_Opt ./ (n_ac_opt_d + eps_val);
n_ev_opt_d = sum(U_ev_down_opt, 1)'; avg_ev_opt_d = EV_Down_Opt ./ (n_ev_opt_d + eps_val);

if sum(n_ac_opt_d) == 0 && sum(n_ev_opt_d) == 0
    SDCI_down_opt = 0; rho_down_opt = 0;
else
    SDCI_down_opt = calculateSDCI(n_ac_opt_d, n_ev_opt_d, avg_ac_opt_d, avg_ev_opt_d);
    rho_down_opt = calculateSpearmanRho(n_ac_opt_d, avg_ac_opt_d, n_ev_opt_d, avg_ev_opt_d);
end

%% 7. 绘图与输出
time = (1:T)'; 
font = 'Microsoft YaHei'; 
lw = 2.0;

% 图1：上调
figure('Position', [100, 100, 1000, 700], 'Color', 'w');
plot(time, Total_Up_Raw, 'b--', 'LineWidth', lw, 'DisplayName', '原始聚合上调'); hold on;
plot(time, Total_Up_Opt, 'r-', 'LineWidth', lw, 'DisplayName', '优化后上调');
plot(time, PV_Up_Opt, 'g-.', 'LineWidth', lw, 'DisplayName', '光伏PV贡献');
plot(time, P_grid_up_demand, 'k:', 'LineWidth', lw, 'DisplayName', '上调需求');
xlabel('时间步'); ylabel('功率 (kW)'); legend('show'); grid on; set(gca,'FontSize',14,'FontName',font);
print('上调功率对比.png', '-dpng', '-r600');

% 图2：下调
figure('Position', [150, 150, 1000, 700], 'Color', 'w');
plot(time, Total_Down_Raw, 'b--', 'LineWidth', lw, 'DisplayName', '原始聚合下调'); hold on;
plot(time, Total_Down_Opt, 'r-', 'LineWidth', lw, 'DisplayName', '优化后下调');
plot(time, PV_Down_Opt, 'g-.', 'LineWidth', lw, 'DisplayName', '光伏PV贡献');
plot(time, P_grid_down_demand, 'k:', 'LineWidth', lw, 'DisplayName', '下调需求');
xlabel('时间步'); ylabel('功率 (kW)'); legend('show'); grid on; set(gca,'FontSize',14,'FontName',font);
print('下调功率对比.png', '-dpng', '-r600');

% 图3：上调指标
figure('Position', [200, 200, 1000, 700], 'Color', 'w');
bar_data = [[SDCI_up_raw; rho_up_raw], [SDCI_up_opt; rho_up_opt]];
b = bar(bar_data);
set(gca, 'XTickLabel', {'SDCI⁺', 'Spearman ρ⁺'}, 'FontSize', 14, 'FontName', font);
legend({'原始', '优化后'}, 'Location', 'north'); grid on; ylim([-1.2, 1.2]);
% 添加标签
for k = 1:2
    for j = 1:2
        text(b(j).XData(k)+b(j).XOffset, bar_data(k,j), sprintf('%.3f', bar_data(k,j)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize', 12);
    end
end
print('上调互补性与相关性.png', '-dpng', '-r600');

% 图4：下调指标
figure('Position', [250, 250, 1000, 700], 'Color', 'w');
bar_data_d = [[SDCI_down_raw; rho_down_raw], [SDCI_down_opt; rho_down_opt]];
b = bar(bar_data_d);
set(gca, 'XTickLabel', {'SDCI⁻', 'Spearman ρ⁻'}, 'FontSize', 14, 'FontName', font);
legend({'原始', '优化后'}, 'Location', 'north'); grid on; ylim([-1.2, 1.2]);
for k = 1:2
    for j = 1:2
        text(b(j).XData(k)+b(j).XOffset, bar_data_d(k,j), sprintf('%.3f', bar_data_d(k,j)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize', 12);
    end
end
print('下调互补性与相关性.png', '-dpng', '-r600');

% 命令行输出
fprintf('\n=== 结果汇总 (含PV) ===\n');
fprintf('上调: 原始 SDCI=%.4f, rho=%.4f | 优化后 SDCI=%.4f, rho=%.4f, Cost=%.2f\n', ...
    SDCI_up_raw, rho_up_raw, SDCI_up_opt, rho_up_opt, cost_up_opt);
fprintf('下调: 原始 SDCI=%.4f, rho=%.4f | 优化后 SDCI=%.4f, rho=%.4f, Cost=%.2f\n', ...
    SDCI_down_raw, rho_down_raw, SDCI_down_opt, rho_down_opt, cost_down_opt);