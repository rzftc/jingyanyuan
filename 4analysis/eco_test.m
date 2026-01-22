%% 测试脚本
clc; close all; clear;

%% 1. 加载数据
load('results_single.mat'); 

% --- 数据提取 ---
P_ac_up_individual = []; P_ev_up_individual = [];
P_ac_down_individual = []; P_ev_down_individual = [];
T = 0; num_ac_total = 0; num_ev_total = 0;

if isfield(results, 'AC_Up_Individual') && ~isempty(results.AC_Up_Individual)
    P_ac_up_individual = results.AC_Up_Individual;
    num_ac_total = size(P_ac_up_individual, 1);
    T = size(P_ac_up_individual, 2);
    if isfield(results, 'AC_Down_Individual') && isequal(size(results.AC_Down_Individual), [num_ac_total, T])
        P_ac_down_individual = abs(results.AC_Down_Individual);
    else
        P_ac_down_individual = zeros(num_ac_total, T); 
    end
else
    if isfield(results, 'EV_Up_Individual') && ~isempty(results.EV_Up_Individual)
        T = size(results.EV_Up_Individual, 2); % 尝试从EV获取T
    end
    P_ac_up_individual = zeros(0,T); P_ac_down_individual = zeros(0,T);
end

if isfield(results, 'EV_Up_Individual') && ~isempty(results.EV_Up_Individual)
    P_ev_up_individual = results.EV_Up_Individual;
    num_ev_total = size(P_ev_up_individual, 1);
    if T == 0; T = size(P_ev_up_individual, 2); end
    if T ~= size(P_ev_up_individual, 2); error('AC和EV数据时间步T不一致'); end
    if isfield(results, 'EV_Down_Individual') && isequal(size(results.EV_Down_Individual), [num_ev_total, T])
        P_ev_down_individual = abs(results.EV_Down_Individual);
    else
        P_ev_down_individual = zeros(num_ev_total, T); 
    end
else
    P_ev_up_individual = zeros(0,T); P_ev_down_individual = zeros(0,T);
end

if T == 0 && (num_ac_total > 0 || num_ev_total > 0)
     error('未能从个体能力数据中确定总时间步数 T，但存在设备数据。');
elseif T == 0 
    disp('警告：没有时间步数据 (T=0)，后续计算将跳过或产生空结果。');
end

AC_Up_raw = results.AC_Up(:)'; EV_Up_raw = results.EV_Up(:)';
AC_Down_raw = abs(results.AC_Down(:)'); EV_Down_raw = abs(results.EV_Down(:)');
if T > 0 % 仅当T有效时检查长度
    if length(AC_Up_raw) ~= T; AC_Up_raw = zeros(1,T); end
    if length(EV_Up_raw) ~= T; EV_Up_raw = zeros(1,T); end
    if length(AC_Down_raw) ~= T; AC_Down_raw = zeros(1,T); end
    if length(EV_Down_raw) ~= T; EV_Down_raw = zeros(1,T); end
end

fprintf('总细分时间步数 T = %d\n', T);
fprintf('总空调数 num_ac_total = %d\n', num_ac_total);
fprintf('总电动汽车数 num_ev_total = %d\n', num_ev_total);

%% 2. 定义仿真参数
dt = 0.05; 
time_axis = [];
if T > 0; time_axis = (0:T-1) * dt; end

c_ac_cost = ones(num_ac_total, 1) * 0.1;  
c_ev_cost = ones(num_ev_total, 1) * 0.08; 

rng('default'); 
P_req_t_series_up = zeros(1,T); P_req_t_series_down = zeros(1,T);
if T > 0
    P_req_base_up = 10 + 5 * sin(2 * pi * time_axis / 24) + randn(1, T) * 2;
    P_req_t_series_up = max(0, P_req_base_up); 
    P_req_base_down = 8 + 4 * cos(2 * pi * time_axis / 24 + pi/2) + randn(1, T) * 2;
    P_req_t_series_down = max(0, P_req_base_down);
end

%% 3. 逐时段优化
u_ac_opt_up_t = zeros(num_ac_total, T); P_ac_dispatched_opt_up_t = zeros(1,T);
u_ev_opt_up_t = zeros(num_ev_total, T); P_ev_dispatched_opt_up_t = zeros(1,T);
optimal_cost_up_t = zeros(1, T); 

u_ac_opt_down_t = zeros(num_ac_total, T); P_ac_dispatched_opt_down_t = zeros(1,T);
u_ev_opt_down_t = zeros(num_ev_total, T); P_ev_dispatched_opt_down_t = zeros(1,T);
optimal_cost_down_t = zeros(1, T); 

if T > 0
    fprintf('\n开始逐时段优化 (最小化成本)...\n');
    for t = 1:T
        % 上调
        p_ac_t_up = []; if num_ac_total > 0; p_ac_t_up = P_ac_up_individual(:, t); end
        p_ev_t_up = []; if num_ev_total > 0; p_ev_t_up = P_ev_up_individual(:, t); end
        if num_ac_total > 0 || num_ev_total > 0
            [u_ac, u_ev, cost_val, flag] = solve_hourly_dispatch(num_ac_total, num_ev_total, p_ac_t_up, p_ev_t_up, c_ac_cost, c_ev_cost, P_req_t_series_up(t));
            if flag > 0
                if num_ac_total > 0; u_ac_opt_up_t(:, t) = u_ac; P_ac_dispatched_opt_up_t(t) = sum(u_ac .* p_ac_t_up); end
                if num_ev_total > 0; u_ev_opt_up_t(:, t) = u_ev; P_ev_dispatched_opt_up_t(t) = sum(u_ev .* p_ev_t_up); end
                optimal_cost_up_t(t) = cost_val;
            else; optimal_cost_up_t(t) = NaN; end
        else; optimal_cost_up_t(t) = 0; end

        % 下调
        p_ac_t_down = []; if num_ac_total > 0; p_ac_t_down = P_ac_down_individual(:, t); end
        p_ev_t_down = []; if num_ev_total > 0; p_ev_t_down = P_ev_down_individual(:, t); end
        if num_ac_total > 0 || num_ev_total > 0
            [u_ac, u_ev, cost_val, flag] = solve_hourly_dispatch(num_ac_total, num_ev_total, p_ac_t_down, p_ev_t_down, c_ac_cost, c_ev_cost, P_req_t_series_down(t));
            if flag > 0
                if num_ac_total > 0; u_ac_opt_down_t(:, t) = u_ac; P_ac_dispatched_opt_down_t(t) = sum(u_ac .* p_ac_t_down); end
                if num_ev_total > 0; u_ev_opt_down_t(:, t) = u_ev; P_ev_dispatched_opt_down_t(t) = sum(u_ev .* p_ev_t_down); end
                optimal_cost_down_t(t) = cost_val;
            else; optimal_cost_down_t(t) = NaN; end
        else; optimal_cost_down_t(t) = 0; end
    end
    fprintf('逐时段优化完成。\n');
end

%% 4. 计算优化前后的SDCI和Rho指标

% --- "优化前" (Raw/Benchmark) 场景指标计算 ---
n_ac_for_raw_calc = ones(T,1) * num_ac_total; 
n_ev_for_raw_calc = ones(T,1) * num_ev_total;
if num_ac_total == 0; n_ac_for_raw_calc = zeros(T,1); end 
if num_ev_total == 0; n_ev_for_raw_calc = zeros(T,1); end

avg_deltaP_ac_raw_up = zeros(T,1); avg_deltaP_ev_raw_up = zeros(T,1);
avg_deltaP_ac_raw_down = zeros(T,1); avg_deltaP_ev_raw_down = zeros(T,1);
if T > 0
    avg_deltaP_ac_raw_up = AC_Up_raw' ./ (n_ac_for_raw_calc + eps);
    avg_deltaP_ev_raw_up = EV_Up_raw' ./ (n_ev_for_raw_calc + eps);
    avg_deltaP_ac_raw_down = AC_Down_raw' ./ (n_ac_for_raw_calc + eps);
    avg_deltaP_ev_raw_down = EV_Down_raw' ./ (n_ev_for_raw_calc + eps);
end

SDCI_up_raw = ensureScalar(calculateSDCI(n_ac_for_raw_calc, n_ev_for_raw_calc, avg_deltaP_ac_raw_up, avg_deltaP_ev_raw_up));
rho_up_raw  = ensureScalar(calculateSpearmanRho(n_ac_for_raw_calc, avg_deltaP_ac_raw_up, n_ev_for_raw_calc, avg_deltaP_ev_raw_up));
SDCI_down_raw = ensureScalar(calculateSDCI(n_ac_for_raw_calc, n_ev_for_raw_calc, avg_deltaP_ac_raw_down, avg_deltaP_ev_raw_down));
rho_down_raw  = ensureScalar(calculateSpearmanRho(n_ac_for_raw_calc, avg_deltaP_ac_raw_down, n_ev_for_raw_calc, avg_deltaP_ev_raw_down));

% --- "优化后"场景的SDCI和Rho计算 ---
n_ac_opt_up_count_t = sum(u_ac_opt_up_t, 1); 
n_ev_opt_up_count_t = sum(u_ev_opt_up_t, 1);
n_ac_opt_down_count_t = sum(u_ac_opt_down_t, 1);
n_ev_opt_down_count_t = sum(u_ev_opt_down_t, 1);

avg_deltaP_ac_opt_up = P_ac_dispatched_opt_up_t ./ (n_ac_opt_up_count_t + eps);
avg_deltaP_ev_opt_up = P_ev_dispatched_opt_up_t ./ (n_ev_opt_up_count_t + eps);
avg_deltaP_ac_opt_down = P_ac_dispatched_opt_down_t ./ (n_ac_opt_down_count_t + eps);
avg_deltaP_ev_opt_down = P_ev_dispatched_opt_down_t ./ (n_ev_opt_down_count_t + eps);

SDCI_up_opt = ensureScalar(calculateSDCI(n_ac_opt_up_count_t', n_ev_opt_up_count_t', avg_deltaP_ac_opt_up', avg_deltaP_ev_opt_up'));
rho_up_opt  = ensureScalar(calculateSpearmanRho(n_ac_opt_up_count_t', avg_deltaP_ac_opt_up', n_ev_opt_up_count_t', avg_deltaP_ev_opt_up'));
SDCI_down_opt = ensureScalar(calculateSDCI(n_ac_opt_down_count_t', n_ev_opt_down_count_t', avg_deltaP_ac_opt_down', avg_deltaP_ev_opt_down'));
rho_down_opt  = ensureScalar(calculateSpearmanRho(n_ac_opt_down_count_t', avg_deltaP_ac_opt_down', n_ev_opt_down_count_t', avg_deltaP_ev_opt_down'));

%% 5. 结果分析与可视化
figure('Name', 'VPP Dispatch & Metrics Comparison');
subplot(2,2,1);
plot(time_axis, optimal_cost_up_t, 'm-^', 'DisplayName', 'Optimal Cost (Up)'); 
hold on;
plot(time_axis, P_ac_dispatched_opt_up_t + P_ev_dispatched_opt_up_t, 'r-s', 'DisplayName', 'Dispatched Power (Up)');
plot(time_axis, P_req_t_series_up, 'k--', 'DisplayName', 'Required Power (Up)');
hold off; xlabel('Time (hours)'); ylabel('Cost / Power (kW)');
title('VPP Up-Regulation Dispatch Results'); legend('show', 'Location', 'best'); grid on;

subplot(2,2,2);
plot(time_axis, optimal_cost_down_t, 'm-^', 'DisplayName', 'Optimal Cost (Down)');
hold on;
plot(time_axis, P_ac_dispatched_opt_down_t + P_ev_dispatched_opt_down_t, 'r-s', 'DisplayName', 'Dispatched Power (Down)');
plot(time_axis, P_req_t_series_down, 'k--', 'DisplayName', 'Required Power (Down)');
hold off; xlabel('Time (hours)'); ylabel('Cost / Power (kW)');
title('VPP Down-Regulation Dispatch Results'); legend('show', 'Location', 'best'); grid on;

subplot(2,2,3)
bar_data_up = [SDCI_up_raw, SDCI_up_opt; rho_up_raw, rho_up_opt];
b_up = bar(bar_data_up, 'grouped');
set(gca, 'XTickLabel', {'SDCI⁺', 'ρ'}); ylabel('Indicator Value');
legend([b_up(1) b_up(2)], {'Raw Aggregated', 'Optimized (MinCost)'}, 'Location', 'northoutside', 'Orientation','horizontal');
title('Up-Regulation Complementarity & Correlation'); grid on; ylim([-1.1 1.1]);

subplot(2,2,4)
bar_data_down = [SDCI_down_raw, SDCI_down_opt; rho_down_raw, rho_down_opt];
b_down = bar(bar_data_down, 'grouped');
set(gca, 'XTickLabel', {'SDCI⁻', 'ρ'}); ylabel('Indicator Value');
legend([b_down(1) b_down(2)],{'Raw Aggregated', 'Optimized (MinCost)'}, 'Location', 'northoutside', 'Orientation','horizontal');
title('Down-Regulation Complementarity & Correlation'); grid on; ylim([-1.1 1.1]);

disp('=== 互补性与相关性指标对比 ===');
disp('【上调】');
fprintf('优化前 (Raw Aggregated): SDCI⁺ = %.4f, ρ⁺ = %.4f\n', SDCI_up_raw, rho_up_raw);
fprintf('优化后 (MinCost Dispatch): SDCI⁺ = %.4f, ρ⁺ = %.4f\n', SDCI_up_opt, rho_up_opt);
disp('【下调】');
fprintf('优化前 (Raw Aggregated): SDCI⁻ = %.4f, ρ⁻ = %.4f\n', SDCI_down_raw, rho_down_raw);
fprintf('优化后 (MinCost Dispatch): SDCI⁻ = %.4f, ρ⁻ = %.4f\n', SDCI_down_opt, rho_down_opt);


%% 辅助函数
function val = ensureScalar(inputVal)
    if isscalar(inputVal)
        val = inputVal;
    elseif isempty(inputVal) 
        val = NaN; 
    else
        val = mean(inputVal(:)); 
    end
end