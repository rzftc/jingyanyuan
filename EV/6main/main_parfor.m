%% main_parfortest.m - (parfor 加速 - 已修正句柄问题 - 最精简注释版)
clc; clear; close all;
rng(2024);

%% ===================== 初始化参数 =====================
excelFile = '../0inputdata/residential_all_models.xlsx';
if ~exist(excelFile, 'file')
    generateEVParameters_real(excelFile, 100, 0.6);
    fprintf('已生成参数模板: %s\n', excelFile);
end
[EVs, t_sim, dt_short, dt_long, P_tar] = initializeFromExcel(excelFile);
fprintf('成功加载%d辆EV数据\n', length(EVs));

%% ===================== 索引预处理 =====================
num_long_steps = t_sim / dt_long;
num_short_per_long = dt_long / dt_short;
total_steps = num_long_steps * num_short_per_long;
assert(mod(num_long_steps, 24) == 0, '长时间步数必须是24的整数倍');

%% ===================== 创建结果存储结构 =====================
num_evs = length(EVs);
results = struct(...
    'P_agg',         zeros(1, total_steps), ...
    'P_base',        zeros(1, total_steps), ...
    'S_agg',         zeros(1, total_steps), ...
    'lambda',        zeros(1, total_steps), ...
    'EV_S_original', zeros(num_evs, total_steps), ...
    'EV_S_mod',      zeros(num_evs, total_steps), ...
    'm3',            zeros(num_evs, 1), ...
    'P_tar',         zeros(1, total_steps), ...
    'P_cu',          zeros(1, total_steps) ...
);

%% ===================== 初始化EV目标充电功率 (向量化) =====================
p_real_vec = [EVs.p_real]';
P_h_max_vec = [EVs.P_h_max]';
P_0_vec = [EVs.P_0]';
P_l_min_vec = [EVs.P_l_min]';
Delta_E_h_max_vec = [EVs.Delta_E_h_max]';
Delta_E_q_max_vec = [EVs.Delta_E_q_max]';
E_tar_set_vec = [EVs.E_tar_set]';
E_ini_vec = [EVs.E_ini]';
P_N_vec = [EVs.P_N]';

Delta_E_vec = calculateDeltaE_vec(p_real_vec, P_h_max_vec, P_0_vec, P_l_min_vec, Delta_E_h_max_vec, Delta_E_q_max_vec);
E_tar_vec = max(E_tar_set_vec - Delta_E_vec, E_ini_vec);

t_ch_vec = zeros(num_evs, 1);
valid_PN_mask = P_N_vec > 0;
if any(valid_PN_mask)
    t_ch_vec(valid_PN_mask) = 60 .* (E_tar_vec(valid_PN_mask) - E_ini_vec(valid_PN_mask)) ./ P_N_vec(valid_PN_mask);
end
t_ch_vec = max(t_ch_vec, 0);

E_tar_cell = num2cell(E_tar_vec);
[EVs.E_tar] = E_tar_cell{:};
tau_rem_cell = num2cell(t_ch_vec);
[EVs.tau_rem] = tau_rem_cell{:};
E_ini_cell = num2cell(E_ini_vec);
[EVs.E_exp] = E_ini_cell{:};
[EVs.E_actual] = E_ini_cell{:};

%% ===================== 初始化聚合SOC =====================
S_agg_current = mean([EVs.S_original]);
[~, P_base_opt] = EVbaseP_aggregate_short(EVs, S_agg_current, 48, dt_long);

%% ===================== 外层循环（长时间步长） =====================
for long_idx = 1:num_long_steps
    t_long = (long_idx - 1) * dt_long;

    %% --------------------- 长时间步处理 ---------------------
    EVs = distributeBasePower(EVs, 1, P_base_opt);
    [lambda_star] = aggregateEVs(EVs, P_tar(long_idx));
    [~, S_agg_next] = calculateVirtualSOC_agg(EVs, dt_long);

    %% ===================== 内层循环（短时间步长） =====================
    for short_idx = 1:num_short_per_long
        step_idx = (long_idx - 1) * num_short_per_long + short_idx;
        t_current = t_long + (short_idx - 1) * dt_short;

        temp_m3 = zeros(num_evs, 1);
        temp_S_original = zeros(num_evs, 1);
        temp_S_mod = zeros(num_evs, 1);
        temp_P_current = zeros(num_evs, 1);

        %% --------------------- 更新EV状态 (并行处理) ---------------------
        parfor i = 1:num_evs % 使用 parfor 加速
            EV = EVs(i);
            EV = updateLockState(EV, t_current);

            EV_for_handle = EV;
            EV_temp_with_handle = generateDemandCurve(EV_for_handle);

            current_P_val = 0;
            if isfield(EV_temp_with_handle, 'demandCurve') && isa(EV_temp_with_handle.demandCurve, 'function_handle')
                current_P_val = EV_temp_with_handle.demandCurve(lambda_star);
            else
                 switch EV.state
                    case 'LockON'
                        current_P_val = EV.P_N;
                    case {'LockOFF', 'OFF'}
                        current_P_val = 0;
                    otherwise
                        current_P_val = 0;
                 end
            end
            EV.P_current = current_P_val;

            EV = calculateVirtualSOC_upgrade(EV, t_current, dt_short);

            if EV.t_dep > EV.t_in
                 m3_val = (EV.E_tar - EV.E_ini) / (EV.eta * ((EV.t_dep - EV.t_in) / 60));
            else
                 m3_val = 0;
            end

            temp_m3(i) = m3_val;
            temp_S_original(i) = EV.S_original;
            temp_S_mod(i) = EV.S_modified;
            temp_P_current(i) = EV.P_current;

            EVs(i) = EV;

        end % 结束 parfor i

        results.EV_S_original(:, step_idx) = temp_S_original;
        results.EV_S_mod(:, step_idx) = temp_S_mod;
        results.m3 = temp_m3;

        results.lambda(step_idx) = lambda_star;
        results.S_agg(step_idx) = S_agg_current;
        results.P_agg(step_idx) = sum(temp_P_current);

        results.P_cu(step_idx) = temp_P_current(10);

    end % 结束 short_idx

    %% --------------------- 更新聚合SOC ---------------------
    S_agg_current = S_agg_next;

end % 结束 long_idx

results.P_tar = repelem(P_tar, num_short_per_long);

%% ===================== 可视化结果 =====================
selected_ev = 10;
time_min = (0:total_steps-1) * dt_short;
time_h = time_min / 60;

S_original = results.EV_S_original(selected_ev, :);
S_mod = results.EV_S_mod(selected_ev, :);
P_current_plot = results.P_cu;

figure('Position', [100 100 1200 600]);

yyaxis left;
plot(time_h, S_original, 'b-', 'LineWidth', 1.5, 'DisplayName', '原始SOC (S_{original})');
hold on;
plot(time_h, S_mod, 'r--', 'LineWidth', 1.5, 'DisplayName', '修正SOC (S_{modified})');
ylabel('SOC');
ylim([-1.1 1.1]);
grid on;

yyaxis right;
plot(time_h, P_current_plot, 'g-.', 'LineWidth', 1.5, 'DisplayName', ['EV ', num2str(selected_ev), ' 功率 (P_{current})']);
ylabel('功率 (kW)');
P_lim = [min(P_current_plot)-1, max(P_current_plot)+1];
ylim(P_lim);

xlabel('时间 (小时)');
title(['第', num2str(selected_ev), '台EV的SOC与当前功率对比']);
legend('Location', 'best');
set(gca, 'FontSize', 12);
hold off;

% plotLambdaAndAggSOC(results, dt_short);
% plotPowerComparison(results, dt_short);