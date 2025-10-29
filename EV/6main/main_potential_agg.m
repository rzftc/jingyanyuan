%% main_potential_agg.m (使用聚合模型和 _new 函数计算潜力 - 简化注释版)
clc; clear; close all;
rng(2024);

%% 初始化参数
excelFile = 'resi_1000.xlsx';
if ~exist(excelFile, 'file')
    generateEVParameters_real(excelFile, 1000, 1.0);
    fprintf('已生成参数模板: %s\n', excelFile);
end
[EVs, t_sim, dt_short, dt_long, P_tar] = initializeFromExcel(excelFile);
fprintf('成功加载%d辆EV数据\n', length(EVs));

%% 时间参数定义
simulation_start_hour = 6;
simulation_end_hour   = 30;
dt = dt_short / 60;       % 短步长 (小时)
dt_minutes = dt_short;    % 短步长 (分钟)
dt_long_minutes = dt_long;% 长步长 (分钟)
t_adj = dt_long / 60;     % 调节时长 (小时)

% 时间轴 (小时)
time_points_absolute = simulation_start_hour : dt : (simulation_start_hour + t_sim/60 - dt);
num_time_points = length(time_points_absolute);
fprintf('仿真时间范围: %.2f 小时 到 %.2f 小时, 步长: %.3f 小时 (%d 点)\n', ...
    time_points_absolute(1), time_points_absolute(end), dt, num_time_points);

%% 索引预处理
num_long_steps = t_sim / dt_long;
num_short_per_long = dt_long / dt_short;
total_steps = num_long_steps * num_short_per_long;
assert(mod(num_long_steps, 24*(60/dt_long)) == 0, '长时间步数必须是24小时对应步数的整数倍');

%% 创建结果存储结构
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
    'P_cu',          zeros(1, total_steps), ...
    'EV_Up',         zeros(1, total_steps), ... % 聚合上调潜力
    'EV_Down',       zeros(1, total_steps) ...  % 聚合下调潜力
    );

%% 初始化EV参数 (向量化)
p_real_vec = [EVs.p_real]';
P_h_max_vec = [EVs.P_h_max]';
P_0_vec = [EVs.P_0]';
P_l_min_vec = [EVs.P_l_min]';
Delta_E_h_max_vec = [EVs.Delta_E_h_max]';
Delta_E_q_max_vec = [EVs.Delta_E_q_max]';
E_tar_set_vec = [EVs.E_tar_set]';
E_ini_vec = [EVs.E_ini]';
P_N_vec = [EVs.P_N]';
C_vec = [EVs.C]';

E_tar_original_vec = E_tar_set_vec;
SOC_original_vec = zeros(num_evs, 1);
E_reg_min_vec = E_tar_original_vec;
E_reg_max_vec = E_tar_original_vec;

p_min = 15; p_max = 50;
p_min_prime = 10; p_max_prime = 40;
base_Price_vec = 30 * ones(num_evs, 1);
p_incentive_vec = [EVs.p_incentive]';
participation_probabilities_vec = calculateParticipation(p_incentive_vec, base_Price_vec);
ptcp_vec = (rand(num_evs, 1) < participation_probabilities_vec);

E_tar_max_flex_vec = 0.2 * C_vec;
[deltaE_up_vec, deltaE_down_vec] = incentiveTempEV_updown(p_incentive_vec, p_min, p_max, p_min_prime, p_max_prime, E_tar_max_flex_vec);

participating_indices = find(ptcp_vec);
if ~isempty(participating_indices)
    E_reg_min_vec(participating_indices) = E_tar_original_vec(participating_indices) - deltaE_down_vec(participating_indices);
    E_reg_min_vec(participating_indices) = max(E_reg_min_vec(participating_indices), E_ini_vec(participating_indices));

    E_reg_max_vec(participating_indices) = E_tar_original_vec(participating_indices) + deltaE_up_vec(participating_indices);
    E_reg_max_vec(participating_indices) = min(E_reg_max_vec(participating_indices), C_vec(participating_indices));
end

% 写回EVs结构体
E_tar_orig_cell = num2cell(E_tar_original_vec);
SOC_orig_cell = num2cell(SOC_original_vec);
ptcp_cell = num2cell(ptcp_vec);
E_reg_min_cell = num2cell(E_reg_min_vec);
E_reg_max_cell = num2cell(E_reg_max_vec);
[EVs.E_tar_original] = E_tar_orig_cell{:};
[EVs.SOC_original] = SOC_orig_cell{:};
[EVs.ptcp] = ptcp_cell{:};
[EVs.E_reg_min] = E_reg_min_cell{:};
[EVs.E_reg_max] = E_reg_max_cell{:};

% 初始化 E_tar, tau_rem, E_exp, E_actual, E_current, P_current
E_tar_vec = max(E_tar_set_vec, E_ini_vec);
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
[EVs.E_current] = E_ini_cell{:};
P_current_init_cell = num2cell(zeros(num_evs, 1));
[EVs.P_current] = P_current_init_cell{:};

%% 初始化聚合SOC 和 基线功率
S_agg_current = mean([EVs.SOC_original]);

fprintf('正在预计算所有EV的基线功率序列...\n');
EVs_for_baseline = EVs;
parfor i = 1:num_evs
    t_dep_h = EVs_for_baseline(i).t_dep / 60;
    t_in_h = EVs_for_baseline(i).t_in / 60;

    EVs_for_baseline(i).P_base_sequence = EVbaseP_ChargeUntilFull(...
        EVs_for_baseline(i).C, EVs_for_baseline(i).eta,...
        EVs_for_baseline(i).E_tar_original, EVs_for_baseline(i).E_ini,...
        t_dep_h, t_in_h, dt, ...
        EVs_for_baseline(i).r, EVs_for_baseline(i).P_N, ...
        EVs_for_baseline(i).SOC_original, num_time_points, time_points_absolute);
end
for i = 1:num_evs
    EVs(i).P_base_sequence = EVs_for_baseline(i).P_base_sequence;
end
clear EVs_for_baseline;
fprintf('基线功率序列计算完成。\n');


%% 外层循环（长时间步长）
for long_idx = 1:num_long_steps
    t_long_start_minute = (long_idx - 1) * dt_long_minutes;

    %% 长时间步处理
    [lambda_star] = aggregateEVs(EVs, P_tar(long_idx));
    [~, S_agg_next] = calculateVirtualSOC_agg(EVs, dt_long_minutes);

    %% 内层循环（短时间步长）
    for short_idx = 1:num_short_per_long
        step_idx = (long_idx - 1) * num_short_per_long + short_idx;
        t_current_minute = t_long_start_minute + (short_idx - 1) * dt_minutes;
        current_absolute_hour = time_points_absolute(step_idx);

        temp_m3 = zeros(num_evs, 1);
        temp_S_original = zeros(num_evs, 1);
        temp_S_mod = zeros(num_evs, 1);
        temp_P_current = zeros(num_evs, 1);

        %% 更新EV状态 (并行处理)
        EVs_in_parfor = EVs;

        parfor i = 1:num_evs
            EV = EVs_in_parfor(i);
            EV = updateLockState(EV, t_current_minute);

            EV_for_handle = EV;
            EV_temp_with_handle = generateDemandCurve(EV_for_handle);
            current_P_val = 0;
            if isfield(EV_temp_with_handle, 'demandCurve') && isa(EV_temp_with_handle.demandCurve, 'function_handle')
                current_P_val = EV_temp_with_handle.demandCurve(lambda_star);
            else
                 switch EV.state
                    case 'LockON'
                        current_P_val = EV.P_N;
                    case {'LockOFF', 'OFF', 'ON'}
                        current_P_val = 0;
                    otherwise
                        current_P_val = 0;
                 end
            end
            EV.P_current = current_P_val;

            EV = calculateVirtualSOC_upgrade(EV, t_current_minute, dt_minutes);

            if EV.t_dep > EV.t_in
                 m3_val = (EV.E_tar - EV.E_ini) / (EV.eta * ((EV.t_dep - EV.t_in) / 60));
            else
                 m3_val = 0;
            end
            temp_m3(i) = m3_val;

            temp_S_original(i) = EV.S_original;
            temp_S_mod(i) = EV.S_modified;
            temp_P_current(i) = EV.P_current;

            EVs_in_parfor(i) = EV;
        end % 结束 parfor i

        EVs = EVs_in_parfor; % 更新主 EVs 数组

        %% (修改) 计算聚合潜力
        active_participating_indices = find(arrayfun(@(ev, t_abs_h) ...
            ev.ptcp && (t_abs_h >= (ev.t_in / 60)) && (t_abs_h < (ev.t_dep / 60)), ...
            EVs, repmat(current_absolute_hour, num_evs, 1)));

        agg_DeltaP_plus = 0;
        agg_DeltaP_minus = 0;

        if ~isempty(active_participating_indices)
            active_EVs = EVs(active_participating_indices);

            E_reg_min_agg = sum([active_EVs.E_reg_min]);
            E_reg_max_agg = sum([active_EVs.E_reg_max]);
            E_current_agg = sum([active_EVs.E_actual]); % 使用更新后的 E_actual
            t_dep_agg_h = mean([active_EVs.t_dep]) / 60; % 平均离网时间 (小时)
            p_on_agg = sum([active_EVs.P_N]);           % 总额定功率
            P_base_agg = sum(arrayfun(@(ev) ev.P_base_sequence(step_idx), active_EVs)); % 总基线功率
            eta_agg = mean([active_EVs.eta]);           % 平均效率

            % 调用 calculateEVAdjustmentPotentia_new 计算聚合潜力
            [agg_DeltaP_plus, agg_DeltaP_minus] = calculateEVAdjustmentPotentia_new(...
                E_reg_min_agg, E_reg_max_agg, E_current_agg, ...
                t_dep_agg_h, current_absolute_hour, ...
                p_on_agg, P_base_agg, eta_agg, t_adj); % 使用小时单位
        end

        %% 记录结果
        results.EV_S_original(:, step_idx) = temp_S_original;
        results.EV_S_mod(:, step_idx) = temp_S_mod;
        results.m3 = temp_m3; % 可能只需要最后的值

        results.lambda(step_idx) = lambda_star;
        results.S_agg(step_idx) = S_agg_current;
        results.P_agg(step_idx) = sum(temp_P_current);
        results.P_cu(step_idx) = temp_P_current(10); % 第10辆车的功率

        % 记录聚合潜力
        results.EV_Up(step_idx)   = agg_DeltaP_plus;
        results.EV_Down(step_idx) = agg_DeltaP_minus;

    end % 结束 short_idx

    %% 更新聚合SOC
    S_agg_current = S_agg_next;

end % 结束 long_idx

results.P_tar = repelem(P_tar, num_short_per_long);

%% 结果保存与可视化
outputFileName = 'main_potential_aggregated_newlogic_results_simplified.mat';
fprintf('\n正在保存结果到 %s ...\n', outputFileName);
try
    save(outputFileName, 'results', '-v7.3');
    fprintf('结果保存成功。\n');
catch ME_save
    fprintf('*** 保存结果文件时出错: %s ***\n', ME_save.message);
end

% 可视化聚合潜力
figure;
plot(time_points_absolute, results.EV_Up, 'r-', 'DisplayName', 'Aggregated Up Potential');
hold on;
plot(time_points_absolute, results.EV_Down, 'b-', 'DisplayName', 'Aggregated Down Potential');
xlabel('Time (hours)');
ylabel('Potential (kW)');
title('Aggregated EV Regulatory Potential (Based on Aggregated Model with _new Logic)');
legend;
grid on;

% 可视化单车SOC和功率对比
selected_ev = 10;
S_original_vis = results.EV_S_original(selected_ev, :);
S_mod_vis = results.EV_S_mod(selected_ev, :);
P_current_vis = results.P_cu;

figure('Position', [100 100 1200 600]);
yyaxis left;
plot(time_points_absolute, S_original_vis, 'b-', 'LineWidth', 1.5, 'DisplayName', '原始SOC');
hold on;
plot(time_points_absolute, S_mod_vis, 'r--', 'LineWidth', 1.5, 'DisplayName', '修正SOC');
ylabel('SOC'); ylim([-1.1 1.1]); grid on;
yyaxis right;
plot(time_points_absolute, P_current_vis, 'g-.', 'LineWidth', 1.5, 'DisplayName', ['EV ', num2str(selected_ev), ' 功率']);
ylabel('功率 (kW)');
P_lim_vis = [min(P_current_vis)-1, max(P_current_vis)+1];
ylim(P_lim_vis);
xlabel('时间 (小时)'); title(['第', num2str(selected_ev), '台EV的SOC与当前功率对比']);
legend('Location', 'best'); set(gca, 'FontSize', 12); hold off;