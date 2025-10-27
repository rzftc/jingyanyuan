%% main_v5_S_modified.m (初始化向量化 - 注释精简)
clc; clear; close all;
rng(2024);

%% ===================== 初始化参数 =====================
% 生成/加载EV参数文件
excelFile = '../0inputdata/residential_all_models.xlsx';
if ~exist(excelFile, 'file')
    generateEVParameters_real(excelFile, 100, 0.6); % 自动生成100辆EV数据
    fprintf('已生成参数模板: %s\n', excelFile);
end

% 从Excel加载参数
[EVs, t_sim, dt_short, dt_long, P_tar] = initializeFromExcel(excelFile);
fprintf('成功加载%d辆EV数据\n', length(EVs));
% plotDeltaEComparison(EVs)
%% ===================== 索引预处理 =====================
num_long_steps = t_sim / dt_long;
num_short_per_long = dt_long / dt_short;
total_steps = num_long_steps * num_short_per_long;
hours_per_step = round(24 / num_long_steps, 2);
assert(mod(num_long_steps, 24) == 0, '长时间步数必须是24的整数倍');
repeat_factor = num_long_steps / 24; % 每小时的步数

%% ===================== 创建结果存储结构 =====================
num_evs = length(EVs); % 在此获取 num_evs
results = struct(...
    'P_agg',        zeros(1, total_steps), ...
    'P_base',       zeros(1, total_steps), ...
    'S_agg',        zeros(1, total_steps), ...
    'lambda',       zeros(1, total_steps), ...
    'EV_S_original', zeros(num_evs, total_steps), ... % 使用 num_evs
    'EV_S_mod',      zeros(num_evs, total_steps), ... % 使用 num_evs
    'm3',            zeros(num_evs, 1), ...
    'P_tar',         zeros(1, total_steps), ...
    'P_cu',          zeros(1, total_steps) ...
);

%% ===================== 初始化EV目标充电功率 (向量化改造) =====================
% --- 提取所需向量 ---
p_real_vec = [EVs.p_real]';
P_h_max_vec = [EVs.P_h_max]';
P_0_vec = [EVs.P_0]';
P_l_min_vec = [EVs.P_l_min]';
Delta_E_h_max_vec = [EVs.Delta_E_h_max]';
Delta_E_q_max_vec = [EVs.Delta_E_q_max]';
E_tar_set_vec = [EVs.E_tar_set]';
E_ini_vec = [EVs.E_ini]';
P_N_vec = [EVs.P_N]';

% --- 调用向量化函数计算 Delta_E ---
Delta_E_vec = calculateDeltaE_vec(p_real_vec, P_h_max_vec, P_0_vec, P_l_min_vec, Delta_E_h_max_vec, Delta_E_q_max_vec);

% --- 向量化计算 E_tar 和 tau_rem ---
E_tar_vec = max(E_tar_set_vec - Delta_E_vec, E_ini_vec); % 电量下限保护

% 计算 t_ch (tau_rem)，注意处理 P_N 可能为0的情况
t_ch_vec = zeros(num_evs, 1);
valid_PN_mask = P_N_vec > 0;
if any(valid_PN_mask)
    t_ch_vec(valid_PN_mask) = 60 .* (E_tar_vec(valid_PN_mask) - E_ini_vec(valid_PN_mask)) ./ P_N_vec(valid_PN_mask);
end
t_ch_vec = max(t_ch_vec, 0); % 确保非负

% --- 使用 cell 和 deal 将向量结果高效写回 EVs 结构体数组 ---
E_tar_cell = num2cell(E_tar_vec);
[EVs.E_tar] = E_tar_cell{:};

tau_rem_cell = num2cell(t_ch_vec);
[EVs.tau_rem] = tau_rem_cell{:};

% --- 初始化 E_exp 和 E_actual (也向量化) ---
E_ini_cell = num2cell(E_ini_vec); % E_ini 已提取
[EVs.E_exp] = E_ini_cell{:};
[EVs.E_actual] = E_ini_cell{:};

%% ===================== 初始化聚合SOC =====================
S_agg_current = mean([EVs.S_original]); % S_original 初始为0
[S_agg_opt, P_base_opt] = EVbaseP_aggregate_short(EVs, S_agg_current, 48, dt_long);

%% ===================== (可选) 预提取常用向量 =====================
% 如果分析器显示循环内访问结构体字段是瓶颈，可以取消注释这些
% eta_vec = [EVs.eta]';
% t_dep_vec = [EVs.t_dep]';
% t_in_vec = [EVs.t_in]';

%% ===================== 外层循环（长时间步长） =====================
for long_idx = 1:num_long_steps
    t_long = (long_idx - 1) * dt_long;
    %% --------------------- 长时间步处理 ---------------------
    EVs = distributeBasePower(EVs, 1,P_base_opt); % t_current 参数实际未使用
    [lambda_star] = aggregateEVs(EVs, P_tar(long_idx));
    [P_agg, S_agg_next] = calculateVirtualSOC_agg(EVs, dt_long); % P_agg 在短步中会被覆盖

    %% ===================== 内层循环（短时间步长） =====================
    for short_idx = 1:num_short_per_long
        step_idx = (long_idx - 1) * num_short_per_long + short_idx;
        t_current = t_long + (short_idx - 1) * dt_short;

        %% --------------------- 更新EV状态 (保持原有调用逻辑) ---------------------
        for i = 1:length(EVs) % 使用 length(EVs) 或 num_evs 均可
            EV = EVs(i);
            % --- 函数调用逻辑保持不变 ---
            EV = updateLockState(EV, t_current);
            EV = generateDemandCurve(EV);
            EV.P_current = EV.demandCurve(lambda_star);
            EV = calculateVirtualSOC_upgrade(EV, t_current, dt_short); % 应调用优化后的版本
            % --- 函数调用逻辑结束 ---

            % 计算 m3 (保持在循环内，因为依赖单 EV 更新后的状态)
            % 注意 P_N, eta, t_dep, t_in 需要是 EV 的字段
            if isfield(EV, 'eta') && isfield(EV, 't_dep') && isfield(EV, 't_in') && (EV.t_dep > EV.t_in)
                 m3 = (EV.E_tar - EV.E_ini) / (EV.eta * ((EV.t_dep - EV.t_in) / 60));
            else
                 m3 = 0;
            end
            results.m3(i) = m3; % 仍然只保留最后的值

            EVs(i) = EV; % 写回更新后的EV
        end

        %% --------------------- 记录结果 ---------------------
        results.lambda(step_idx) = lambda_star;
        results.S_agg(step_idx) = S_agg_current; % 使用长步开始时的 S_agg
        results.P_agg(step_idx) = sum([EVs.P_current]); % 根据本短步最终 P_current 计算

        results.P_tar = repelem(P_tar, num_short_per_long);
        results.P_cu(step_idx) = EVs(10).P_current;

        % 记录所有EV的SOC (保持不变)
        for ev_idx = 1:length(EVs)
            results.EV_S_original(ev_idx, step_idx) = EVs(ev_idx).S_original;
            results.EV_S_mod(ev_idx, step_idx) = EVs(ev_idx).S_modified;
        end
    end

    %% --------------------- 更新聚合SOC ---------------------
    S_agg_current = S_agg_next;
end
%% ===================== 可视化结果 =====================
% (可视化代码保持不变)
selected_ev = 10;
time_min = (0:total_steps-1) * dt_short;
time_h = time_min / 60;

S_original = results.EV_S_original(selected_ev, :);
S_mod = results.EV_S_mod(selected_ev, :);
P_current_vis = results.P_cu; % 使用记录的 P_cu (EV10的功率)

figure('Position', [100 100 1200 600]);

yyaxis left;
plot(time_h, S_original, 'b-', 'LineWidth', 1.5, 'DisplayName', '原始SOC (S_{original})');
hold on;
plot(time_h, S_mod, 'r--', 'LineWidth', 1.5, 'DisplayName', '修正SOC (S_{modified})');
ylabel('SOC');
ylim([-1.1 1.1]);
grid on;

yyaxis right;
plot(time_h, P_current_vis, 'g-.', 'LineWidth', 1.5, 'DisplayName', ['EV ', num2str(selected_ev), ' 功率 (P_{current})']);
ylabel('功率 (kW)');
P_lim = [min(P_current_vis)-1, max(P_current_vis)+1];
ylim(P_lim);

xlabel('时间 (小时)');
title(['第', num2str(selected_ev), '台EV的SOC与当前功率对比']);
legend('Location', 'best');
set(gca, 'FontSize', 12);
hold off;

% 可选调用其他绘图函数
% plotLambdaAndSOC(results, dt_short, selected_ev);
% plotPowerComparison(results, dt_short);
% plotLambdaAndAggSOC(results, dt_short);