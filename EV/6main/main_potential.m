%% main_parfor_incentive.m (在此基础上添加EV潜力计算)
clc; clear; close all;
rng(2024);

%% ===================== 初始化参数 =====================
excelFile = 'resi_1000.xlsx'; % 您的 EV 数据文件名
if ~exist(excelFile, 'file')
    % --- 您原有的文件生成逻辑 ---
    generateEVParameters_real(excelFile, 1000, 1.0); % 示例：生成1000辆居民区EV
    fprintf('已生成参数模板: %s\n', excelFile);
end
[EVs, t_sim, dt_short, dt_long, P_tar] = initializeFromExcel(excelFile);
fprintf('成功加载%d辆EV数据\n', length(EVs));

%%% --- (新增/修改) 时间参数定义 (确保与潜力计算函数匹配) --- %%%
simulation_start_hour = 6; % 假设您的仿真从早上6点开始 (对应 t=360 分钟)
simulation_end_hour   = 30; % 假设您的仿真到第二天早上6点结束 (对应 t=1800 分钟)
dt = dt_short / 60;       % 使用短时间步长 dt_short，转换为小时
t_adj = dt_long / 60;     % 调节时长，这里假设等于长时间步长 dt_long，转换为小时

% 生成绝对时间轴 (小时)
time_points_absolute = simulation_start_hour : dt : (simulation_start_hour + t_sim/60 - dt); % 确保时间点数量正确
num_time_points = length(time_points_absolute);
fprintf('仿真时间范围: %.2f 小时 到 %.2f 小时, 步长: %.3f 小时 (%d 点)\n', ...
    time_points_absolute(1), time_points_absolute(end), dt, num_time_points);
%%% --- (新增/修改 结束) --- %%%


%% ===================== 索引预处理 =====================
num_long_steps = t_sim / dt_long;
num_short_per_long = dt_long / dt_short;
total_steps = num_long_steps * num_short_per_long; % total_steps 现在等于 num_time_points
assert(mod(num_long_steps, 24*(60/dt_long)) == 0, '长时间步数必须是24小时对应步数的整数倍'); % 调整断言

%% ===================== 创建结果存储结构 =====================
num_evs = length(EVs);
results = struct(... % 保留您原有的结果字段
    'P_agg',         zeros(1, total_steps), ...
    'P_base',        zeros(1, total_steps), ... % 这个P_base似乎是聚合后的？
    'S_agg',         zeros(1, total_steps), ...
    'lambda',        zeros(1, total_steps), ...
    'EV_S_original', zeros(num_evs, total_steps), ...
    'EV_S_mod',      zeros(num_evs, total_steps), ...
    'm3',            zeros(num_evs, 1), ...
    'P_tar',         zeros(1, total_steps), ... % 这个P_tar似乎是每短步的？
    'P_cu',          zeros(1, total_steps) ...
    );

%%% --- (新增) 预分配EV调节潜力存储 --- %%%
results.EV_Up   = zeros(1, total_steps); % 聚合上调潜力
results.EV_Down = zeros(1, total_steps); % 聚合下调潜力
% (可选) 个体潜力
% results.EV_Up_Individual = zeros(num_evs, total_steps);
% results.EV_Down_Individual = zeros(num_evs, total_steps);
%%% --- (新增 结束) --- %%%

%% ===================== 初始化EV目标充电功率 & 状态 (向量化) =====================
% --- 您原有的向量化初始化代码 ---
p_real_vec = [EVs.p_real]';
P_h_max_vec = [EVs.P_h_max]';
P_0_vec = [EVs.P_0]';
P_l_min_vec = [EVs.P_l_min]';
Delta_E_h_max_vec = [EVs.Delta_E_h_max]';
Delta_E_q_max_vec = [EVs.Delta_E_q_max]';
E_tar_set_vec = [EVs.E_tar_set]';
E_ini_vec = [EVs.E_ini]';
P_N_vec = [EVs.P_N]';
C_vec = [EVs.C]'; % 需要电池容量

% --- (新增/修改) 添加 E_tar_original, SOC_original, ptcp, E_reg_min/max 等字段 ---
E_tar_original_vec = E_tar_set_vec; % 假设 E_tar_set 就是原始目标
SOC_original_vec = zeros(num_evs, 1); % 假设初始虚拟SOC为0 (initializeFromExcel似乎没给)
E_reg_min_vec = E_tar_original_vec; % 初始化灵活性窗口
E_reg_max_vec = E_tar_original_vec;

% --- (新增/修改) 激励响应参数 和 参与度计算 ---
p_min = 15; p_max = 50;
p_min_prime = 10; p_max_prime = 40;
base_Price_vec = 30 * ones(num_evs, 1); % 假设基准电价为30
p_incentive_vec = [EVs.p_incentive]'; % 从EVs结构体获取激励价格
participation_probabilities_vec = calculateParticipation(p_incentive_vec, base_Price_vec);
ptcp_vec = (rand(num_evs, 1) < participation_probabilities_vec);

% --- (新增/修改) 计算灵活性窗口 E_reg_min/max ---
E_tar_max_flex_vec = 0.2 * C_vec; % 最大灵活性为容量的20%
[deltaE_up_vec, deltaE_down_vec] = incentiveTempEV_updown(p_incentive_vec, p_min, p_max, p_min_prime, p_max_prime, E_tar_max_flex_vec);

participating_indices = find(ptcp_vec);
if ~isempty(participating_indices)
    E_reg_min_vec(participating_indices) = E_tar_original_vec(participating_indices) - deltaE_down_vec(participating_indices);
    E_reg_min_vec(participating_indices) = max(E_reg_min_vec(participating_indices), E_ini_vec(participating_indices)); % 不低于初始电量

    E_reg_max_vec(participating_indices) = E_tar_original_vec(participating_indices) + deltaE_up_vec(participating_indices);
    E_reg_max_vec(participating_indices) = min(E_reg_max_vec(participating_indices), C_vec(participating_indices)); % 不超过电池容量
end

% --- 将计算出的向量写回 EVs 结构体 ---
E_tar_orig_cell = num2cell(E_tar_original_vec);
SOC_orig_cell = num2cell(SOC_original_vec);
ptcp_cell = num2cell(ptcp_vec);
E_reg_min_cell = num2cell(E_reg_min_vec);
E_reg_max_cell = num2cell(E_reg_max_vec);
[EVs.E_tar_original] = E_tar_orig_cell{:};
[EVs.SOC_original] = SOC_orig_cell{:}; % 确保 initializeEVsFromExcel 也创建了这个字段
[EVs.ptcp] = ptcp_cell{:};
[EVs.E_reg_min] = E_reg_min_cell{:};
[EVs.E_reg_max] = E_reg_max_cell{:};
% --- (新增/修改 结束) ---

% --- 您原有的 E_tar, tau_rem, E_exp, E_actual 初始化代码 ---
% 注意：这里的 E_tar 可能需要根据 ptcp_vec 和 deltaE_vec 调整，
% 但您提供的 main_parfor_incentive 代码似乎没有这样做，而是直接用了 incentiveTempEV 的第三个输出来调整 E_tar。
% 为了保持原有逻辑，我暂时注释掉 deltaE_vec 的应用。
% 如果需要严格按照新逻辑，需要修改下面的 E_tar_vec 计算。
% E_tar_vec = max(E_tar_set_vec - deltaE_vec, E_ini_vec); % 这是原代码逻辑
% ... (后面 tau_rem, E_exp, E_actual 的初始化不变) ...
E_tar_vec = max(E_tar_set_vec, E_ini_vec); % 保持原代码逻辑

t_ch_vec = zeros(num_evs, 1);
valid_PN_mask = P_N_vec > 0;
if any(valid_PN_mask)
    t_ch_vec(valid_PN_mask) = 60 .* (E_tar_vec(valid_PN_mask) - E_ini_vec(valid_PN_mask)) ./ P_N_vec(valid_PN_mask);
end
t_ch_vec = max(t_ch_vec, 0);

E_tar_cell = num2cell(E_tar_vec);
[EVs.E_tar] = E_tar_cell{:}; % 这个 E_tar 是经过激励调整的（根据原脚本逻辑）
tau_rem_cell = num2cell(t_ch_vec);
[EVs.tau_rem] = tau_rem_cell{:};
E_ini_cell = num2cell(E_ini_vec);
[EVs.E_exp] = E_ini_cell{:}; % E_exp 应该初始化为 E_ini
[EVs.E_actual] = E_ini_cell{:}; % E_actual 初始化为 E_ini
% --- (新增) 初始化 E_current 和 P_current ---
[EVs.E_current] = E_ini_cell{:};
P_current_init_cell = num2cell(zeros(num_evs, 1));
[EVs.P_current] = P_current_init_cell{:};
% --- (新增 结束) ---


%% ===================== 初始化聚合SOC 和 基线功率 =====================
S_agg_current = mean([EVs.S_original]); % 应该使用 SOC_original? 您的原代码似乎用了 S_original (虚拟SOC)
% [~, P_base_opt] = EVbaseP_aggregate_short(EVs, S_agg_current, 48, dt_long); % 聚合基线功率计算

%%% --- (修改/新增) 计算个体基线功率序列 P_base_sequence --- %%%
fprintf('正在预计算所有EV的基线功率序列...\n');
EVs_for_baseline = EVs; % 复制一份用于计算
parfor i = 1:num_evs
    % 注意：EVbaseP_ChargeUntilFull 需要 C_EV, eta, E_tar_original, E_in, t_dep(小时), t_in(小时), dt(小时), r, p_on, SOC_original, num_time_points, time_points_absolute(小时)
    t_dep_h = EVs_for_baseline(i).t_dep / 60; % 分钟转小时
    t_in_h = EVs_for_baseline(i).t_in / 60;   % 分钟转小时

    EVs_for_baseline(i).P_base_sequence = EVbaseP_ChargeUntilFull(...
        EVs_for_baseline(i).C, EVs_for_baseline(i).eta,...
        EVs_for_baseline(i).E_tar_original, EVs_for_baseline(i).E_ini,...
        t_dep_h, t_in_h, dt, ... % 使用小时单位的时间
        EVs_for_baseline(i).r, EVs_for_baseline(i).P_N, ... % p_on 使用 P_N
        EVs_for_baseline(i).SOC_original, num_time_points, time_points_absolute); % 使用绝对时间轴
end
% 将计算好的基线序列赋回给主 EVs 结构体
for i = 1:num_evs
    EVs(i).P_base_sequence = EVs_for_baseline(i).P_base_sequence;
end
clear EVs_for_baseline; % 释放内存
fprintf('基线功率序列计算完成。\n');
%%% --- (修改/新增 结束) --- %%%


%% ===================== 外层循环（长时间步长） =====================
for long_idx = 1:num_long_steps
    t_long_start_minute = (long_idx - 1) * dt_long; % 长步开始的分钟数 (相对0)

    %% --------------------- 长时间步处理 ---------------------
    % EVs = distributeBasePower(EVs, 1, P_base_opt); % 分配聚合基线功率 (可能不再需要)
    [lambda_star] = aggregateEVs(EVs, P_tar(long_idx)); % 计算 lambda*
    [P_agg_long, S_agg_next] = calculateVirtualSOC_agg(EVs, dt_long); % 计算下一个长步的聚合SOC

    %% ===================== 内层循环（短时间步长） =====================
    for short_idx = 1:num_short_per_long
        step_idx = (long_idx - 1) * num_short_per_long + short_idx;
        t_current_minute = t_long_start_minute + (short_idx - 1) * dt_short; % 当前短步开始的分钟数 (相对0)
        current_absolute_hour = time_points_absolute(step_idx); % 当前绝对时间 (小时)

        %%% --- (新增) 初始化当前短步的聚合潜力 --- %%%
        current_t_total_DeltaP_plus = 0;
        current_t_total_DeltaP_minus = 0;
        %%% --- (新增 结束) --- %%%

        temp_m3 = zeros(num_evs, 1);
        temp_S_original = zeros(num_evs, 1);
        temp_S_mod = zeros(num_evs, 1);
        temp_P_current = zeros(num_evs, 1);

        %% --------------------- 更新EV状态 (并行处理) ---------------------
        % 使用 parfor 加速时需要确保 EVs 结构体在循环内部被正确处理
        % 创建一个临时的 EVs 副本以在 parfor 内部修改
        EVs_in_parfor = EVs;

        parfor i = 1:num_evs
            EV = EVs_in_parfor(i); % 获取当前 EV 副本

            % 注意：updateLockState 需要相对于 0 点的分钟数，并且处理跨天
            % 您的 t_in/t_dep 似乎已经是绝对分钟数，无需取模
            EV = updateLockState(EV, t_current_minute);

            % --- 您原有的 demandCurve 计算逻辑 ---
            EV_for_handle = EV;
            EV_temp_with_handle = generateDemandCurve(EV_for_handle);
            current_P_val = 0;
            if isfield(EV_temp_with_handle, 'demandCurve') && isa(EV_temp_with_handle.demandCurve, 'function_handle')
                current_P_val = EV_temp_with_handle.demandCurve(lambda_star);
            else % 后备逻辑
                 switch EV.state
                    case 'LockON'
                        current_P_val = EV.P_N;
                    case {'LockOFF', 'OFF'}
                        current_P_val = 0;
                    otherwise % 可能的 'ON' 状态或其他未定义状态
                        current_P_val = 0; % 默认为0
                 end
            end
            EV.P_current = current_P_val;
            % --- 结束 demandCurve 计算 ---

            % --- 调用 calculateVirtualSOC_upgrade 更新 S_original 和 S_modified ---
            EV = calculateVirtualSOC_upgrade(EV, t_current_minute, dt_short);

            % --- 计算 m3 (与潜力计算无关，但原脚本有) ---
             if EV.t_dep > EV.t_in % 避免除以零
                 m3_val = (EV.E_tar - EV.E_ini) / (EV.eta * ((EV.t_dep - EV.t_in) / 60));
             else
                 m3_val = 0;
             end
             temp_m3(i) = m3_val; % 存储 m3

             %%% --- (新增) 计算调节潜力 --- %%%
             DeltaP_plus_i = 0;
             DeltaP_minus_i = 0;
             is_online_h = (current_absolute_hour >= (EV.t_in / 60)) && (current_absolute_hour < (EV.t_dep / 60)); % 判断是否在线 (小时)

             if EV.ptcp && is_online_h % 只有参与且在线的 EV 才有潜力
                 P_base_i = EV.P_base_sequence(step_idx); % 获取当前短步的基线功率
                 t_dep_h = EV.t_dep / 60; % 离网时间 (小时)

                 [DeltaP_plus_i, DeltaP_minus_i] = calculateEVAdjustmentPotentia_new(...
                     EV.E_reg_min, EV.E_reg_max, EV.E_actual, ... % 使用 E_actual 作为当前电量近似
                     t_dep_h, current_absolute_hour, ...
                     EV.P_N, P_base_i, EV.eta, t_adj); % 使用小时单位
             end
             %%% --- (新增 结束) --- %%%

            % --- 存储临时结果 ---
            temp_S_original(i) = EV.S_original;
            temp_S_mod(i) = EV.S_modified;
            temp_P_current(i) = EV.P_current;

            %%% --- (新增) 累加个体潜力到聚合值 (使用原子操作或临时变量) --- %%%
            % 注意：直接在 parfor 中累加到外部变量是不安全的
            % 这里我们将个体潜力存储起来，在 parfor 结束后累加
            temp_delta_p_plus(i) = DeltaP_plus_i;
            temp_delta_p_minus(i) = DeltaP_minus_i;
            %%% --- (新增 结束) --- %%%

            % --- 将更新后的 EV 写回临时数组 ---
            EVs_in_parfor(i) = EV;

        end % 结束 parfor i

        % --- 在 parfor 之后更新主 EVs 数组 ---
        EVs = EVs_in_parfor;

        % --- 在 parfor 之后累加聚合潜力 ---
        current_t_total_DeltaP_plus = sum(temp_delta_p_plus);
        current_t_total_DeltaP_minus = sum(temp_delta_p_minus);
        clear temp_delta_p_plus temp_delta_p_minus; % 清理临时变量

        % --- 记录结果 ---
        results.EV_S_original(:, step_idx) = temp_S_original;
        results.EV_S_mod(:, step_idx) = temp_S_mod;
        results.m3 = temp_m3; % m3 似乎只需要最后一个值?

        results.lambda(step_idx) = lambda_star;
        results.S_agg(step_idx) = S_agg_current; % 记录长步开始时的聚合SOC
        results.P_agg(step_idx) = sum(temp_P_current); % 记录当前短步的实际聚合功率

        results.P_cu(step_idx) = temp_P_current(10); % 记录第10辆车的功率

        %%% --- (新增) 记录聚合潜力 --- %%%
        results.EV_Up(step_idx)   = current_t_total_DeltaP_plus;
        results.EV_Down(step_idx) = current_t_total_DeltaP_minus;
        %%% --- (新增 结束) --- %%%

        % --- (可选新增) 记录个体潜力 ---
        % results.EV_Up_Individual(:, step_idx) = temp_delta_p_plus;
        % results.EV_Down_Individual(:, step_idx) = temp_delta_p_minus;

    end % 结束 short_idx

    %% --------------------- 更新聚合SOC ---------------------
    S_agg_current = S_agg_next; % 更新为下一个长步开始时的聚合SOC

end % 结束 long_idx

results.P_tar = repelem(P_tar, num_short_per_long); % 扩展目标功率以匹配短步长

%% ===================== 结果保存与可视化 =====================
% --- 您原有的结果保存代码 ---
outputFileName = 'main_parfor_incentive_results_with_potential.mat'; % 定义新的输出文件名
fprintf('\n正在保存包含EV潜力的结果到 %s ...\n', outputFileName);
try
    % 保存整个 results 结构体，现在它包含了新增的 EV_Up 和 EV_Down
    save(outputFileName, 'results', '-v7.3');
    fprintf('结果保存成功。\n');
catch ME_save
    fprintf('*** 保存结果文件时出错: %s ***\n', ME_save.message);
end

% --- 您原有的可视化代码 ---
% (可以添加绘制 EV_Up 和 EV_Down 曲线的代码)
figure;
plot(time_points_absolute, results.EV_Up, 'r-', 'DisplayName', 'Aggregated Up Potential');
hold on;
plot(time_points_absolute, results.EV_Down, 'b-', 'DisplayName', 'Aggregated Down Potential');
xlabel('Time (hours)');
ylabel('Potential (kW)');
title('Aggregated EV Regulatory Potential');
legend;
grid on;