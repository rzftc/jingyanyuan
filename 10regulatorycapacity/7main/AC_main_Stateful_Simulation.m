% AC_main_Stateful_Simulation.m
% 功能：AC集群状态化仿真主程序，包含详细绘图与数据导出

clear; close all; clc;
tic;

%% 1. 系统初始化
rng(2023, 'Threefry');          % 固定随机种子
T_total = 24;                   % 总时长
dt = 5/60;                      % 时间分辨率
time_points = 0:dt:T_total;     
T_steps_total = length(time_points);
steps_per_hour = round(1/dt);
num_hours = floor(T_steps_total / steps_per_hour);
base_price = 30;                % 基础电价

%% 2. 读取与初始化设备参数
acFile = 'AC_template2.xlsx';
fprintf('正在初始化空调参数...\n');
try
    ACs = initializeACsFromExcel(acFile);
catch ME
    error('文件加载失败: %s', ME.message);
end
num_AC = length(ACs);

% 备份原始参数
for i = 1:num_AC
    ACs(i).Tset_original = ACs(i).Tset;
    ACs(i).Tmax_original = ACs(i).Tmax;
    ACs(i).Tmin_original = ACs(i).Tmin;
    if ~isfield(ACs(i), 'p_incentive')
        ACs(i).p_incentive = round(60*rand(), 1);
    end
end
fprintf('加载完成: %d 台空调。\n', num_AC);

%% 3. 激励响应参数
p_min = 15; p_max = 50; 
p_min_prime = 10; p_max_prime = 40; 
T_set_max = 3;
current_p = 25.0;
fprintf('\n== 当前仿真电价: %.1f 元 ==\n', current_p);

%% 4. 预计算模块

% 4.1 单体参数与环境温度设定
fprintf('  4.1 预计算单体参数与环境温度...\n');
temp_ACs = ACs;

% 构造环境温度曲线
max_Tset_all = max([ACs.Tset_original]);
T_ja_min = max_Tset_all + 0.1;
T_ja_peak = max_Tset_all + 6.0;
T_ja_mean = (T_ja_min + T_ja_peak) / 2;
T_ja_amp = (T_ja_peak - T_ja_min) / 2;

base_trend = T_ja_mean + T_ja_amp * cos(2*pi*(time_points - 15)/24);
window_size = 2 * steps_per_hour;
noise_padding = ceil(window_size / 2);
raw_noise = randn(1, T_steps_total + 2 * noise_padding);
smooth_noise = movmean(raw_noise, window_size);
fluctuations = (smooth_noise(noise_padding+1:end-noise_padding) / std(smooth_noise)) * (T_ja_amp * 0.2);
base_T_ja_unified = max(base_trend + fluctuations, T_ja_min);

parfor i = 1:num_AC
    participation = calculateParticipation(current_p, base_price);
    [~, ~, deltaT] = incentiveTempAC(current_p, p_min, p_max, p_min_prime, p_max_prime, T_set_max);
    temp_ACs(i).ptcp = (rand() < participation);

    if temp_ACs(i).ptcp
        temp_ACs(i).Tmax = temp_ACs(i).Tset_original + deltaT;
        temp_ACs(i).Tmin = temp_ACs(i).Tset_original - deltaT;
    end
    temp_ACs(i).T_ja = base_T_ja_unified;

    [alpha, beta, gamma] = calculateACABC_single(...
        temp_ACs(i).R, temp_ACs(i).C, temp_ACs(i).eta,...
        temp_ACs(i).Tmax, temp_ACs(i).Tmin, temp_ACs(i).Tset, dt);
    temp_ACs(i).alpha = alpha;
    temp_ACs(i).beta = beta;
    temp_ACs(i).gamma = gamma;
end
ACs = temp_ACs;

% 4.2 聚合模型参数计算
fprintf('  4.2 计算聚合模型参数 (A, B, C)...\n');
ACs_participating = ACs([ACs.ptcp]);
num_AC_participating = length(ACs_participating);

if num_AC_participating == 0
    error('无空调参与，仿真终止。');
end

AggParams = calculateAggregatedACParams(ACs_participating);
fprintf('  聚合参数: A=%.4f, B=%.4f, C=%.4f\n', AggParams.A, AggParams.B, AggParams.C);

% 提取矩阵用于后续计算
% [修正]: 确保 T_ja_matrix 是 [T_steps, num_AC_p]
T_ja_matrix = cat(1, ACs_participating.T_ja)'; 

% [修正]: 确保这里的 repmat 生成 [T_steps, num_AC_p]
% 注意：[ACs...].Tset_original 得到行向量 [1, N]，repmat(..., T, 1) 得到 [T, N]
Tset_matrix = repmat([ACs_participating.Tset_original], T_steps_total, 1);
R_matrix = repmat([ACs_participating.R], T_steps_total, 1);
eta_matrix = repmat([ACs_participating.eta], T_steps_total, 1);

% 4.3 生成电网指令
fprintf('  4.3 生成电网调节指令...\n');
P_grid_command_series = generate_hourly_regulation_signal(T_steps_total, steps_per_hour, num_hours, num_AC_participating);

%% 5. 状态化时序仿真

CURRENT_SOC_AC = [ACs_participating.SOC]';

% 结果存储
Agg_SOC_History = zeros(T_steps_total, 1);
Individual_SOC_History = zeros(T_steps_total, num_AC_participating);
Agg_P_Command_History = zeros(T_steps_total, 1);
Agg_P_Achieved_History = zeros(T_steps_total, 1);
Agg_P_Potential_Up_History = zeros(T_steps_total, 1);
Agg_P_Potential_Down_History = zeros(T_steps_total, 1);
Individual_Power_History = zeros(T_steps_total, num_AC_participating);
Agg_Model_Potential_Up_History = zeros(T_steps_total, 1);
Agg_Model_Potential_Down_History = zeros(T_steps_total, 1);

fprintf('  5. 开始仿真...\n');

for t_idx = 1:T_steps_total
    % 1. 状态记录
    SOC_agg_t = mean(CURRENT_SOC_AC, 'omitnan');
    Agg_SOC_History(t_idx) = SOC_agg_t;
    Individual_SOC_History(t_idx, :) = CURRENT_SOC_AC';

    % 2. 计算临时基线功率总和
    P_base_instant_agg = 0;
    for k = 1:num_AC_participating
         P_base_k = ACbaseP_single(ACs_participating(k).T_ja(t_idx), ...
                                   ACs_participating(k).Tset, ...
                                   ACs_participating(k).R, ...
                                   ACs_participating(k).eta);
         P_base_instant_agg = P_base_instant_agg + P_base_k;
    end

    % 3. 计算聚合模型理论潜力
    if abs(AggParams.B) > 1e-9
        P_agg_energy_up = (1 - AggParams.A * SOC_agg_t - AggParams.C) / (AggParams.B * dt);
        P_agg_energy_down = (0 - AggParams.A * SOC_agg_t - AggParams.C) / (AggParams.B * dt);
    else
        P_agg_energy_up = 0; P_agg_energy_down = 0;
    end
    P_agg_power_up = P_base_instant_agg;
    P_agg_power_down = -P_base_instant_agg;
    
    Agg_Model_Potential_Up_History(t_idx) = min(P_agg_energy_up, P_agg_power_up);
    Agg_Model_Potential_Down_History(t_idx) = max(P_agg_energy_down, P_agg_power_down);

    % 4. 指令获取与目标预测
    Delta_P_S_command = P_grid_command_series(t_idx);
    Agg_P_Command_History(t_idx) = Delta_P_S_command;
    
    SOC_target_next = AggParams.A * SOC_agg_t + AggParams.B * Delta_P_S_command + AggParams.C;
    SOC_target_next = max(0, min(1, SOC_target_next));

    % 5. 单体响应
    temp_AC_Up_agg = 0;
    temp_AC_Down_agg = 0;
    temp_SOC_next = zeros(num_AC_participating, 1);
    temp_P_achieved = zeros(num_AC_participating, 1);

    parfor i = 1:num_AC_participating
        ac_i = ACs_participating(i);
        soc_curr = CURRENT_SOC_AC(i);

        P_base_i = ACbaseP_single(ac_i.T_ja(t_idx), ac_i.Tset, ac_i.R, ac_i.eta);
        [P_plus, P_minus] = calculateACAdjustmentPotentia(...
            P_base_i, 2*abs(P_base_i), 0, ...
            ac_i.alpha, ac_i.beta, ac_i.gamma, soc_curr, dt);

        temp_AC_Up_agg = temp_AC_Up_agg + P_plus;
        temp_AC_Down_agg = temp_AC_Down_agg + P_minus;

        delta_Pj = 0;
        if abs(ac_i.beta) > 1e-9
            delta_Pj = (SOC_target_next - ac_i.alpha * soc_curr - ac_i.gamma) / ac_i.beta;
        end
        delta_Pj_clipped = max(P_minus, min(P_plus, delta_Pj));

        soc_next = updateACSOC_single(soc_curr, delta_Pj_clipped, ...
            ac_i.alpha, ac_i.beta, ac_i.gamma);

        temp_SOC_next(i) = soc_next;
        temp_P_achieved(i) = delta_Pj_clipped;
    end

    % 6. 记录与更新
    Agg_P_Potential_Up_History(t_idx) = temp_AC_Up_agg;
    Agg_P_Potential_Down_History(t_idx) = temp_AC_Down_agg;
    Agg_P_Achieved_History(t_idx) = sum(temp_P_achieved);
    Individual_Power_History(t_idx, :) = temp_P_achieved';
    CURRENT_SOC_AC = temp_SOC_next;
end
fprintf('  仿真完成。\n');

% --- 数据后处理与维度修正 ---

% [修正点]：Tmax_mat 的维度生成
% ACs_participating.Tmax 是标量或，[ACs_participating.Tmax] 构成行向量 [1, N]
% repmat( [1, N], T, 1 ) -> [T, N]，与 Individual_SOC_History 维度一致
Tmax_mat = repmat([ACs_participating.Tmax], T_steps_total, 1);
Tmin_mat = repmat([ACs_participating.Tmin], T_steps_total, 1);

% 计算反推温度
Individual_Temp_History = Tmax_mat - Individual_SOC_History .* (Tmax_mat - Tmin_mat);

P_standby = 0.05;
% Baseline_Power_History 也是 [T, N]
Baseline_Power_History = (T_ja_matrix - Tset_matrix) ./ (R_matrix .* eta_matrix);
Total_Power_History = Baseline_Power_History + Individual_Power_History;
Total_Power_History(Total_Power_History < P_standby) = P_standby;

Agg_Baseline_Power = sum(Baseline_Power_History, 2);
Agg_Total_Power = sum(Total_Power_History, 2);

%% 6. 结果可视化

% 图1: 功率跟踪
figure('Name', '功率跟踪对比', 'Position', [100 100 1000 450]);
plot(time_points, Agg_P_Command_History, 'k:', 'LineWidth', 2.5, 'DisplayName', '电网指令'); hold on;
plot(time_points, Agg_P_Achieved_History, 'r-', 'LineWidth', 1.5, 'DisplayName', '响应功率');
xlabel('时间 (h)'); ylabel('功率 (kW)'); title('图1: 聚合响应功率 vs 电网指令'); legend; grid on; xlim([0, 24]);

% 图2: SOC对比
figure('Name', 'SOC状态对比', 'Position', [100 550 1000 450]);
plot(time_points, Individual_SOC_History, 'LineWidth', 0.5, 'Color', [0.8 0.8 0.8]); hold on;
plot(time_points, Agg_SOC_History, 'k--', 'LineWidth', 2, 'DisplayName', '聚合SOC');
xlabel('时间 (h)'); ylabel('SOC'); title('图2: 单体SOC vs 聚合SOC'); xlim([0, 24]); ylim([-0.1, 1.1]); grid on;

% 图3: 室内温度
figure('Name', '室内温度变化', 'Position', [100 300 1000 450]);
plot(time_points, Individual_Temp_History, 'LineWidth', 0.5);
xlabel('时间 (h)'); ylabel('温度 (°C)'); title('图3: 单体空调室内温度变化'); xlim([0, 24]); grid on;

% 图4: 调节功率
figure('Name', '调节功率', 'Position', [100 400 1000 450]);
plot(time_points, Individual_Power_History, 'LineWidth', 0.5); yline(0, 'k--');
xlabel('时间 (h)'); ylabel('功率 (kW)'); title('图4: 单体调节功率'); xlim([0, 24]); grid on;

% 图5: 总制冷功率
figure('Name', '总制冷功率', 'Position', [100 500 1000 450]);
plot(time_points, Total_Power_History, 'LineWidth', 0.5); yline(P_standby, 'k--');
xlabel('时间 (h)'); ylabel('功率 (kW)'); title('图5: 单体总制冷功率'); xlim([0, 24]); grid on;

% 图6: 聚合功率
figure('Name', '聚合功率对比', 'Position', [100 600 1000 450]);
plot(time_points, Agg_Baseline_Power, 'b--', 'LineWidth', 2, 'DisplayName', '基线功率'); hold on;
plot(time_points, Agg_Total_Power, 'r-', 'LineWidth', 2, 'DisplayName', '总制冷功率');
xlabel('时间 (h)'); ylabel('功率 (kW)'); title('图6: 基线 vs 总功率'); legend; grid on; xlim([0, 24]);

% 图7: 调节潜力对比
figure('Name', '调节潜力对比', 'Position', [100 700 1000 450]);
plot(time_points, Agg_P_Potential_Up_History, 'b-', 'LineWidth', 1.5, 'DisplayName', '单体累加上调'); hold on;
plot(time_points, Agg_P_Potential_Down_History, 'r-', 'LineWidth', 1.5, 'DisplayName', '单体累加下调');
plot(time_points, Agg_Model_Potential_Up_History, 'b--', 'LineWidth', 2, 'DisplayName', '聚合模型上调');
plot(time_points, Agg_Model_Potential_Down_History, 'r--', 'LineWidth', 2, 'DisplayName', '聚合模型下调');
xlabel('时间 (h)'); ylabel('潜力 (kW)'); title('图7: 调节潜力对比'); legend; grid on; xlim([0, 24]);

fprintf('所有任务完成。\n');
toc;