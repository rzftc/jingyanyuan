clear; close all; clc;

tic; % 启动一个总计时器
%% 1. 系统初始化
rng(2023, 'Threefry'); % 固定随机种子
T_total = 24; % 总时长（小时）
dt = 5/60;    % 时间分辨率（小时）
time_points = 0:dt:T_total; % 仿真时间点
T_steps_total = length(time_points);
steps_per_hour = round(1/dt);
num_hours = floor(T_steps_total / steps_per_hour);

base_price = 30; % 基础电价（元/kWh）

%% 2. 初始化 AC 参数
acFile = 'AC_template2.xlsx'; % 确保 AC_template2.xlsx 在路径中
fprintf('正在初始化空调参数...\n');
try
    % 依赖 1initialize/initializeACsFromExcel.m
    ACs = initializeACsFromExcel(acFile);
catch ME
    error('无法加载 %s。请确保 initializeACsFromExcel.m 在路径中。\n错误: %s', acFile, ME.message);
end
num_AC = length(ACs);

% [新增] 检查并处理 P_max_cold 参数
fprintf('正在配置空调最大功率限制 (P_max_cold)...\n');
for i = 1:num_AC
    ACs(i).Tset_original = ACs(i).Tset;
    ACs(i).Tmax_original = ACs(i).Tmax;
    ACs(i).Tmin_original = ACs(i).Tmin;
    if ~isfield(ACs(i), 'p_incentive')
        ACs(i).p_incentive = round(60*rand(), 1);
    end
    
    % [修改关键点]：确保每台空调都有确定的 P_max_cold
    if isfield(ACs(i), 'P_max_cold') && ~isnan(ACs(i).P_max_cold)
        % 使用 Excel 中读取的值 (单位应为 kW)
        ACs(i).P_max_limit = ACs(i).P_max_cold;
    else
        % 兼容性回退：如果没有 P_max_cold 列，使用旧逻辑 (1.1倍最大基线)
        % 这里的 T_ja 还没生成，暂时无法精确计算，先给个标记或警告
        warning('AC %d 缺少 P_max_cold 参数，将在后续步骤动态估算。建议更新 Excel 模板。', i);
        ACs(i).P_max_limit = -1; % 标记为未定义
    end
end
fprintf('加载了 %d 台空调。\n', num_AC);

%% 3. 激励响应参数 (选择一个场景)
p_min = 15; p_max = 50; p_min_prime = 10; p_max_prime = 40; T_set_max = 3;

% --- 选择一个代表性价格进行仿真 ---
current_p = 25.0;
fprintf('\n== 仿真价格场景 (Price: %.1f 元) ==\n', current_p);

%% 4. 预计算 (流程图 步骤1 & 2)

% 4.1 预计算 (alpha, beta, gamma) 和 T_ja
fprintf('  Step 4.1: 预计算单体参数...\n');
temp_ACs = ACs;

% --- [V11 用户请求修改] ---
% 1. 找到所有空调中的最高Tset
max_Tset_all = max([ACs.Tset_original]);
fprintf('  [V11 修正] 检测到最高 Tset 为: %.2f C\n', max_Tset_all);

% 2. 定义 T_ja 曲线参数，确保 T_ja > max_Tset_all
T_ja_min_ambient = max_Tset_all + 0.1; % 保证夜间 T_ja 至少比 Tset 高 0.1 度
T_ja_peak_ambient = max_Tset_all + 2.0; % 白天峰值高 6 度

T_ja_mean = (T_ja_min_ambient + T_ja_peak_ambient) / 2;
T_ja_amplitude = (T_ja_peak_ambient - T_ja_min_ambient) / 2;

% 3. 生成基准趋势 (Trend)，峰值在中午 15：00
base_trend = T_ja_mean + T_ja_amplitude * cos(2*pi*(time_points - 15)/24); 

% 4. 生成波动 (Fluctuations)
window_size = 2 * steps_per_hour; 
noise_padding = ceil(window_size / 2);
white_noise = randn(1, T_steps_total + 2 * noise_padding);
fluctuations_raw = movmean(white_noise, window_size);
fluctuations_centered = fluctuations_raw(noise_padding + 1 : noise_padding + T_steps_total);

% 5. 缩放波动
fluctuation_scale = T_ja_amplitude * 0.2; 
scaled_fluctuations = (fluctuations_centered / std(fluctuations_centered, 'omitnan')) * fluctuation_scale;

% 6. 合并并强制执行约束
base_ambient_temp_unified = base_trend + scaled_fluctuations;
base_ambient_temp_unified = max(base_ambient_temp_unified, T_ja_min_ambient); 

fprintf('  [V11 修正] 生成带波动的 T_ja 曲线: 峰值 approx %.2f C\n', max(base_ambient_temp_unified));
% --- [V11 用户请求修改] 结束 ---


parfor i = 1:num_AC
    % 依赖 5userUncertainties/calculateParticipation.m
    participation = calculateParticipation(current_p, base_price);
    % 依赖 2AC/incentiveTempAC.m
    [~, ~, deltaT_flex_magnitude] = incentiveTempAC(...
        current_p, p_min, p_max, p_min_prime, p_max_prime, T_set_max);
    temp_ACs(i).ptcp = (rand() < participation);

    if temp_ACs(i).ptcp
        temp_ACs(i).Tmax = temp_ACs(i).Tset_original + deltaT_flex_magnitude;
        temp_ACs(i).Tmin = temp_ACs(i).Tset_original - deltaT_flex_magnitude;
    end

    % --- [V11 修正] ---
    % 4. 将统一的 T_ja 分配给所有空调
    temp_ACs(i).T_ja = base_ambient_temp_unified;
    % --- [V11 修正] 结束 ---
    
    % [新增] 如果 Excel 中没填 P_max_cold，在这里根据 T_ja 补救计算
    if temp_ACs(i).P_max_limit < 0
        % 计算全天最大可能基线
        P_base_profile = ACbaseP_single(base_ambient_temp_unified, temp_ACs(i).Tset_original, temp_ACs(i).R, temp_ACs(i).eta);
        temp_ACs(i).P_max_limit = 1.1 * max(abs(P_base_profile));
    end

    % 依赖 2AC/calculateACABC_single.m
    [alpha, beta, gamma] = calculateACABC_single(...
        temp_ACs(i).R, temp_ACs(i).C, temp_ACs(i).eta,...
        temp_ACs(i).Tmax, temp_ACs(i).Tmin, temp_ACs(i).Tset, dt);
    temp_ACs(i).alpha = alpha;
    temp_ACs(i).beta = beta;
    temp_ACs(i).gamma = gamma;
end
ACs = temp_ACs;
fprintf('  Step 4.1: 完成。\n');

% 4.2 计算聚合参数 (流程图 步骤2: "面向SOC状态一致的聚合模型构建")
fprintf('  Step 4.2: 计算聚合模型参数 (A, B, C)...\n');
ACs_participating = ACs([ACs.ptcp]);
num_AC_participating = length(ACs_participating);

if num_AC_participating == 0
    error('该价格下无空调参与，仿真停止。\n');
end

% 依赖 2AC/calculateAggregatedACParams.m
AggParams = calculateAggregatedACParams(ACs_participating);
fprintf('  聚合参数: A=%.4f, B=%.4f, C=%.4f (共 %d 台空调)\n', ...
    AggParams.A, AggParams.B, AggParams.C, num_AC_participating);

% --- [V6] 提取用于计算P_base的参数 ---
T_ja_participating_T = cat(1, ACs_participating.T_ja)';     % (T_steps x N_p)
Tset_vec_p = [ACs_participating.Tset_original]';            % (N_p x 1)
R_vec_p = [ACs_participating.R]';                          % (N_p x 1)
eta_vec_p = [ACs_participating.eta]';                      % (N_p x 1)
P_max_limit_vec = [ACs_participating.P_max_limit]';        % (N_p x 1) [新增]

% 复制为 (T_steps x N_p)
Tset_matrix_p = repmat(Tset_vec_p', T_steps_total, 1);
R_matrix_p = repmat(R_vec_p', T_steps_total, 1);
eta_matrix_p = repmat(eta_vec_p', T_steps_total, 1);
P_max_limit_matrix = repmat(P_max_limit_vec', T_steps_total, 1); % [新增]
% --- [V6] 结束 ---

% 4.3 生成电网指令 (流程图 步骤3)
fprintf('  Step 4.3: 生成电网指令...\n');
% 依赖 1initialize/generate_hourly_regulation_signal.m
P_grid_command_series = generate_hourly_regulation_signal(T_steps_total, steps_per_hour, num_hours, num_AC_participating);


%% 5. 主时间循环 (状态化仿真)

% 5.1 初始化状态向量 (t=0)
CURRENT_SOC_AC = [ACs_participating.SOC]'; % (N_participating x 1)

% 5.2 结果存储 (用于绘图)
Agg_SOC_History = zeros(T_steps_total, 1);
Individual_SOC_History = zeros(T_steps_total, num_AC_participating);
Agg_P_Command_History = zeros(T_steps_total, 1);
Agg_P_Achieved_History = zeros(T_steps_total, 1);
Agg_P_Potential_Up_History = zeros(T_steps_total, 1);   
Agg_P_Potential_Down_History = zeros(T_steps_total, 1); 
Individual_Power_History = zeros(T_steps_total, num_AC_participating); 

% --- [新增 V12] 聚合模型理论潜力存储 ---
Agg_Model_Potential_Up_History = zeros(T_steps_total, 1);   
Agg_Model_Potential_Down_History = zeros(T_steps_total, 1);
% --- [新增 V12] 结束 ---

% --- [新增] 个体潜力存储 (用于保存到 MAT) ---
AC_Up_Individual = zeros(num_AC_participating, T_steps_total);
AC_Down_Individual = zeros(num_AC_participating, T_steps_total);
% --- [新增] 结束 ---

% [修改] 计算聚合物理最大功率 (使用 P_max_cold)
Agg_P_Physical_Max_Total = 0;
for i = 1:num_AC_participating
    % 直接使用 P_max_limit (即 P_max_cold)
    Agg_P_Physical_Max_Total = Agg_P_Physical_Max_Total + ACs_participating(i).P_max_limit;
end
fprintf('  聚合体物理总容量 (Sum of P_max_cold): %.2f kW\n', Agg_P_Physical_Max_Total);

fprintf('  Step 5: 开始 %d 步的状态化仿真...\n', T_steps_total);

for t_idx = 1:T_steps_total

    % 1. 获取当前聚合状态 SOC(t)
    SOC_agg_t = mean(CURRENT_SOC_AC, 'omitnan');
    Agg_SOC_History(t_idx) = SOC_agg_t;
    Individual_SOC_History(t_idx, :) = CURRENT_SOC_AC'; 
    
    % 2. 获取电网指令 ΔP_S (流程图 步骤4)
    Delta_P_S_command = P_grid_command_series(t_idx);
    Agg_P_Command_History(t_idx) = Delta_P_S_command;

    % 3. 预测目标聚合 SOC(t+1) (流程图 步骤5)
    SOC_target_next = AggParams.A * SOC_agg_t + AggParams.B * Delta_P_S_command + AggParams.C;
    SOC_target_next = max(0, min(1, SOC_target_next)); % 约束

    % 4. 临时变量 (用于 parfor)
    temp_AC_Up_agg = 0;
    temp_AC_Down_agg = 0;
    temp_P_base_agg = 0; % 用于聚合模型约束的基线功率总和
    temp_SOC_for_next_step = zeros(num_AC_participating, 1);
    temp_P_achieved_this_step = zeros(num_AC_participating, 1);
    
    temp_AC_Up_Ind = zeros(num_AC_participating, 1);
    temp_AC_Down_Ind = zeros(num_AC_participating, 1);

    % 5. 【核心】状态转移 (parfor)
    parfor i = 1:num_AC_participating

        ac_i = ACs_participating(i);
        soc_current_i = CURRENT_SOC_AC(i); 
        
        % [修改] 重新计算 P_base_i，并实施 P_max_cold 限制
        % 原始热力学基线
        P_base_raw = ACbaseP_single(ac_i.T_ja(t_idx), ac_i.Tset, ac_i.R, ac_i.eta);
        % 强制基线不超过 P_max_cold
        P_base_i = min(P_base_raw, ac_i.P_max_limit);
        
        temp_P_base_agg = temp_P_base_agg + P_base_i;

        % A. 计算当前单体物理潜力 (用于记录和约束)
        % [修改] 传入 ac_i.P_max_limit (P_max_cold) 作为物理上限
        [P_plus, P_minus] = calculateACAdjustmentPotentia(...
            P_base_i, ac_i.P_max_limit, 0, ... 
            ac_i.alpha, ac_i.beta, ac_i.gamma,...
            soc_current_i, dt);

        temp_AC_Up_agg = temp_AC_Up_agg + P_plus;
        temp_AC_Down_agg = temp_AC_Down_agg + P_minus;
        
        temp_AC_Up_Ind(i) = P_plus;
        temp_AC_Down_Ind(i) = P_minus;

        % B. 反解理论功率 ΔP_j (流程图 步骤6)
        delta_Pj_theory = 0;
        if abs(ac_i.beta) > 1e-9
            delta_Pj_theory = (SOC_target_next - ac_i.alpha * soc_current_i - ac_i.gamma) / ac_i.beta;
        end

        % C. 裁剪指令至物理可行域
        delta_Pj_clipped = max(P_minus, min(P_plus, delta_Pj_theory));

        % D. 更新状态 (实现 式 2-10)
        soc_next_i = updateACSOC_single(soc_current_i, delta_Pj_clipped, ...
            ac_i.alpha, ac_i.beta, ac_i.gamma);

        temp_SOC_for_next_step(i) = soc_next_i;
        temp_P_achieved_this_step(i) = delta_Pj_clipped; 
    end
    
    % 6. 存储当前时间步 t 的单体累加潜力
    Agg_P_Potential_Up_History(t_idx) = temp_AC_Up_agg;
    Agg_P_Potential_Down_History(t_idx) = temp_AC_Down_agg;
    
    AC_Up_Individual(:, t_idx) = temp_AC_Up_Ind;
    AC_Down_Individual(:, t_idx) = temp_AC_Down_Ind;
    
    % --- [修正 V12] 计算聚合模型理论潜力 (与单体逻辑对齐) ---
    % 1. 能量约束
    if abs(AggParams.B) > 1e-9
        P_agg_energy_up = (1 - AggParams.A * SOC_agg_t - AggParams.C) / (AggParams.B * dt);
        P_agg_energy_down = (0 - AggParams.A * SOC_agg_t - AggParams.C) / (AggParams.B * dt);
    else
        P_agg_energy_up = 0; P_agg_energy_down = 0;
    end
    
    % 2. 物理功率约束: 
    % [关键修改]：上调约束 = 总 P_max_cold - 当前总(受限)基线功率
    P_agg_power_up = Agg_P_Physical_Max_Total - temp_P_base_agg; 
    
    % 下调约束 = -当前总(受限)基线功率
    P_agg_power_down = -temp_P_base_agg;
    
    % 3. 取交集 (Min/Max)
    Agg_Model_Potential_Up_History(t_idx) = min(P_agg_energy_up, P_agg_power_up);
    Agg_Model_Potential_Down_History(t_idx) = max(P_agg_energy_down, P_agg_power_down);
    
    % 7. 存储实际响应功率
    Agg_P_Achieved_History(t_idx) = sum(temp_P_achieved_this_step);

    % [V5] 存储单体功率历史
    Individual_Power_History(t_idx, :) = temp_P_achieved_this_step';

    % 8. 更新状态向量用于下一循环
    CURRENT_SOC_AC = temp_SOC_for_next_step; 

end % 结束 t_idx 循环

fprintf('  Step 5: 仿真完成。\n');

% --- 步骤 5.5: 反推室内温度 ---
fprintf('  Step 5.5: 正在反推室内温度历史...\n');
Tmax_vec_p = [ACs_participating.Tmax]';
Tmin_vec_p = [ACs_participating.Tmin]';
TRange_vec_p = Tmax_vec_p - Tmin_vec_p;
TRange_vec_p(abs(TRange_vec_p) < 1e-6) = 1e-6;

Tmax_matrix_p = repmat(Tmax_vec_p', T_steps_total, 1);
Tmin_matrix_p = repmat(Tmin_vec_p', T_steps_total, 1);
TRange_matrix_p = repmat(TRange_vec_p', T_steps_total, 1);

Individual_Temp_History = Tmax_matrix_p - Individual_SOC_History .* TRange_matrix_p;
fprintf('  Step 5.5: 温度反推完成。\n');
% --- 结束 ---

% --- [V6] 步骤 5.6: 计算总功率 ---
fprintf('  Step 5.6: 正在计算总制冷功率 (应用 P_max_cold 约束)...\n');

P_standby = 0.05; % 假设 50W 最小功率

% 1. 计算热力学基线 (原始值)
Baseline_Power_History_Raw = (T_ja_participating_T - Tset_matrix_p) ./ (R_matrix_p .* eta_matrix_p);

% [修改] 2. 对基线功率应用 P_max_cold 限制
Baseline_Power_History = min(Baseline_Power_History_Raw, P_max_limit_matrix);

% 3. 计算调度后的总功率 (基线 + 调节)
Total_Power_History = Baseline_Power_History + Individual_Power_History;

% [修改] 4. 对总功率应用 P_max_cold 限制 (双重保险)
Total_Power_History = min(Total_Power_History, P_max_limit_matrix);

% 5. 施加最小功率约束 (待机功率)
Total_Power_History(Total_Power_History < P_standby) = P_standby;

fprintf('  Step 5.6: 总功率计算完成。\n');
% --- [V6] 结束 ---

% --- [新增 V11] 步骤 5.7: 聚合图6所需的数据 ---
fprintf('  Step 5.7: 聚合基线功率和总制冷功率...\n');
% 聚合基线功率 (T_steps x 1)
Agg_Baseline_Power = sum(Baseline_Power_History, 2);
% 聚合总制冷功率 (T_steps x 1)
Agg_Total_Power = sum(Total_Power_History, 2);
fprintf('  Step 5.7: 聚合完成。\n');
% --- [新增 V11] 结束 ---

% --- [Step 5.8 修改版]：保存完整仿真数据到 MAT 文件 ---
fprintf('  Step 5.8: 正在保存完整仿真数据到 MAT 文件...\n');

results = struct();
results.dt = dt;
results.time_points = time_points; 

% 1. 聚合潜力数据
results.Agg_P_Potential_Up_History = Agg_P_Potential_Up_History;     
results.Agg_P_Potential_Down_History = Agg_P_Potential_Down_History; 
results.Agg_Model_Potential_Up_History = Agg_Model_Potential_Up_History;     
results.Agg_Model_Potential_Down_History = Agg_Model_Potential_Down_History; 

% 2. 功率跟踪数据
results.Agg_P_Command_History = Agg_P_Command_History;   
results.Agg_P_Achieved_History = Agg_P_Achieved_History; 

% 3. 功率分析数据
results.Agg_Baseline_Power = Agg_Baseline_Power; 
results.Agg_Total_Power = Agg_Total_Power;       
results.Total_Power_History = Total_Power_History; 
results.Individual_Power_History = Individual_Power_History; 

% 4. 状态与温度数据
results.Individual_SOC_History_Transposed = Individual_SOC_History'; 
results.Individual_Temp_History_Transposed = Individual_Temp_History';

% 5. 其他个体数据
results.AC_Up_Individual = AC_Up_Individual;
results.AC_Down_Individual = AC_Down_Individual;

% 保存文件
output_mat_name = 'AC_Stateful_Simulation_Results_5min_limit.mat';
save(output_mat_name, 'results', '-v7.3');
fprintf('  完整数据已保存至: %s\n', output_mat_name);