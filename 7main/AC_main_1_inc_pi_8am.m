
clear; close all; clc;

tic; % 启动一个总计时器
%% 1. 系统初始化
rng(2023, 'Threefry'); % 固定随机种子
T_total = 24; % 总时长（小时）
dt = 5/60;    % 时间分辨率（小时）

% [修改] 设定起始时间为早上 8:00
T_start_hour = 8; 
% [修改] 时间轴平移：从 8:00 到 32:00 (次日 8:00)
time_points = T_start_hour : dt : (T_start_hour + T_total); 

T_steps_total = length(time_points);
steps_per_hour = round(1/dt);
num_hours = floor(T_steps_total / steps_per_hour);

base_price = 30; % 基础电价（元/kWh）

%% 2. 初始化 AC 参数
acFile = 'actest.xlsx'; % 确保 AC_template1.xlsx 在路径中
fprintf('正在初始化空调参数...\n');
try
    % 依赖 1initialize/initializeACsFromExcel.m
    ACs = initializeACsFromExcel(acFile);
catch ME
    error('无法加载 %s。请确保 initializeACsFromExcel.m 在路径中。\n错误: %s', acFile, ME.message);
end
num_AC = length(ACs);

% 备份原始设置
for i = 1:num_AC
    ACs(i).Tset_original = ACs(i).Tset;
    ACs(i).Tmax_original = ACs(i).Tmax;
    ACs(i).Tmin_original = ACs(i).Tmin;
    ACs(i).SOC = 0.5;
    if ~isfield(ACs(i), 'p_incentive')
        ACs(i).p_incentive = round(60*rand(), 1);
    end
end
fprintf('加载了 %d 台空调。\n', num_AC);

%% 3. 激励响应参数
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
T_ja_peak_ambient = max_Tset_all + 6.0; % 白天峰值高 6 度

T_ja_mean = (T_ja_min_ambient + T_ja_peak_ambient) / 2;
T_ja_amplitude = (T_ja_peak_ambient - T_ja_min_ambient) / 2;

% 3. 生成基准趋势 (Trend)，峰值在中午 15：00
% [说明] time_points 现已变为 8~32，余弦函数会自动计算对应的物理时间温度
% 例如 t=8时 cos((8-15)/24*2pi) 处于升温段；t=15时达到峰值
base_trend = T_ja_mean + T_ja_amplitude * cos(2*pi*(time_points - 15)/24); 

% 4. 生成波动 (Fluctuations)
window_size = 2 * steps_per_hour; 
noise_padding = ceil(window_size / 2);
white_noise = randn(1, T_steps_total + 2 * noise_padding);
fluctuations_raw = movmean(white_noise, window_size);
fluctuations_centered = fluctuations_raw(noise_padding + 1 : noise_padding + T_steps_total);
fluctuation_scale = T_ja_amplitude * 0.2; 
scaled_fluctuations = (fluctuations_centered / std(fluctuations_centered, 'omitnan')) * fluctuation_scale;

% 6. 合并并强制执行约束
base_ambient_temp_unified = base_trend + scaled_fluctuations;
base_ambient_temp_unified = max(base_ambient_temp_unified, T_ja_min_ambient); 

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
    temp_ACs(i).T_ja = base_ambient_temp_unified;
    % --- [V11 修正] 结束 ---

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

% 4.2 计算聚合参数
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

% =========================================================================
% [系统校正] 变频器 PI 参数 - 二阶欠阻尼稳定模式
% =========================================================================
% 目标：产生明显的超调(Up/Down波动)并最终严格相等(消除静差)。
% 原理：
% 1. 降低 Kp：避免系统发散/不稳定 (Kp*beta < 2.0)。
% 2. 提高 Ki：强积分作用会引入相移，导致系统产生"浪涌"般的震荡收敛，
%    并保证最终误差为0。

all_betas = [ACs_participating.beta];
mean_beta = mean(all_betas(all_betas > 1e-6)); 
if isempty(mean_beta) || isnan(mean_beta), mean_beta = 0.02; end

% 参数设定：
% Kp = 0.8/beta: 提供适度的瞬态响应，留有余地给积分项。
% Ki = 0.5/beta: 非常强的积分增益（通常Ki远小于Kp，这里故意调大以制造震荡）。
Kp_inverter = 1.01 / mean_beta;  
Ki_inverter = 0.7 / mean_beta;   

fprintf('\n*** [控制模式] 欠阻尼震荡收敛模式 (Underdamped) ***\n');
fprintf('    平均 Beta: %.5f\n', mean_beta);
fprintf('    设定 Kp: %.2f (稳定基准)\n', Kp_inverter);
fprintf('    设定 Ki: %.2f (强积分->制造波动+消除误差)\n', Ki_inverter);
fprintf('    预期效果: 功率曲线应呈现"冲过头再拉回"的波浪状，并最终重合。\n');
% =========================================================================

% --- [V6] 提取用于计算P_base的参数 ---
T_ja_participating_T = cat(1, ACs_participating.T_ja)';     
Tset_vec_p = [ACs_participating.Tset_original]';            
R_vec_p = [ACs_participating.R]';                          
eta_vec_p = [ACs_participating.eta]';                      
Tset_matrix_p = repmat(Tset_vec_p', T_steps_total, 1);
R_matrix_p = repmat(R_vec_p', T_steps_total, 1);
eta_matrix_p = repmat(eta_vec_p', T_steps_total, 1);
% --- [V6] 结束 ---

% =========================================================================
% 4.3 [修改] 生成用户指定的阶跃指令
% =========================================================================
fprintf('  Step 4.3: 生成用户指定的阶跃电网指令...\n');

P_grid_command_series = zeros(T_steps_total, 1);
steps_4h = 4 * steps_per_hour; 

% [用户指定]
% command_levels = [0, 100, 0, -110, 100, -90]; 
user_commands = [100, -200, 0, 100, 250, 0]; 

for k = 1:length(user_commands)
    start_idx = (k-1)*steps_4h + 1;
    end_idx = min(k*steps_4h, T_steps_total);
    if start_idx > T_steps_total, break; end
    
    P_grid_command_series(start_idx:end_idx) = user_commands(k);
end
% =========================================================================


%% 5. 主时间循环 (状态化仿真)

% 5.1 初始化状态向量 (t=0)
CURRENT_SOC_AC = [ACs_participating.SOC]'; 

% --- [新增] 初始化 PI 控制器积分误差 ---
CURRENT_INTEGRAL_ERROR = zeros(num_AC_participating, 1);

% 5.2 结果存储
Agg_SOC_History = zeros(T_steps_total, 1);
Individual_SOC_History = zeros(T_steps_total, num_AC_participating);
Agg_P_Command_History = zeros(T_steps_total, 1);
Agg_P_Achieved_History = zeros(T_steps_total, 1);
Agg_P_Potential_Up_History = zeros(T_steps_total, 1);   
Agg_P_Potential_Down_History = zeros(T_steps_total, 1); 
Individual_Power_History = zeros(T_steps_total, num_AC_participating); 

Agg_Model_Potential_Up_History = zeros(T_steps_total, 1);   
Agg_Model_Potential_Down_History = zeros(T_steps_total, 1);

AC_Up_Individual = zeros(num_AC_participating, T_steps_total);
AC_Down_Individual = zeros(num_AC_participating, T_steps_total);

% 计算物理最大功率
Agg_P_Physical_Max_Total = 0;
for i = 1:num_AC_participating
    ac_i = ACs_participating(i);
    P_base_profile_i = ACbaseP_single(ac_i.T_ja, ac_i.Tset, ac_i.R, ac_i.eta);
    P_max_i = 1.1*max(abs(P_base_profile_i));
    Agg_P_Physical_Max_Total = Agg_P_Physical_Max_Total + P_max_i;
end

fprintf('  Step 5: 开始 %d 步的状态化仿真 (8:00 - 次日8:00)...\n', T_steps_total);

for t_idx = 1:T_steps_total

    % 1. 获取当前聚合状态 SOC(t)
    SOC_agg_t = mean(CURRENT_SOC_AC, 'omitnan');
    Agg_SOC_History(t_idx) = SOC_agg_t;
    Individual_SOC_History(t_idx, :) = CURRENT_SOC_AC'; 
    
    % 2. 获取电网指令
    Delta_P_S_command = P_grid_command_series(t_idx);
    Agg_P_Command_History(t_idx) = Delta_P_S_command;

    % 3. 预测目标聚合 SOC(t+1)
    SOC_target_next = AggParams.A * SOC_agg_t + AggParams.B * Delta_P_S_command + AggParams.C;
    SOC_target_next = max(0, min(1, SOC_target_next)); 

    % 4. 临时变量
    temp_AC_Up_agg = 0;
    temp_AC_Down_agg = 0;
    temp_P_base_agg = 0; 
    temp_SOC_for_next_step = zeros(num_AC_participating, 1);
    temp_P_achieved_this_step = zeros(num_AC_participating, 1);
    temp_AC_Up_Ind = zeros(num_AC_participating, 1);
    temp_AC_Down_Ind = zeros(num_AC_participating, 1);
    temp_Integral_Error_Next = zeros(num_AC_participating, 1);

    % 5. 【核心】状态转移 (parfor)
    parfor i = 1:num_AC_participating

        ac_i = ACs_participating(i);
        soc_current_i = CURRENT_SOC_AC(i); 
        int_error_i = CURRENT_INTEGRAL_ERROR(i); 
        
        % 重新计算P_base
        P_base_i = ACbaseP_single(ac_i.T_ja(t_idx), ac_i.Tset, ac_i.R, ac_i.eta);
        temp_P_base_agg = temp_P_base_agg + P_base_i;

        % A. 计算当前单体物理潜力 (上/下限)
        [P_plus, P_minus] = calculateACAdjustmentPotentia(...
            P_base_i, 1.1*max(abs(ACbaseP_single(ac_i.T_ja, ac_i.Tset, ac_i.R, ac_i.eta))), 0, ... 
            ac_i.alpha, ac_i.beta, ac_i.gamma,...
            soc_current_i, dt);

        temp_AC_Up_agg = temp_AC_Up_agg + P_plus;
        temp_AC_Down_agg = temp_AC_Down_agg + P_minus;
        
        temp_AC_Up_Ind(i) = P_plus;
        temp_AC_Down_Ind(i) = P_minus;

        % ==========================================================
        % B. [PI 控制核心]
        % ==========================================================
        % 传入 P_plus/P_minus 以便在底层做抗饱和检查
        [delta_Pj_PI_raw, int_error_next] = calculatePIPower_single(...
            soc_current_i, ...
            SOC_target_next, ...
            int_error_i, ...
            Kp_inverter, ...
            Ki_inverter, ...
            dt, ...
            P_plus, ...  
            P_minus);    
        
        temp_Integral_Error_Next(i) = int_error_next;
        % ==========================================================

        % C. 裁剪指令至物理可行域
        % (注意：实际执行的功率 delta_Pj_clipped 必须在 P_minus 和 P_plus 之间)
        delta_Pj_clipped = max(P_minus, min(P_plus, delta_Pj_PI_raw));

        % D. 更新状态
        soc_next_i = updateACSOC_single(soc_current_i, delta_Pj_clipped, ...
            ac_i.alpha, ac_i.beta, ac_i.gamma);

        temp_SOC_for_next_step(i) = soc_next_i;
        temp_P_achieved_this_step(i) = delta_Pj_clipped; 
    end
    
    % 6. 存储聚合结果
    Agg_P_Potential_Up_History(t_idx) = temp_AC_Up_agg;
    Agg_P_Potential_Down_History(t_idx) = temp_AC_Down_agg;
    
    AC_Up_Individual(:, t_idx) = temp_AC_Up_Ind;
    AC_Down_Individual(:, t_idx) = temp_AC_Down_Ind;
    
    if abs(AggParams.B) > 1e-9
        P_agg_energy_up = (1 - AggParams.A * SOC_agg_t - AggParams.C) / (AggParams.B * dt);
        P_agg_energy_down = (0 - AggParams.A * SOC_agg_t - AggParams.C) / (AggParams.B * dt);
    else
        P_agg_energy_up = 0; P_agg_energy_down = 0;
    end
    
    P_agg_power_up = Agg_P_Physical_Max_Total - temp_P_base_agg; 
    P_agg_power_down = -temp_P_base_agg;
    
    Agg_Model_Potential_Up_History(t_idx) = min(P_agg_energy_up, P_agg_power_up);
    Agg_Model_Potential_Down_History(t_idx) = max(P_agg_energy_down, P_agg_power_down);
    
    % 7. 存储实际响应
    Agg_P_Achieved_History(t_idx) = sum(temp_P_achieved_this_step);
    Individual_Power_History(t_idx, :) = temp_P_achieved_this_step';

    % 8. 更新状态
    CURRENT_SOC_AC = temp_SOC_for_next_step; 
    
    % [新增] 更新积分误差
    CURRENT_INTEGRAL_ERROR = temp_Integral_Error_Next;

end % 结束 t_idx 循环

fprintf('  Step 5: 仿真完成。\n');

% --- 后处理代码 ---
Tmax_vec_p = [ACs_participating.Tmax]';
Tmin_vec_p = [ACs_participating.Tmin]';
TRange_vec_p = Tmax_vec_p - Tmin_vec_p;
TRange_vec_p(abs(TRange_vec_p) < 1e-6) = 1e-6;

Tmax_matrix_p = repmat(Tmax_vec_p', T_steps_total, 1);
Tmin_matrix_p = repmat(Tmin_vec_p', T_steps_total, 1);
TRange_matrix_p = repmat(TRange_vec_p', T_steps_total, 1);

Individual_Temp_History = Tmax_matrix_p - Individual_SOC_History .* TRange_matrix_p;

P_standby = 0.05; 
Baseline_Power_History = (T_ja_participating_T - Tset_matrix_p) ./ (R_matrix_p .* eta_matrix_p);
Total_Power_History = Baseline_Power_History + Individual_Power_History;
Total_Power_History(Total_Power_History < P_standby) = P_standby;

Agg_Baseline_Power = sum(Baseline_Power_History, 2);
Agg_Total_Power = sum(Total_Power_History, 2);

SOC_Final = mean(CURRENT_SOC_AC, 'omitnan');
SOC_Sequence = [Agg_SOC_History; SOC_Final]; 

Agg_Model_Dev_Power = zeros(T_steps_total, 1);
if abs(AggParams.B) > 1e-9
    for t = 1:T_steps_total
        s_t = SOC_Sequence(t);
        s_next = SOC_Sequence(t+1);
        Agg_Model_Dev_Power(t) = (s_next - AggParams.A * s_t - AggParams.C) / AggParams.B;
    end
else
    warning('AggParams.B is close to 0, cannot inverse model.');
end

Agg_Model_Total_Power = Agg_Baseline_Power + Agg_Model_Dev_Power;

results = struct();
results.dt = dt;
results.time_points = time_points; 
results.Agg_P_Potential_Up_History = Agg_P_Potential_Up_History;     
results.Agg_P_Potential_Down_History = Agg_P_Potential_Down_History; 
results.Agg_Model_Potential_Up_History = Agg_Model_Potential_Up_History;     
results.Agg_Model_Potential_Down_History = Agg_Model_Potential_Down_History; 
results.Agg_P_Command_History = Agg_P_Command_History;   
results.Agg_P_Achieved_History = Agg_P_Achieved_History; 
results.Agg_Baseline_Power = Agg_Baseline_Power; 
results.Agg_Total_Power = Agg_Total_Power;       
results.Agg_Model_Total_Power = Agg_Model_Total_Power; 
results.Total_Power_History = Total_Power_History; 
results.Individual_Power_History = Individual_Power_History; 
results.Individual_SOC_History_Transposed = Individual_SOC_History'; 
results.Individual_Temp_History_Transposed = Individual_Temp_History';
results.AC_Up_Individual = AC_Up_Individual;
results.AC_Down_Individual = AC_Down_Individual;

% [修改] 输出文件添加 8am 后缀
output_mat_name = 'AC_Stateful_Simulation_Results_5min_pi_8am.mat';
save(output_mat_name, 'results', '-v7.3');
fprintf('  完整数据已保存至: %s\n', output_mat_name);