%% AC_main_Yearly_Sim_Heating_Cooling.m
% 全年(8760小时)空调集群仿真：自动兼容制冷(正功率)与制热(负功率)
%
% 主要修改点：
% 1. T_total = 8760 小时。
% 2. 模拟了全年的环境温度变化 (T_ja)。
% 3. 在主循环中根据 P_base 的符号自动判断运行模式 (制冷/制热)。

clear; close all; clc;

tic; % 启动计时器

%% 1. 系统初始化
rng(2024, 'Threefry'); % 固定随机种子

% --- [修改 1] 时间设置为全年 ---
T_total = 8760;      % 365天 * 24小时
dt = 60/60;          % 时间步长：15分钟 (考虑到全年数据量大，建议用15min或60min，5min内存可能吃紧)
% dt = 5/60;         % 如果内存足够(>16GB)，可以使用 5min
time_points = 0:dt:T_total; 
T_steps_total = length(time_points);
steps_per_hour = round(1/dt);
num_hours = floor(T_steps_total / steps_per_hour);

base_price = 30; % 基础电价

%% 2. 初始化 AC 参数
acFile = 'AC_template_year.xlsx'; 
fprintf('正在初始化空调参数...\n');
try
    % 假设 initializeACsFromExcel 在路径中
    ACs = initializeACsFromExcel(acFile);
catch ME
    error('加载 %s 失败。请确保文件存在且 helper 函数在路径中。\n错误: %s', acFile, ME.message);
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
fprintf('加载了 %d 台空调进行全年仿真。\n', num_AC);

%% 3. 激励响应参数
p_min = 15; p_max = 50; p_min_prime = 10; p_max_prime = 40; T_set_max = 3;
current_p = 25.0; % 固定一个电价进行仿真

%% 4. 预计算：生成全年环境温度与单体参数

fprintf('  Step 4.1: 生成全年环境温度曲线...\n');
temp_ACs = ACs;

% --- [修改 2] 生成全年环境温度曲线 (模拟四季) ---
% 假设：
%   - 年均温 15℃
%   - 季节波动幅度 15℃ (夏天最高 ~30℃, 冬天最低 ~0℃)
%   - 日夜温差波动幅度 5℃
%   - 最冷在 1月中旬 (t=0附近)，最热在 7月中旬 (t=4380附近)
%   - 公式：T(t) = T_avg - T_season_amp * cos(2*pi*t/8760) + T_day_amp * sin(...)

T_avg_year = 15;
T_season_amp = 18; % 季节振幅
T_day_amp = 5;     % 日夜振幅

% 1. 季节趋势 (余弦波，t=0时为谷底-冬天，t=4380为峰值-夏天)
season_trend = -T_season_amp * cos(2*pi * time_points / 8760);

% 2. 日夜波动 (正弦波，峰值在下午)
% (time_points - 14) 使得最高温出现在下午14:00左右
day_fluctuation = T_day_amp * cos(2*pi * (time_points - 14) / 24);

% 3. 随机天气波动 (移动平均噪声)
white_noise = randn(1, T_steps_total + 100);
smooth_noise = movmean(white_noise, 4*steps_per_hour); % 4小时平滑
weather_noise = 2.0 * smooth_noise(1:T_steps_total); % 幅度2度

% 4. 合成全年温度
base_ambient_temp_unified = T_avg_year + season_trend + day_fluctuation + weather_noise;

% 绘图预览温度曲线 (可选)
% figure; plot(time_points, base_ambient_temp_unified); title('全年环境温度模拟'); xlabel('Hour'); ylabel('Temp (C)'); drawnow;

fprintf('  温度生成完毕: 最低 %.1f℃, 最高 %.1f℃\n', min(base_ambient_temp_unified), max(base_ambient_temp_unified));

% --- 4.2 计算单体参数 ---
parfor i = 1:num_AC
    % 简单的激励响应
    participation = calculateParticipation(current_p, base_price);
    [~, ~, deltaT_flex] = incentiveTempAC(current_p, p_min, p_max, p_min_prime, p_max_prime, T_set_max);
    temp_ACs(i).ptcp = (rand() < participation);

    if temp_ACs(i).ptcp
        temp_ACs(i).Tmax = temp_ACs(i).Tset_original + deltaT_flex;
        temp_ACs(i).Tmin = temp_ACs(i).Tset_original - deltaT_flex;
    end
    
    % 将全年的温度赋值给每台空调 (假设同一地区)
    temp_ACs(i).T_ja = base_ambient_temp_unified;

    % 计算热力学参数 alpha, beta, gamma
    % 注意：虽然 beta, gamma 在不同温度下可能微调，这里假设线性模型参数在全年常数化近似
    [alpha, beta, gamma] = calculateACABC_single(...
        temp_ACs(i).R, temp_ACs(i).C, temp_ACs(i).eta,...
        temp_ACs(i).Tmax, temp_ACs(i).Tmin, temp_ACs(i).Tset, dt);
    temp_ACs(i).alpha = alpha;
    temp_ACs(i).beta = beta;
    temp_ACs(i).gamma = gamma;
end
ACs = temp_ACs;

% --- 4.3 聚合参数 ---
ACs_participating = ACs([ACs.ptcp]);
num_AC_participating = length(ACs_participating);
if num_AC_participating == 0, error('无空调参与'); end

AggParams = calculateAggregatedACParams(ACs_participating);
fprintf('  聚合模型构建完成 (参与数: %d)\n', num_AC_participating);

% --- 4.4 生成全年电网指令 ---
% 这里简单生成一个随机指令，实际中可以是新能源波动
fprintf('  Step 4.3: 生成全年电网指令...\n');
P_grid_command_series = generate_hourly_regulation_signal(T_steps_total, steps_per_hour, 8760, num_AC_participating);


%% 5. 主时间循环 (状态化仿真)

% 5.1 初始化
CURRENT_SOC_AC = [ACs_participating.SOC]'; 

% 5.2 结果存储 (考虑到内存，仅存储聚合数据，不存储 Individual 大矩阵)
Agg_SOC_History = zeros(T_steps_total, 1);
Agg_P_Base_Total = zeros(T_steps_total, 1);     % 总基线功率 (有正有负)
Agg_P_Potential_Up = zeros(T_steps_total, 1);
Agg_P_Potential_Down = zeros(T_steps_total, 1);
Agg_P_Response = zeros(T_steps_total, 1);

fprintf('  Step 5: 开始 %d 步 (%.1f年) 的仿真...\n', T_steps_total, T_total/8760);

% 进度条显示辅助
print_interval = floor(T_steps_total / 20); 

for t_idx = 1:T_steps_total
    if mod(t_idx, print_interval) == 0
        fprintf('    进度: %.0f%% (时间: %.0f h)\n', t_idx/T_steps_total*100, time_points(t_idx));
    end

    % 1. 聚合状态
    SOC_agg_t = mean(CURRENT_SOC_AC, 'omitnan');
    Agg_SOC_History(t_idx) = SOC_agg_t;
    
    % 2. 获取指令
    Delta_P_S = P_grid_command_series(t_idx);
    
    % 3. 预测目标 SOC
    SOC_target_next = AggParams.A * SOC_agg_t + AggParams.B * Delta_P_S + AggParams.C;
    SOC_target_next = max(0, min(1, SOC_target_next));

    % 4. 并行计算单体响应
    temp_Up_sum = 0;
    temp_Down_sum = 0;
    temp_P_base_sum = 0;
    temp_P_resp_sum = 0;
    temp_SOC_next = zeros(num_AC_participating, 1);

    parfor i = 1:num_AC_participating
        ac_i = ACs_participating(i);
        soc_curr = CURRENT_SOC_AC(i);
        
        % A. 计算基线功率 P_base
        % 公式: P = (T_ja - T_set) / (R * eta)
        % 若 T_ja > T_set (夏), P_base > 0 (制冷)
        % 若 T_ja < T_set (冬), P_base < 0 (制热)
        T_amb = ac_i.T_ja(t_idx);
        P_base_i = ACbaseP_single(T_amb, ac_i.Tset, ac_i.R, ac_i.eta);
        
        temp_P_base_sum = temp_P_base_sum + P_base_i;

        % --- [修改 3] 动态确定功率边界 ---
        % 根据 P_base 的符号判断运行模式
        if P_base_i > 0 
            % [制冷模式] 功率范围 [0, 2*P_base]
            P_phys_max = 2 * P_base_i;
            P_phys_min = 0;
        else
            % [制热模式] 功率范围 [2*P_base, 0] (注意 2*P_base 是更小的负数)
            P_phys_max = 0;
            P_phys_min = 2 * P_base_i;
        end
        % (死区/过渡季处理: 如果 P_base 极小，则 max=min=0，自然处理)

        % B. 计算调节潜力
        [P_plus, P_minus] = calculateACAdjustmentPotentia(...
            P_base_i, P_phys_max, P_phys_min, ... % 传入动态边界
            ac_i.alpha, ac_i.beta, ac_i.gamma,...
            soc_curr, dt);

        temp_Up_sum = temp_Up_sum + P_plus;
        temp_Down_sum = temp_Down_sum + P_minus;

        % C. 计算响应与状态更新
        delta_Pj = 0;
        if abs(ac_i.beta) > 1e-9
            delta_Pj = (SOC_target_next - ac_i.alpha * soc_curr - ac_i.gamma) / ac_i.beta;
        end
        
        % 裁剪指令 (限制在物理潜力范围内)
        delta_Pj_clipped = max(P_minus, min(P_plus, delta_Pj));

        % 更新状态
        soc_next = updateACSOC_single(soc_curr, delta_Pj_clipped, ...
            ac_i.alpha, ac_i.beta, ac_i.gamma);

        temp_SOC_next(i) = soc_next;
        temp_P_resp_sum = temp_P_resp_sum + delta_Pj_clipped;
    end

    % 5. 存储聚合结果
    Agg_P_Base_Total(t_idx)     = temp_P_base_sum;
    Agg_P_Potential_Up(t_idx)   = temp_Up_sum;
    Agg_P_Potential_Down(t_idx) = temp_Down_sum;
    Agg_P_Response(t_idx)       = temp_P_resp_sum;
    
    % 更新状态
    CURRENT_SOC_AC = temp_SOC_next;
end

fprintf('  仿真完成。\n');

% %% 6. 结果可视化 (针对全年数据的简化绘图)
% 
% % 定义绘图用的时间轴 (小时)
% t_axis = time_points;
% 
% % 图1: 全年基线功率 (展示制冷/制热切换)
% figure('Name', '全年基线功率', 'Position', [100 100 1000 400]);
% plot(t_axis, Agg_P_Base_Total, 'k-', 'LineWidth', 1);
% yline(0, 'r--', 'Zero');
% xlim([0, T_total]);
% xlabel('时间 (小时)'); ylabel('集群总基线功率 (kW)');
% title('全年空调集群基线功率 (正=制冷, 负=制热)');
% grid on;
% % 标注季节
% text(100, max(Agg_P_Base_Total)*0.8, '冬季(制热)', 'Color','b');
% text(4380, max(Agg_P_Base_Total)*0.8, '夏季(制冷)', 'Color','r');
% text(8000, max(Agg_P_Base_Total)*0.8, '冬季(制热)', 'Color','b');
% 
% % 图2: 调节潜力包络
% figure('Name', '全年调节潜力', 'Position', [100 550 1000 400]);
% hold on;
% % 绘制上调潜力区域
% area(t_axis, Agg_P_Potential_Up, 'FaceColor', [0.6 0.8 1], 'EdgeColor', 'none', 'DisplayName', '上调潜力');
% % 绘制下调潜力区域
% area(t_axis, Agg_P_Potential_Down, 'FaceColor', [1 0.6 0.6], 'EdgeColor', 'none', 'DisplayName', '下调潜力');
% yline(0, 'k-');
% xlabel('时间 (小时)'); ylabel('调节潜力 (kW)');
% title('全年空调集群调节潜力 (Up/Down)');
% legend; xlim([0, T_total]); grid on;
% hold off;
% 
% % 图3: 局部缩放 (查看春秋季/过渡季情况)
% % 选取春季某段 (例如 2000-2200小时)
% zoom_range = (t_axis >= 2000) & (t_axis <= 2200);
% if any(zoom_range)
%     figure('Name', '过渡季细节', 'Position', [150 300 800 400]);
%     plot(t_axis(zoom_range), Agg_P_Base_Total(zoom_range), 'k-', 'LineWidth', 1.5);
%     xlabel('时间 (小时)'); ylabel('功率 (kW)');
%     title('过渡季(春季)基线功率细节：自动启停');
%     grid on;
% end

%% 7. 保存结果
results_struct = struct();
results_struct.P_base = Agg_P_Base_Total;
results_struct.P_up = Agg_P_Potential_Up;
results_struct.P_down = Agg_P_Potential_Down;
results_struct.SOC = Agg_SOC_History;
results_struct.Temperature_Profile = base_ambient_temp_unified;

save('AC_Yearly_Sim_Results.mat', 'results_struct', '-v7.3');
fprintf('  结果已保存至 AC_Yearly_Sim_Results.mat\n');
toc;