%% 全功能虚拟电厂调节潜力分析系统 (内存优化版 + EV逻辑同步main_diff_delt_48_updown)
clear; close all; clc;

tic; % 开始计时

%% 1. 系统初始化
rng(2023);                                      % 固定随机种子，保证结果可重复
% T_total = 24; % 不再直接使用 T_total 定义时间跨度

% --- 时间参数定义 (与 main_diff_delt_48_updown.m 保持一致) ---
simulation_start_hour = 6;  % 第一天早上 6 点
simulation_end_hour   = 30; % 第二天早上 6 点 (24 + 6)

dt = 5/60;                                      % 时间分辨率（小时）
% T_total = simulation_end_hour - simulation_start_hour; % 计算得出总时长 (现在是 24 小时)

% 基于开始和结束时间生成仿真绝对时间点
time_points_absolute = simulation_start_hour:dt:simulation_end_hour;
num_time_points = length(time_points_absolute);

base_price = 30;                                % 基础电价（元/kWh）
t_adj = 60/60;                                  % 调节时长（小时）

%% 2. 初始化参数
evFile = '2EV_residential.xlsx';                % (使用您最新的文件名)
acFile = 'AC_template1.xlsx';                   % (使用您最新的文件名)

%% 3. 读取设备参数
fprintf('正在读取设备参数...\n');
ACs = initializeACsFromExcel(acFile);
EVs = initializeEVsFromExcel(evFile);
num_AC = length(ACs);
num_EV = length(EVs);
fprintf('读取完成: %d 台 AC, %d 台 EV。\n', num_AC, num_EV);

% --- 为 EV 结构体添加必要的原始备份字段 (如果 initializeEVsFromExcel 没有做) ---
% 确保 SOC_original 和 E_tar_original 存在
for i = 1:num_EV
    if ~isfield(EVs(i), 'SOC_original') && isfield(EVs(i), 'SOC')
        EVs(i).SOC_original = EVs(i).SOC;
    elseif ~isfield(EVs(i), 'SOC_original') && ~isfield(EVs(i), 'SOC')
        EVs(i).SOC_original = EVs(i).E_in / EVs(i).C_EV; % 估算一个初始SOC
        EVs(i).SOC = EVs(i).SOC_original;
    end
    if ~isfield(EVs(i), 'E_tar_original')
        EVs(i).E_tar_original = EVs(i).E_tar;
    end
    % 初始化运行状态
    EVs(i).E_current = EVs(i).E_in;
    EVs(i).E_exp     = EVs(i).E_in;
    EVs(i).P_current = 0;
end

%% 4. 激励响应模块
fprintf('正在计算激励响应 (并行)...\n');
%% 4.1 参数设定
p_min = 15; p_max = 50;                         % 原始电价范围
p_min_prime = 10; p_max_prime = 40; T_set_max = 3; % 调整后电价范围

%% AC温度设定调整 (与原 slow 版本逻辑一致)
temp_ACs = ACs;
parfor i = 1:num_AC
    participation = calculateParticipation(temp_ACs(i).p_incentive, base_price);
    [~, ~, deltaT] = incentiveTempAC(...
        temp_ACs(i).p_incentive, p_min, p_max, p_min_prime, p_max_prime, T_set_max);
    temp_ACs(i).ptcp = (rand() < participation);
    if temp_ACs(i).ptcp
        temp_ACs(i).Tmax = temp_ACs(i).Tset + deltaT;
        temp_ACs(i).Tmin = temp_ACs(i).Tset - deltaT;
    else % 确保未参与的也有 Tmax/Tmin
        temp_ACs(i).Tmax = temp_ACs(i).Tset;
        temp_ACs(i).Tmin = temp_ACs(i).Tset;
    end
    % --- 使用绝对时间生成环境温度 ---
    base_temp = temp_ACs(i).Tset + 4*sin(2*pi*(time_points_absolute-simulation_start_hour)/24); % T_ja 基于绝对时间
    temp_range = temp_ACs(i).Tmax - temp_ACs(i).Tmin;
    if abs(temp_range) < 1e-6; temp_range = 0.1; end % 防止除零
    noise = 0.2 * temp_range * randn(size(time_points_absolute));
    temp_ACs(i).T_ja = min(max(base_temp + noise, temp_ACs(i).Tmin), temp_ACs(i).Tmax);
end
ACs = temp_ACs;
clear temp_ACs;

%% 4.2 EV目标电量调整 (!!! 同步 main_diff_delt_48_updown.m !!!)
temp_EVs = EVs;
parfor i = 1:num_EV
    temp_EVs(i).p_incentive = 20; % 或根据需要设置激励电价
    participation = calculateParticipation(temp_EVs(i).p_incentive, base_price);
    temp_EVs(i).ptcp = (rand() < participation);

    % 使用 E_tar_original (已在步骤3中确保存在)
    % temp_EVs(i).E_tar_original = temp_EVs(i).E_tar; % 不再需要，假设已初始化

    E_flex_max = 0.2 * temp_EVs(i).C_EV;
    [deltaE_up, deltaE_down] = incentiveTempEV_updown(temp_EVs(i).p_incentive, p_min, p_max, p_min_prime, p_max_prime, E_flex_max);

    if temp_EVs(i).ptcp
        temp_EVs(i).E_reg_min = temp_EVs(i).E_tar_original - deltaE_down;
         if temp_EVs(i).E_reg_min <= temp_EVs(i).E_in
             temp_EVs(i).E_reg_min = temp_EVs(i).E_in;
         end
        temp_EVs(i).E_reg_max = temp_EVs(i).E_tar_original + deltaE_up;
         if temp_EVs(i).E_reg_max >= temp_EVs(i).C_EV
             temp_EVs(i).E_reg_max = temp_EVs(i).C_EV;
         end
    else
        temp_EVs(i).E_reg_min = temp_EVs(i).E_tar_original;
        temp_EVs(i).E_reg_max = temp_EVs(i).E_tar_original;
    end
    % E_tar 本身不再修改，使用 E_reg_min/max 和 E_tar_original
end
EVs = temp_EVs;
clear temp_EVs;

%% 5. 预计算模块
fprintf('正在执行预计算 (并行)...\n');
%% 5.1 EV基线功率计算 (!!! 同步 main_diff_delt_48_updown.m !!!)
temp_EVs = EVs;
parfor i = 1:num_EV
    % !!! 调用 EVbaseP_ChargeUntilFull 时传入绝对时间序列 !!!
    temp_EVs(i).P_base_sequence = EVbaseP_ChargeUntilFull(...
        temp_EVs(i).C_EV, temp_EVs(i).eta,...
        temp_EVs(i).E_tar_original, temp_EVs(i).E_in,...
        temp_EVs(i).t_dep, temp_EVs(i).t_in, dt, ...
        temp_EVs(i).r, temp_EVs(i).p_on, temp_EVs(i).SOC_original, num_time_points, time_points_absolute); % 传入绝对时间
end
EVs = temp_EVs;
clear temp_EVs;

%% 5.2 AC 预计算 (与原 slow 版本逻辑一致)
temp_ACs = ACs;
parfor i = 1:num_AC
    [alpha, beta, gamma] = calculateACABC_single(temp_ACs(i).R, temp_ACs(i).C, temp_ACs(i).eta, temp_ACs(i).Tmax, temp_ACs(i).Tmin, temp_ACs(i).Tset, dt);
    temp_ACs(i).alpha = alpha;
    temp_ACs(i).beta = beta;
    temp_ACs(i).gamma = gamma;
end
ACs = temp_ACs;
clear temp_ACs;

%% 5.5 结构体解构 (!!! 关键内存优化：解决序列化错误 !!!)
fprintf('正在解构结构体以优化内存...\n');

% --- AC ---
AC_ptcp     = logical([ACs.ptcp]);
AC_T_ja     = cat(1, ACs.T_ja); % M x T 矩阵
AC_Tset     = [ACs.Tset]';
AC_R        = [ACs.R]';
AC_eta      = [ACs.eta]';
AC_Tmax     = [ACs.Tmax]';
AC_Tmin     = [ACs.Tmin]';
AC_alpha    = [ACs.alpha]';
AC_beta     = [ACs.beta]';
AC_gamma    = [ACs.gamma]';
clear ACs; % !!! 释放原始结构体内存 !!!

% --- EV ---
EV_ptcp       = logical([EVs.ptcp]);
EV_P_base_seq = cat(2, EVs.P_base_sequence)'; % 转置为 M x T 矩阵
EV_C_EV       = [EVs.C_EV]';
EV_eta        = [EVs.eta]';
EV_E_tar_orig = [EVs.E_tar_original]';
EV_E_in       = [EVs.E_in]';
EV_t_dep      = [EVs.t_dep]';
EV_t_in       = [EVs.t_in]';
EV_r          = [EVs.r]';
EV_p_on       = [EVs.p_on]';
EV_E_reg_min  = [EVs.E_reg_min]';
EV_E_reg_max  = [EVs.E_reg_max]';
EV_E_exp      = [EVs.E_exp]';     % 初始 E_exp
EV_E_current  = [EVs.E_current]'; % 初始 E_current
EV_P_current  = [EVs.P_current]'; % 初始 P_current
EV_SOC_orig   = [EVs.SOC_original]'; % 初始 SOC

% --- 计算 m3 (只需计算一次) ---
m3 = zeros(num_EV,1);
temp_m3 = m3;
parfor i = 1:num_EV
    % m3 的计算仅依赖于原始参数，与 ptcp 无关
    [~, ~, m3_val] = calculateEVABC_single(EV_C_EV(i), EV_eta(i), EV_E_tar_orig(i), EV_E_in(i), EV_t_dep(i), EV_t_in(i), dt, EV_r(i));
    temp_m3(i) = m3_val;
end
m3 = temp_m3;
clear temp_m3 EVs; % !!! 释放原始结构体内存 !!!
fprintf('解构完成。\n');

%% 6. 主时间循环
%% 6.1 结果预分配 (!!! 硬盘存储修改 !!!)

% 聚合结果 (保留在内存中)
AC_Up = zeros(num_time_points,1);
AC_Down = zeros(num_time_points,1);
EV_Up = zeros(num_time_points,1);
EV_Down = zeros(num_time_points,1);
EV_Total_Charge_Power = zeros(num_time_points, 1); % 新增：总充电功率
Online_EV_Count = zeros(num_time_points, 1);       % 新增：在线车辆数

% --- 个体结果：创建 matfile 对象 (硬盘映射) ---
individual_results_file = 'individual_results_slow_synced.mat';
fprintf('创建硬盘映射文件: %s\n', individual_results_file);
try
    if exist(individual_results_file, 'file'), delete(individual_results_file); end
    m = matfile(individual_results_file, 'Writable', true);

    fprintf('正在硬盘上预分配空间 (这可能需要几分钟)...\n');
    m.AC_Up_Individual = zeros(num_AC, num_time_points, 'single');
    m.AC_Down_Individual = zeros(num_AC, num_time_points, 'single');
    m.SOC_AC = zeros(num_AC, num_time_points, 'single');

    m.EV_Up_Individual = zeros(num_EV, num_time_points, 'single');
    m.EV_Down_Individual = zeros(num_EV, num_time_points, 'single');
    m.SOC_EV = zeros(num_EV, num_time_points, 'single');
    m.E_current_EV = zeros(num_EV, num_time_points, 'single'); % 存储 E_current 变化

    % 存储仿真元数据
    m.time_points_absolute = time_points_absolute;
    m.dt = dt;
    m.simulation_start_hour = simulation_start_hour;
    m.t_adj = t_adj;

catch E
    fprintf('创建 matfile 失败! 请检查硬盘权限或空间。\n');
    rethrow(E);
end
fprintf('硬盘预分配完成。\n');

%% 6.2 时间步进循环 (!!! 硬盘存储修改 !!!)

% --- 在主循环外初始化 EV 状态数组 ---
% 这些数组将在每个时间步被 parfor 更新
temp_EV_E_exp     = EV_E_exp;     % M x 1
temp_EV_E_current = EV_E_current; % M x 1
temp_EV_P_current = EV_P_current; % M x 1
temp_SOC_EV       = EV_SOC_orig;  % M x 1, 初始化为原始SOC

fprintf('开始主时间循环...\n');
for t_idx = 1:num_time_points
    current_absolute_hour = time_points_absolute(t_idx); % 当前绝对时间

    if mod(t_idx-1, round(1/dt)) == 0 % 每小时打印一次进度
        fprintf('== 正在处理绝对时间 %.2f/%.2f 小时 (%.1f%%) ==\n', ...
            current_absolute_hour, simulation_end_hour, (t_idx/num_time_points)*100);
    end

    % --- 为当前时间步准备临时存储 (在 parfor 内部赋值) ---
    temp_EV_Up_Ind   = zeros(num_EV, 1, 'single');
    temp_EV_Down_Ind = zeros(num_EV, 1, 'single');
    temp_SOC_EV_Ind  = zeros(num_EV, 1, 'single');
    temp_AC_Up_Ind   = zeros(num_AC, 1, 'single');
    temp_AC_Down_Ind = zeros(num_AC, 1, 'single');
    temp_SOC_AC_Ind  = zeros(num_AC, 1, 'single');

    % --- 准备下一时间步的状态 (将在 parfor 中被计算和赋值) ---
    next_EV_E_exp     = temp_EV_E_exp;     % 先拷贝当前状态
    next_EV_E_current = temp_EV_E_current;
    next_EV_P_current = temp_EV_P_current;
    next_SOC_EV       = temp_SOC_EV;

    %% 电动汽车集群处理 (!!! 同步 main_diff_delt_48_updown.m !!!)
    temp_EV_Up = 0;
    temp_EV_Down = 0;
    temp_EV_Total_Charge_Power = 0; % 当前时间步总充电功率
    temp_Online_EV_Count = 0;       % 当前时间步在线车辆数

    % --- 获取当前时间步的所有 EV 基线功率 (M x 1 向量) ---
    current_P_base_col_vector = EV_P_base_seq(:, t_idx);

    parfor i = 1:num_EV
        online = (current_absolute_hour >= EV_t_in(i)) && (current_absolute_hour < EV_t_dep(i));

        if EV_ptcp(i)
            if online
                temp_Online_EV_Count = temp_Online_EV_Count + 1; % 原子操作或在 parfor 外累加更安全，但这里影响不大

                P_base_i = current_P_base_col_vector(i);

                % 调用 calculateEVS_single 更新状态
                % 输入的是当前步开始时的状态 temp_EV_...
                [E_exp_i, E_current_i, P_current_i, SOC_i] = ...
                    calculateEVS_single(m3(i), temp_EV_E_exp(i), EV_E_tar_orig(i),...
                    EV_eta(i), temp_EV_E_current(i), temp_EV_P_current(i), EV_C_EV(i),...
                    EV_r(i), EV_p_on(i), dt, EV_t_dep(i), current_absolute_hour); % 使用绝对时间

                % 累加充电功率
                if P_current_i > 1e-3
                     temp_EV_Total_Charge_Power = temp_EV_Total_Charge_Power + P_current_i;
                end

                % 调用 calculateEVAdjustmentPotentia_new 计算潜力
                [DeltaP_plus, DeltaP_minus] = calculateEVAdjustmentPotentia_new(...
                    EV_E_reg_min(i), EV_E_reg_max(i), ...
                    E_current_i, EV_t_dep(i), current_absolute_hour, ... % 使用绝对时间
                    EV_p_on(i), P_base_i, EV_eta(i), t_adj); % 使用 t_adj

                temp_EV_Up = temp_EV_Up + DeltaP_plus;
                temp_EV_Down = temp_EV_Down + DeltaP_minus;

                % 存储个体结果到临时数组
                temp_EV_Up_Ind(i)   = DeltaP_plus;
                temp_EV_Down_Ind(i) = DeltaP_minus;
                temp_SOC_EV_Ind(i)  = SOC_i;

                % --- 存储计算出的下一时间步状态 ---
                next_EV_E_exp(i)     = E_exp_i;
                next_EV_E_current(i) = E_current_i;
                next_EV_P_current(i) = P_current_i;
                next_SOC_EV(i)       = SOC_i;

            else % 离线
                % 保持上一时刻的 SOC
                 if t_idx > 1
                     % 从硬盘读取上一时刻的SOC (效率低，最好在内存中传递)
                     % temp_SOC_EV_Ind(i) = m.SOC_EV(i, t_idx-1); % 直接访问matfile可能慢
                     temp_SOC_EV_Ind(i) = temp_SOC_EV(i); % 直接用上一时刻内存中的值
                 else
                     temp_SOC_EV_Ind(i) = EV_SOC_orig(i); % 初始SOC
                 end
                 % next_SOC_EV(i) = temp_SOC_EV_Ind(i); % 确保 next_SOC_EV 也更新
                 next_SOC_EV(i) = temp_SOC_EV(i); % 使用上一时刻的值填充
            end
        else % 不参与
             % 保持上一时刻的 SOC
             if t_idx > 1
                 temp_SOC_EV_Ind(i) = temp_SOC_EV(i);
             else
                 temp_SOC_EV_Ind(i) = EV_SOC_orig(i);
             end
             next_SOC_EV(i) = temp_SOC_EV(i);
        end
    end % 结束 parfor EV

    % --- 在 parfor 外部更新聚合结果和状态变量 ---
    EV_Up(t_idx) = temp_EV_Up;
    EV_Down(t_idx) = temp_EV_Down;
    EV_Total_Charge_Power(t_idx) = temp_EV_Total_Charge_Power;
    Online_EV_Count(t_idx) = temp_Online_EV_Count;

    % --- 将个体结果写入硬盘 ---
    m.EV_Up_Individual(:, t_idx) = temp_EV_Up_Ind;
    m.EV_Down_Individual(:, t_idx) = temp_EV_Down_Ind;
    m.SOC_EV(:, t_idx) = next_SOC_EV; % 存储更新后的 SOC
    m.E_current_EV(:, t_idx) = next_EV_E_current; % 存储更新后的 E_current

    % --- 更新下一时间步的状态 ---
    temp_EV_E_exp     = next_EV_E_exp;
    temp_EV_E_current = next_EV_E_current;
    temp_EV_P_current = next_EV_P_current;
    temp_SOC_EV       = next_SOC_EV; % 更新内存中的 SOC 状态

    %% 空调分析 - (AC 逻辑保持不变, 使用绝对时间)
    temp_AC_Up = 0;
    temp_AC_Down = 0;

    % AC_T_ja 是 M x T 矩阵, 提取第 t_idx 列 (M x 1 向量)
    current_T_ja_AC = AC_T_ja(:, t_idx);

    parfor i = 1:num_AC
        if AC_ptcp(i)
            P_base_i = ACbaseP_single(current_T_ja_AC(i), AC_Tset(i), AC_R(i), AC_eta(i));
            SOC_i = calculateACS_single(current_T_ja_AC(i), AC_Tmax(i), AC_Tmin(i));

            [DeltaP_plus, DeltaP_minus] = calculateACAdjustmentPotentia(...
                P_base_i, 2*abs(P_base_i), 0,...
                AC_alpha(i), AC_beta(i), AC_gamma(i),...
                SOC_i, dt); % 注意这里的 dt 是时间步长，不是 t_adj

            temp_AC_Up = temp_AC_Up + DeltaP_plus;
            temp_AC_Down = temp_AC_Down + DeltaP_minus;

            temp_AC_Up_Ind(i)   = DeltaP_plus;
            temp_AC_Down_Ind(i) = DeltaP_minus;
            temp_SOC_AC_Ind(i)  = SOC_i;
        else
            % 记录不参与 AC 的 SOC (例如，基于上一时刻或初始值)
            if t_idx > 1
                % temp_SOC_AC_Ind(i) = m.SOC_AC(i, t_idx-1); % 效率低
                % 需要在内存中维护 AC 的 SOC 状态，类似于 EV
                % 简化处理：假设不参与的AC SOC保持不变（可能不准确）
                 temp_SOC_AC_Ind(i) = NaN; % 或其他标记
            else
                 temp_SOC_AC_Ind(i) = calculateACS_single(current_T_ja_AC(i), AC_Tmax(i), AC_Tmin(i)); % 计算初始SOC
            end
        end
    end % 结束 parfor AC

    AC_Up(t_idx) = temp_AC_Up;
    AC_Down(t_idx) = temp_AC_Down;

    m.AC_Up_Individual(:, t_idx) = temp_AC_Up_Ind;
    m.AC_Down_Individual(:, t_idx) = temp_AC_Down_Ind;
    m.SOC_AC(:, t_idx) = temp_SOC_AC_Ind;

end % 结束 for t_idx

fprintf('所有时间步处理完毕。\n');

%% 7. 将结果组装回结构体

% `results` 结构体只包含聚合数据
results = struct(...
    'AC_Up', AC_Up,...
    'AC_Down', AC_Down,...
    'EV_Up', EV_Up,...
    'EV_Down', EV_Down,...
    'EV_Total_Charge_Power', EV_Total_Charge_Power, ... % 添加总充电功率
    'Online_EV_Count', Online_EV_Count, ...           % 添加在线车辆数
    'm3', m3... % m3 是否仍需要？它在 EV 计算中用到，但可能不需要最终输出
    );

% 将聚合结果保存到单独的文件
aggregate_results_file = 'aggregate_results_slow_synced.mat';
save(aggregate_results_file, 'results');

fprintf('仿真完成。\n');
fprintf('======================================================\n');
fprintf('聚合结果已保存到: %s\n', aggregate_results_file);
fprintf('全量个体数据已保存到: %s\n', individual_results_file);
fprintf('======================================================\n');
toc; % 结束计时