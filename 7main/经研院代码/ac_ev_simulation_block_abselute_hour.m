% ac_ev_simulation_block_abselute_hour.m
% 功能：全功能虚拟电厂调节潜力仿真系统（分块处理 + 绝对时间轴）
% 说明：
% 1. 支持大规模数据的分块加载与处理。
% 2. 使用绝对时间轴（例如 06:00 至次日 06:00）。
% 3. 保存详细的 EV 状态（E_current, P_current, SOC）及时间信息。

clear; close all; clc;

tic; % 启动总计时器

%% 1. 系统初始化 (全局参数)
rng(2023);                                      % 固定随机种子

% --- 设定绝对时间轴 ---
simulation_start_hour = 6;  % 仿真开始时间 (06:00)
simulation_end_hour   = 30; % 仿真结束时间 (次日 06:00)
dt = 1/60;                  % 时间分辨率 (小时)
time_points_absolute = simulation_start_hour:dt:simulation_end_hour;

base_price = 30;            % 基础电价 (元/kWh)
t_adj = 1;                  % 调节时长 (小时)

%% 1.5. 模块运行控制
runAC = 0;     % 是否运行 AC 仿真
runEV = true;  % 是否运行 EV 仿真

%% 1.6. 分块范围控制
startChunk = 12;    % 起始分块编号
endChunk = 12;      % 结束分块编号 (设为 inf 可运行至文件末尾)

if ~runAC && ~runEV
    error('AC 和 EV 仿真均被禁用。请至少启用一个模块。');
end

%% 2. 初始化文件路径
if runEV
    evFile = '2EV_residential.xlsx';
end
if runAC
    acFile = 'AC_template1.xlsx';
end

outputDir = 'chunk_results_abs_hour_test';      % 结果输出目录
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% 3. 主分块循环
chunkSize = 10;  % 每次处理的设备数量 (需匹配Excel数据结构)
chunkIndex = startChunk;

while true
    fprintf('\n====================================================\n');
    fprintf('== 准备处理分块 %d (设备 %d 到 %d) ==\n', ...
        chunkIndex, (chunkIndex-1)*chunkSize + 1, chunkIndex*chunkSize);
    fprintf('====================================================\n');

    % 检查是否超出设定的结束分块
    if chunkIndex > endChunk
        fprintf('已达到设定的结束块 %d，仿真停止。\n', endChunk);
        break;
    end

    %% 3.1. 分块读取设备参数

    % 计算读取范围
    startRow_Excel = 1 + (chunkIndex - 1) * chunkSize + 1;
    endRow_Excel = startRow_Excel + chunkSize - 1;
    dataRange = sprintf('%d:%d', startRow_Excel, endRow_Excel);

    ACs = [];
    EVs = [];

    % --- 读取 AC 数据 ---
    if runAC
        try
            acOpts = detectImportOptions(acFile, 'ReadVariableNames', true);
            acOpts.DataRange = dataRange;
            acTable = readtable(acFile, acOpts);
            if ~isempty(acTable)
                ACs = table2struct(acTable);
            end
        catch ME
            fprintf('读取 AC 文件出错或到达末尾: %s\n', ME.message);
            ACs = [];
        end
    end

    % --- 读取 EV 数据 ---
    if runEV
        try
            evOpts = detectImportOptions(evFile, 'ReadVariableNames', true);
            evOpts.DataRange = dataRange;
            evTable = readtable(evFile, evOpts);
            if ~isempty(evTable)
                EVs = table2struct(evTable);
            end
        catch ME
            fprintf('读取 EV 文件出错或到达末尾: %s\n', ME.message);
            EVs = [];
        end
    end

    % --- 检查循环终止条件 ---
    if runAC && runEV && (isempty(ACs) || isempty(EVs))
        fprintf('数据读取完毕，处理结束。\n'); break;
    elseif runAC && ~runEV && isempty(ACs)
        fprintf('AC 数据读取完毕，处理结束。\n'); break;
    elseif ~runAC && runEV && isempty(EVs)
        fprintf('EV 数据读取完毕，处理结束。\n'); break;
    end

    if runAC, num_AC = length(ACs); else, num_AC = 0; end
    if runEV, num_EV = length(EVs); else, num_EV = 0; end

    fprintf('分块 %d: 加载 %d 台 AC 和 %d 台 EV。\n', chunkIndex, num_AC, num_EV);

    %% 4. 激励响应模块
    p_min = 15; p_max = 50;
    p_min_prime = 10; p_max_prime = 40; T_set_max = 3;

    if runAC
        fprintf('分块 %d: 计算 AC 激励响应...\n', chunkIndex);
        temp_ACs = ACs;
        parfor i = 1:num_AC
            participation = calculateParticipation(temp_ACs(i).p_incentive, base_price);
            [~, ~, deltaT] = incentiveTempAC(...
                temp_ACs(i).p_incentive, p_min, p_max, p_min_prime, p_max_prime, T_set_max);

            temp_ACs(i).ptcp = (rand() < participation);

            if temp_ACs(i).ptcp
                temp_ACs(i).Tmax = temp_ACs(i).Tset + deltaT;
                temp_ACs(i).Tmin = temp_ACs(i).Tset - deltaT;
            else
                temp_ACs(i).Tmax = temp_ACs(i).Tset;
                temp_ACs(i).Tmin = temp_ACs(i).Tset;
            end

            % 生成环境温度 (基于绝对时间偏移)
            base_temp = temp_ACs(i).Tset + 4*sin(2*pi*(time_points_absolute-simulation_start_hour)/24);
            temp_range = temp_ACs(i).Tmax - temp_ACs(i).Tmin;

            if temp_range == 0, temp_range = 0.1; end

            noise = 0.2 * temp_range * randn(size(time_points_absolute));

            if temp_ACs(i).ptcp
                temp_ACs(i).T_ja = min(max(base_temp + noise, temp_ACs(i).Tmin), temp_ACs(i).Tmax);
            else
                temp_ACs(i).T_ja = base_temp + noise;
            end
        end
        ACs = temp_ACs;
    end

    if runEV
        fprintf('分块 %d: 计算 EV 激励响应...\n', chunkIndex);
        temp_EVs = EVs;
        parfor i = 1:num_EV
            temp_EVs(i).p_incentive = 35;
            participation = calculateParticipation(temp_EVs(i).p_incentive, base_price);
            temp_EVs(i).ptcp = (rand() < participation);
            temp_EVs(i).E_tar_original = temp_EVs(i).E_tar;

            E_flex_max = 0.2 * temp_EVs(i).C_EV;
            [deltaE_up, deltaE_down] = incentiveTempEV_updown(temp_EVs(i).p_incentive, p_min, p_max, p_min_prime, p_max_prime, E_flex_max);

            if temp_EVs(i).ptcp
                temp_EVs(i).E_reg_min = max(temp_EVs(i).E_in, temp_EVs(i).E_tar_original - deltaE_down);
                temp_EVs(i).E_reg_max = min(temp_EVs(i).C_EV, temp_EVs(i).E_tar_original + deltaE_up);
            else
                temp_EVs(i).E_reg_min = temp_EVs(i).E_tar_original;
                temp_EVs(i).E_reg_max = temp_EVs(i).E_tar_original;
            end
        end
        EVs = temp_EVs;
    end

    %% 5. 预计算模块
    num_time_points_for_baseline = length(time_points_absolute);

    if runEV
        fprintf('分块 %d: 预计算 EV 基线...\n', chunkIndex);
        temp_EVs = EVs;
        parfor i = 1:num_EV
            temp_EVs(i).P_base_sequence = EVbaseP_ChargeUntilFull(...
                temp_EVs(i).C_EV, temp_EVs(i).eta,...
                temp_EVs(i).E_tar_original, temp_EVs(i).E_in,...
                temp_EVs(i).t_dep, temp_EVs(i).t_in, dt, ...
                temp_EVs(i).r, temp_EVs(i).p_on, temp_EVs(i).SOC, num_time_points_for_baseline, time_points_absolute);
        end
        EVs = temp_EVs;
    end

    if runAC
        fprintf('分块 %d: 预计算 AC 热力参数...\n', chunkIndex);
        temp_ACs = ACs;
        parfor i = 1:num_AC
            [alpha, beta, gamma] = calculateACABC_single(...
                temp_ACs(i).R, temp_ACs(i).C, temp_ACs(i).eta,...
                temp_ACs(i).Tmax, temp_ACs(i).Tmin, temp_ACs(i).Tset, dt);
            temp_ACs(i).alpha = alpha;
            temp_ACs(i).beta = beta;
            temp_ACs(i).gamma = gamma;
        end
        ACs = temp_ACs;
    end

    %% 6. 主时间循环
    % 6.1 结果预分配
    num_steps = length(time_points_absolute);

    if runAC
        AC_Up = zeros(num_steps,1);
        AC_Down = zeros(num_steps,1);
        SOC_AC = zeros(num_AC, num_steps);
        AC_Up_Individual = zeros(num_AC, num_steps);
        AC_Down_Individual = zeros(num_AC, num_steps);
    end

    if runEV
        EV_Up = zeros(num_steps,1);
        EV_Down = zeros(num_steps,1);
        SOC_EV = zeros(num_EV, num_steps);
        m3 = zeros(num_EV,1);
        EV_Up_Individual = zeros(num_EV, num_steps);
        EV_Down_Individual = zeros(num_EV, num_steps);

        % 额外状态存储
        E_current_EV = zeros(num_EV, num_steps);
        P_current_EV = zeros(num_EV, num_steps);

        if ~isfield(EVs, 'E_exp'), [EVs.E_exp] = deal([EVs.E_in]); end
        if ~isfield(EVs, 'E_current'), [EVs.E_current] = deal([EVs.E_in]); end
        if ~isfield(EVs, 'P_current'), [EVs.P_current] = deal(0); end
    end

    % 6.2 时间步进循环
    fprintf('开始分块 %d 的时序仿真...\n', chunkIndex);

    for t_idx = 1:length(time_points_absolute)
        current_absolute_hour = time_points_absolute(t_idx);

        if mod(t_idx, round(length(time_points_absolute)/20)) == 1 || t_idx == length(time_points_absolute)
            fprintf('  分块 %d: 仿真时间 %.2f h (%.1f%%)\n', ...
                chunkIndex, current_absolute_hour, 100*t_idx/length(time_points_absolute));
        end

        %% EV 集群处理
        if runEV
            temp_EVs = EVs;
            temp_EV_Up = 0;
            temp_EV_Down = 0;
            temp_SOC_EV = zeros(num_EV,1);

            temp_E_current_EV = zeros(num_EV, 1);
            temp_P_current_EV = zeros(num_EV, 1);
            temp_m3 = m3;

            parfor i = 1:num_EV
                % 判断是否在线 (使用绝对时间)
                online = (current_absolute_hour >= temp_EVs(i).t_in) && (current_absolute_hour < temp_EVs(i).t_dep);

                if temp_EVs(i).ptcp
                    if online
                        % 在线：计算状态与潜力
                        temp_EVs(i).P_base = temp_EVs(i).P_base_sequence(t_idx);

                        [~, ~, m3_val] = calculateEVABC_single(temp_EVs(i).C_EV, temp_EVs(i).eta, temp_EVs(i).E_tar_original, temp_EVs(i).E_in, temp_EVs(i).t_dep, temp_EVs(i).t_in, dt, temp_EVs(i).r);
                        temp_m3(i) = m3_val;

                        [temp_EVs(i).E_exp, temp_EVs(i).E_current, temp_EVs(i).P_current, temp_EVs(i).SOC] = ...
                            calculateEVS_single(m3_val, temp_EVs(i).E_exp, temp_EVs(i).E_tar_original, temp_EVs(i).eta, temp_EVs(i).E_current, temp_EVs(i).P_current, temp_EVs(i).C_EV, temp_EVs(i).r, temp_EVs(i).p_on, dt, temp_EVs(i).t_dep, current_absolute_hour);

                        [DeltaP_plus, DeltaP_minus] = calculateEVAdjustmentPotentia_new(temp_EVs(i).E_reg_min, temp_EVs(i).E_reg_max, temp_EVs(i).E_current, temp_EVs(i).t_dep, current_absolute_hour, temp_EVs(i).p_on, temp_EVs(i).P_base, temp_EVs(i).eta, t_adj);

                        temp_EV_Up = temp_EV_Up + DeltaP_plus;
                        temp_EV_Down = temp_EV_Down + DeltaP_minus;

                        EV_Up_Individual(i, t_idx) = DeltaP_plus;
                        EV_Down_Individual(i, t_idx) = DeltaP_minus;
                        temp_SOC_EV(i) = temp_EVs(i).SOC;

                        temp_E_current_EV(i) = temp_EVs(i).E_current;
                        temp_P_current_EV(i) = temp_EVs(i).P_current;

                    else
                        % 离线：保持状态
                        if t_idx > 1
                            temp_SOC_EV(i) = SOC_EV(i, t_idx-1);
                            temp_E_current_EV(i) = E_current_EV(i, t_idx-1);
                            temp_P_current_EV(i) = P_current_EV(i, t_idx-1);
                        else
                            temp_SOC_EV(i) = temp_EVs(i).SOC;
                            temp_E_current_EV(i) = temp_EVs(i).E_in;
                            temp_P_current_EV(i) = 0;
                        end
                        EV_Up_Individual(i, t_idx)   = 0;
                        EV_Down_Individual(i, t_idx) = 0;

                        if t_idx > 1
                            temp_EVs(i).E_current = EVs(i).E_current;
                            temp_EVs(i).P_current = EVs(i).P_current;
                            temp_EVs(i).E_exp     = EVs(i).E_exp;
                        end
                    end
                else
                    % 不参与
                    if t_idx > 1
                        temp_SOC_EV(i) = SOC_EV(i, t_idx-1);
                        temp_E_current_EV(i) = E_current_EV(i, t_idx-1);
                        temp_P_current_EV(i) = P_current_EV(i, t_idx-1);
                    else
                        temp_SOC_EV(i) = temp_EVs(i).SOC;
                        temp_E_current_EV(i) = temp_EVs(i).E_in;
                        temp_P_current_EV(i) = 0;
                    end
                    EV_Up_Individual(i, t_idx)   = 0;
                    EV_Down_Individual(i, t_idx) = 0;

                    if t_idx > 1
                        temp_EVs(i).E_current = EVs(i).E_current;
                        temp_EVs(i).P_current = EVs(i).P_current;
                        temp_EVs(i).E_exp     = EVs(i).E_exp;
                    end
                end
            end

            EVs = temp_EVs;
            m3 = temp_m3;
            EV_Up(t_idx) = temp_EV_Up;
            EV_Down(t_idx) = temp_EV_Down;
            SOC_EV(:, t_idx) = temp_SOC_EV;
            E_current_EV(:, t_idx) = temp_E_current_EV;
            P_current_EV(:, t_idx) = temp_P_current_EV;
        end

        %% AC 集群处理
        if runAC
            temp_ACs = ACs;
            temp_AC_Up = 0;
            temp_AC_Down = 0;
            temp_SOC_AC = zeros(num_AC,1);

            parfor i = 1:num_AC
                if temp_ACs(i).ptcp
                    % 在线：计算状态与潜力
                    temp_ACs(i).P_base = ACbaseP_single(temp_ACs(i).T_ja(t_idx), temp_ACs(i).Tset, temp_ACs(i).R, temp_ACs(i).eta);
                    temp_ACs(i).SOC = calculateACS_single(temp_ACs(i).T_ja(t_idx), temp_ACs(i).Tmax, temp_ACs(i).Tmin);

                    [DeltaP_plus, DeltaP_minus] = calculateACAdjustmentPotentia(temp_ACs(i).P_base, 2*abs(temp_ACs(i).P_base), 0, temp_ACs(i).alpha, temp_ACs(i).beta, temp_ACs(i).gamma, temp_ACs(i).SOC, dt);

                    temp_AC_Up = temp_AC_Up + DeltaP_plus;
                    temp_AC_Down = temp_AC_Down + DeltaP_minus;

                    AC_Up_Individual(i, t_idx) = DeltaP_plus;
                    AC_Down_Individual(i, t_idx) = DeltaP_minus;
                    temp_SOC_AC(i) = temp_ACs(i).SOC;
                else
                    % 不参与
                    if t_idx > 1
                        temp_SOC_AC(i) = SOC_AC(i, t_idx-1);
                    else
                        temp_SOC_AC(i) = calculateACS_single(temp_ACs(i).T_ja(1), temp_ACs(i).Tmax, temp_ACs(i).Tmin);
                    end
                    AC_Up_Individual(i, t_idx) = 0;
                    AC_Down_Individual(i, t_idx) = 0;
                end
            end

            ACs = temp_ACs;
            AC_Up(t_idx) = temp_AC_Up;
            AC_Down(t_idx) = temp_AC_Down;
            SOC_AC(:, t_idx) = temp_SOC_AC;
        end
    end
    fprintf('分块 %d 仿真完成。\n', chunkIndex);

    %% 7. 结果组装与保存
    results = struct();

    % 保存元数据
    results.dt = dt;
    results.time_points_absolute = time_points_absolute;

    if runAC
        results.AC_Up = AC_Up;
        results.AC_Down = AC_Down;
        results.SOC_AC = SOC_AC;
        results.AC_Up_Individual = AC_Up_Individual;
        results.AC_Down_Individual = AC_Down_Individual;
    end

    if runEV
        results.EV_Up = EV_Up;
        results.EV_Down = EV_Down;
        results.SOC_EV = SOC_EV;
        results.m3 = m3;
        results.EV_Up_Individual = EV_Up_Individual;
        results.EV_Down_Individual = EV_Down_Individual;
        results.E_current_EV = E_current_EV;
        results.P_current_EV = P_current_EV;
    end

    outputFileName = fullfile(outputDir, sprintf('results_chunk_%d.mat', chunkIndex));
    try
        fprintf('正在保存分块 %d 到 %s ...\n', chunkIndex, outputFileName);
        save(outputFileName, 'results', '-v7.3');
        fprintf('保存成功。\n');
    catch ME_save
        fprintf('*** 保存出错: %s ***\n尝试保存到本地 LAST_CHUNK_FAILED_SAVE.mat\n', ME_save.message);
        try
            save('LAST_CHUNK_FAILED_SAVE.mat', 'results', '-v7.3');
        catch ME_save2
            fprintf('*** 保存彻底失败: %s ***\n', ME_save2.message);
        end
    end

    %% 9. 清理内存并准备下一循环
    chunkIndex = chunkIndex + 1;

    clear results num_AC num_EV temp_ACs temp_EVs acTable evTable acOpts evOpts ...
        temp_AC_Up temp_AC_Down temp_EV_Up temp_EV_Down ...
        temp_SOC_AC temp_SOC_EV temp_m3;

    if runAC
        clear ACs AC_Up AC_Down SOC_AC AC_Up_Individual AC_Down_Individual;
    end
    if runEV
        clear EVs EV_Up EV_Down SOC_EV m3 EV_Up_Individual EV_Down_Individual E_current_EV P_current_EV;
    end

end

fprintf('\n====================================================\n');
fprintf('== 所有分块处理完毕。结果保存在 %s ==\n', outputDir);
fprintf('====================================================\n');

toc; % 结束总计时