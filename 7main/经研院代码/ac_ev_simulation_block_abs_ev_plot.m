%% 全功能虚拟电厂调节潜力分析系统 (*** 分块处理版 + 模块可选 + 范围可选 ***)
% 这是一个主程序脚本 (Script)，不是函数。
% 您可以直接运行此文件。
% (*** 已按要求修改时间轴：仿真从第6小时运行到第30小时 ***)
%
% [!!! V2 修改：增加了 E_current_EV, P_current_EV, dt, 和 time_points_absolute 的保存 !!!]

clear; close all; clc;

tic; % 启动一个总计时器

%% 1. 系统初始化 (全局参数)
rng(2023);                                      % 固定随机种子，保证结果可重复

% --- 修改：定义绝对时间范围 (仿照 main_diff_delt_48_updown.m) ---
simulation_start_hour = 6;  % 仿真开始于 6:00
simulation_end_hour   = 30; % 仿真结束于次日 6:00 (24 + 6)
dt = 1/60;                                      % 时间分辨率（小时）
% 基于绝对时间生成时间序列
time_points_absolute = simulation_start_hour:dt:simulation_end_hour; 
% --- 修改结束 ---

base_price = 30;                                % 基础电价（元/kWh）
t_adj = 1;                                      % 调节时长（小时），从 new main file 引入

%% 1.5. *** 模块运行控制 ***
% --- 请在这里配置您想运行的模块 ---
runAC = 0;  % <--- 设置 true 以运行 AC 仿真
runEV = true;  % <--- 设置 true 以运行 EV 仿真
% ---------------------------------

%% 1.6. *** 新增：分块范围控制 (您的新需求) ***
% --- 请在这里配置您想运行的分块范围 ---
startChunk = 12;    % <--- 设置起始块编号 (例如: 1)
endChunk = 12;    % <--- [!!! 示例：仅运行第1块 !!!] (设置为 inf 可运行到文件末尾)
% ---------------------------------

if ~runAC && ~runEV
    error('AC 和 EV 仿真均被禁用。请至少启用一个模块 (runAC 或 runEV)。');
end

%% 2. 初始化参数 (全局参数)
if runEV
    evFile = '2EV_residential.xlsx';                 % 电动汽车数据文件
end
if runAC
    acFile = 'AC_template1.xlsx';      
end

outputDir = 'chunk_results_abs_hour_test';                   % *** 新增：用于存储结果的文件夹 ***
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
             
%% 3. *** 新增：主分块循环 ***
chunkSize = 10;  % *** 每次处理10000台设备 *** [!!! 修正：匹配您的Excel大小 !!!]

% *** 修改：从您指定的 startChunk 开始 ***
chunkIndex = startChunk;   

while true
    fprintf('\n====================================================\n');
    fprintf('== 准备处理分块 %d (设备 %d 到 %d) ==\n', ...
             chunkIndex, (chunkIndex-1)*chunkSize + 1, chunkIndex*chunkSize);
    fprintf('====================================================\n');
    
    % *** 新增：检查是否超出结束块范围 ***
    if chunkIndex > endChunk
        fprintf('已达到设定的结束块 %d，仿真停止。\n', endChunk);
        break; % 退出 while 循环
    end
    % *** 检查结束 ***
    
    
    %% 3.1. 分块读取设备参数 (*** 替换原始的 Section 3 ***)
    
    % --- 计算要读取的 Excel 行范围 ---
    startRow_Excel = 1 + (chunkIndex - 1) * chunkSize + 1; % +1 是为了跳过标题行
    endRow_Excel = startRow_Excel + chunkSize - 1;
    dataRange = sprintf('%d:%d', startRow_Excel, endRow_Excel);
    
    ACs = []; % 初始化为空
    EVs = []; % 初始化为空
    
    % --- 读取 AC 数据 ---
    if runAC
        try
            % fprintf('正在从 %s 读取行 %s...\n', acFile, dataRange); % (减少打印)
            acOpts = detectImportOptions(acFile, 'ReadVariableNames', true);
            acOpts.DataRange = dataRange; % 指定读取范围
            
            acTable = readtable(acFile, acOpts);
            
            if ~isempty(acTable)
                ACs = table2struct(acTable);
            end
            
        catch ME
            fprintf('读取空调文件时出错或到达文件末尾: %s\n', ME.message);
            ACs = []; % 确保出错时为空
        end
    end

    % --- 读取 EV 数据 ---
    if runEV
        try
            % fprintf('正在从 %s 读取行 %s...\n', evFile, dataRange); % (减少打印)
            evOpts = detectImportOptions(evFile, 'ReadVariableNames', true);
            evOpts.DataRange = dataRange; % 指定读取范围
            
            evTable = readtable(evFile, evOpts);
            
            if ~isempty(evTable)
                EVs = table2struct(evTable);
            end
            
        catch ME
            fprintf('读取EV文件时出错或到达文件末尾: %s\n', ME.message);
            EVs = []; % 确保出错时为空
        end
    end

    % --- 检查数据与循环退出条件 ---
    % (如果读不到数据，也要退出)
    if runAC && runEV && (isempty(ACs) || isempty(EVs))
         fprintf('空调或电动汽车数据为空 (已到达文件末尾或范围 %s 无数据)，处理结束。\n', dataRange);
         break;
    elseif runAC && ~runEV && isempty(ACs)
        fprintf('未在 %s 中找到更多AC数据 (已到达文件末尾或范围 %s 无数据)，处理结束。\n', acFile, dataRange);
        break; % 退出 while 循环
    elseif ~runAC && runEV && isempty(EVs)
        fprintf('未在 %s 中找到更多EV数据 (已到达文件末尾或范围 %s 无数据)，处理结束。\n', evFile, dataRange);
        break; % 退出 while 循环
    end
    
    % 获取当前分块的设备数量
    if runAC, num_AC = length(ACs); else, num_AC = 0; end
    if runEV, num_EV = length(EVs); else, num_EV = 0; end
                     
    fprintf('分块 %d: 成功加载 %d 台 AC 和 %d 台 EV。\n', chunkIndex, num_AC, num_EV);

    
    % =================================================================
    % == 以下是您原始脚本的核心逻辑 (Section 4, 5, 6) ==
    % =================================================================

    %% 4. 激励响应模块
    p_min = 15; p_max = 50;                         
    p_min_prime = 10; p_max_prime = 40; T_set_max = 3;                

    if runAC
        fprintf('分块 %d: 正在计算 AC 激励响应...\n', chunkIndex);
        %% AC温度设定调整
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
            
            % --- 修改：使用绝对时间轴和偏移量计算 T_ja ---
            % (仿照 ac_ev_simulation_improve.m)
            base_temp = temp_ACs(i).Tset + 4*sin(2*pi*(time_points_absolute-simulation_start_hour)/24); 
            temp_range = temp_ACs(i).Tmax - temp_ACs(i).Tmin;
            % --- 修改结束 ---
            
            if temp_range == 0 
                temp_range = 0.1; 
            end
            
            % --- 修改：使用绝对时间轴大小 ---
            noise = 0.2 * temp_range * randn(size(time_points_absolute));
            % --- 修改结束 ---
            
            if temp_ACs(i).ptcp
                temp_ACs(i).T_ja = min(max(base_temp + noise, temp_ACs(i).Tmin), temp_ACs(i).Tmax);
            else
                temp_ACs(i).T_ja = base_temp + noise; 
            end
        end
        ACs = temp_ACs;
    end

    if runEV
        fprintf('分块 %d: 正在计算 EV 激励响应...\n', chunkIndex);
        %% 4.2 EV目标电量调整
        temp_EVs = EVs; 
        parfor i = 1:num_EV
            temp_EVs(i).p_incentive = 35;
            participation = calculateParticipation(temp_EVs(i).p_incentive, base_price);
            ptcp_result = (rand() < participation);
            temp_EVs(i).ptcp = ptcp_result;
            temp_EVs(i).E_tar_original = temp_EVs(i).E_tar;
            
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
        end
        EVs = temp_EVs; 
    end

    %% 5. 预计算模块
    % --- 修改：使用绝对时间轴 ---
    % H = T_total; % 不再需要
    % H_steps = H / dt; % 不再需要
    num_time_points_for_baseline = length(time_points_absolute); 
    % --- 修改结束 ---

    if runEV
        fprintf('分块 %d: 正在预计算 EV 基线...\n', chunkIndex);
        %% 5.1 EV基线功率计算
        temp_EVs = EVs;
        parfor i = 1:num_EV
            % --- 修改：传入绝对时间轴 ---
            temp_EVs(i).P_base_sequence = EVbaseP_ChargeUntilFull(...
                temp_EVs(i).C_EV, temp_EVs(i).eta,...
                temp_EVs(i).E_tar_original, temp_EVs(i).E_in,... 
                temp_EVs(i).t_dep, temp_EVs(i).t_in, dt, ... 
                temp_EVs(i).r, temp_EVs(i).p_on, temp_EVs(i).SOC, num_time_points_for_baseline, time_points_absolute);
            % --- 修改结束 ---
        end
        EVs = temp_EVs;
    end

    if runAC
        fprintf('分块 %d: 正在预计算 AC 参数...\n', chunkIndex);
        %% 5.2 AC 预计算
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
    %% 6.1 结果预分配 (*** 根据模块动态分配 ***)
    
    % --- 修改：使用绝对时间轴的长度 ---
    num_steps = length(time_points_absolute);
    % --- 修改结束 ---

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
        
        % [!!! 新增 !!!] 增加 E_current 和 P_current 的存储
        E_current_EV = zeros(num_EV, num_steps);
        P_current_EV = zeros(num_EV, num_steps);
        
        if ~isfield(EVs, 'E_exp')
            for i=1:num_EV
                 EVs(i).E_exp = EVs(i).E_in;
             end
    
        end
        if ~isfield(EVs, 'E_current')
             for i=1:num_EV
                 EVs(i).E_current = EVs(i).E_in;
             end
        end
        if ~isfield(EVs, 'P_current')
            [EVs.P_current] = deal(0);
        end
    end

    %% 6.2 时间步进循环
    fprintf('开始分块 %d 的时间步进仿真...\n', chunkIndex);
    
    % --- 修改：使用绝对时间轴循环 ---
    for t_idx = 1:length(time_points_absolute)
        current_absolute_hour = time_points_absolute(t_idx); % 获取当前绝对时间
        
        if mod(t_idx, round(length(time_points_absolute)/20)) == 1 || t_idx == length(time_points_absolute)
             fprintf('  分块 %d: 仿真绝对时间 %.2f小时 (%.1f%%)\n', ...
                 chunkIndex, current_absolute_hour, 100*t_idx/length(time_points_absolute));
        end
    % --- 修改结束 ---
    
        %% 电动汽车集群处理
        if runEV
            temp_EVs = EVs; 
            temp_EV_Up = 0;
            temp_EV_Down = 0;
            temp_SOC_EV = zeros(num_EV,1);
            
            % [!!! 新增 !!!] 增加 E_current 和 P_current 的临时存储
            temp_E_current_EV = zeros(num_EV, 1);
            temp_P_current_EV = zeros(num_EV, 1);
            
            temp_m3 = m3;
            
            parfor i = 1:num_EV
                % --- 修改：使用绝对时间检查 online 状态 ---
                % (假设 t_in 和 t_dep 是绝对时间, e.g., 8.0, 30.0)
                online = (current_absolute_hour >= temp_EVs(i).t_in) && (current_absolute_hour < temp_EVs(i).t_dep); 
                % --- 修改结束 ---
                
                if temp_EVs(i).ptcp 
                    if online
                        temp_EVs(i).P_base = temp_EVs(i).P_base_sequence(t_idx);
                        
                        [~, ~, m3_val] = calculateEVABC_single(temp_EVs(i).C_EV, temp_EVs(i).eta, temp_EVs(i).E_tar_original, temp_EVs(i).E_in, temp_EVs(i).t_dep, temp_EVs(i).t_in, dt, temp_EVs(i).r);
                        temp_m3(i) = m3_val;
                        
                        % --- 修改：传入 current_absolute_hour ---
                        [temp_EVs(i).E_exp, temp_EVs(i).E_current, temp_EVs(i).P_current, temp_EVs(i).SOC] = ...
                            calculateEVS_single(m3_val, temp_EVs(i).E_exp, temp_EVs(i).E_tar_original, temp_EVs(i).eta, temp_EVs(i).E_current, temp_EVs(i).P_current, temp_EVs(i).C_EV, temp_EVs(i).r, temp_EVs(i).p_on, dt, temp_EVs(i).t_dep, current_absolute_hour);
                        
                        [DeltaP_plus, DeltaP_minus] = calculateEVAdjustmentPotentia_new(temp_EVs(i).E_reg_min, temp_EVs(i).E_reg_max, temp_EVs(i).E_current, temp_EVs(i).t_dep, current_absolute_hour, temp_EVs(i).p_on, temp_EVs(i).P_base, temp_EVs(i).eta, t_adj);
                        % --- 修改结束 ---
                        
                        temp_EV_Up = temp_EV_Up + DeltaP_plus;
                        temp_EV_Down = temp_EV_Down + DeltaP_minus;
                        
                        EV_Up_Individual(i, t_idx) = DeltaP_plus;
                        EV_Down_Individual(i, t_idx) = DeltaP_minus;
                        temp_SOC_EV(i) = temp_EVs(i).SOC;
                        
                        % [!!! 新增 !!!] 存储 E_current 和 P_current
                        temp_E_current_EV(i) = temp_EVs(i).E_current;
                        temp_P_current_EV(i) = temp_EVs(i).P_current;
                        
                    else
                        if t_idx > 1 
                            temp_SOC_EV(i) = SOC_EV(i, t_idx-1); 
                            % [!!! 新增 !!!] 继承上一时刻的 E 和 P
                            temp_E_current_EV(i) = E_current_EV(i, t_idx-1);
                            temp_P_current_EV(i) = P_current_EV(i, t_idx-1);
                        else 
                            temp_SOC_EV(i) = temp_EVs(i).SOC; 
                            % [!!! 新增 !!!] 使用初始 E 和 P
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
                else % 如果不参与 (ptcp == false)
                     if t_idx > 1 
                        temp_SOC_EV(i) = SOC_EV(i, t_idx-1);
                        % [!!! 新增 !!!] 继承上一时刻的 E 和 P
                        temp_E_current_EV(i) = E_current_EV(i, t_idx-1);
                        temp_P_current_EV(i) = P_current_EV(i, t_idx-1);
                     else 
                        temp_SOC_EV(i) = temp_EVs(i).SOC; 
                        % [!!! 新增 !!!] 使用初始 E 和 P
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
            % [!!! 新增 !!!] 将临时向量存入主矩阵
            E_current_EV(:, t_idx) = temp_E_current_EV;
            P_current_EV(:, t_idx) = temp_P_current_EV;
        end
        
        %% 空调分析
        if runAC
            temp_ACs = ACs; 
            temp_AC_Up = 0;
            temp_AC_Down = 0;
            temp_SOC_AC = zeros(num_AC,1);
            
            parfor i = 1:num_AC
                if temp_ACs(i).ptcp
                    % (AC逻辑不变，T_ja已在第4节基于绝对时间预算好)
                    temp_ACs(i).P_base = ACbaseP_single(temp_ACs(i).T_ja(t_idx), temp_ACs(i).Tset, temp_ACs(i).R, temp_ACs(i).eta);
                    temp_ACs(i).SOC = calculateACS_single(temp_ACs(i).T_ja(t_idx), temp_ACs(i).Tmax, temp_ACs(i).Tmin);
                    
                    [DeltaP_plus, DeltaP_minus] = calculateACAdjustmentPotentia(temp_ACs(i).P_base, 2*abs(temp_ACs(i).P_base), 0, temp_ACs(i).alpha, temp_ACs(i).beta, temp_ACs(i).gamma, temp_ACs(i).SOC, dt);
                    
                    temp_AC_Up = temp_AC_Up + DeltaP_plus;
                    temp_AC_Down = temp_AC_Down + DeltaP_minus;
                    
                    AC_Up_Individual(i, t_idx) = DeltaP_plus;
                    AC_Down_Individual(i, t_idx) = DeltaP_minus;
                    temp_SOC_AC(i) = temp_ACs(i).SOC;
                else
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

    %% 7. 将结果组装回结构体 (*** 动态组装 ***)
    results = struct(); % 创建一个空结构体
    
    % [!!! 新增 !!!] 保存 dt 和 time_points_absolute
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
        % [!!! 新增 !!!]
        results.E_current_EV = E_current_EV;
        results.P_current_EV = P_current_EV;
    end
    
    %% 8. *** 新增：保存当前分块的结果 ***
    outputFileName = fullfile(outputDir, sprintf('results_chunk_%d.mat', chunkIndex));
    try
        fprintf('正在保存分块 %d 的结果到 %s ...\n', chunkIndex, outputFileName);
        save(outputFileName, 'results', '-v7.3'); % 使用 -v7.3 以支持大文件
        fprintf('保存成功。\n');
    catch ME_save
        fprintf('*** 保存文件时出错: %s ***\n', ME_save.message);
        fprintf('将尝试在工作区保存... \n');
        try
            save('LAST_CHUNK_FAILED_SAVE.mat', 'results', '-v7.3');
            fprintf('已保存到 LAST_CHUNK_FAILED_SAVE.mat\n');
        catch ME_save2
             fprintf('*** 彻底保存失败: %s ***\n', ME_save2.message);
        end
    end
    
    %% 9. *** 新增：为下一次循环做准备 (动态清理) ***
    chunkIndex = chunkIndex + 1; % 移动到下一个分块
    
    % *** 清理本次循环的大型变量以释放内存 ***
    clear results num_AC num_EV ...
          temp_ACs temp_EVs acTable evTable acOpts evOpts ...
          temp_AC_Up temp_AC_Down temp_EV_Up temp_EV_Down ...
          temp_SOC_AC temp_SOC_EV temp_m3;
      
    if runAC
        clear ACs AC_Up AC_Down SOC_AC AC_Up_Individual AC_Down_Individual;
    end
    if runEV
        % [!!! 修改 !!!]
        clear EVs EV_Up EV_Down SOC_EV m3 EV_Up_Individual EV_Down_Individual E_current_EV P_current_EV;
    end
      
end % *** 结束主 while 循环 ***

fprintf('\n====================================================\n');
fprintf('== 所有分块处理完毕。 ==\n');
fprintf('== 结果已保存于 %s 文件夹中。 ==\n', outputDir);
fprintf('====================================================\n');

toc; % 结束总计时