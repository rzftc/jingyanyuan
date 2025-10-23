%% 全功能虚拟电厂调节潜力分析系统 (*** 分块处理版 ***)
clear; close all; clc;

%% 1. 系统初始化 (全局参数)
rng(2023);                                      % 固定随机种子，保证结果可重复
T_total = 24;                                   % 总时长（小时）
dt = 0.05;                                      % 时间分辨率（小时）
time_points = 0:dt:T_total;                     % 生成时间序列
base_price = 30;                                % 基础电价（元/kWh）
t_adj = 1;                                      % 调节时长（小时），从 new main file 引入

%% 2. 初始化参数 (全局参数)
evFile = '2EV_residential.xlsx';                 % 电动汽车数据文件
acFile = 'AC_template1.xlsx';       
outputDir = 'chunk_results';                   % *** 新增：用于存储结果的文件夹 ***
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
             
%% 3. *** 新增：主分块循环 ***
chunkSize = 1000;  % *** 每次处理1000台设备 ***
chunkIndex = 1;   % *** 当前分块索引 ***

while true
    fprintf('\n====================================================\n');
    fprintf('== 开始处理分块 %d (设备 %d 到 %d) ==\n', ...
             chunkIndex, (chunkIndex-1)*chunkSize + 1, chunkIndex*chunkSize);
    fprintf('====================================================\n');
    
    %% 3.1. 分块读取设备参数 (*** 替换原始的 Section 3 ***)
    
    % --- 计算要读取的 Excel 行范围 ---
    % 假设第1行是标题，数据从第2行开始
    startRow_Excel = 1 + (chunkIndex - 1) * chunkSize + 1; % +1 是为了跳过标题行
    endRow_Excel = startRow_Excel + chunkSize - 1;
    dataRange = sprintf('%d:%d', startRow_Excel, endRow_Excel);
    
    % --- 读取 AC 数据 ---
    ACs = []; % 清空上一个分块的数据
    try
        fprintf('正在从 %s 读取行 %s...\n', acFile, dataRange);
        acOpts = detectImportOptions(acFile, 'ReadVariableNames', true);
        acOpts.DataRange = dataRange; % 指定读取范围
        
        acTable = readtable(acFile, acOpts);
        
        if isempty(acTable)
            fprintf('未在 %s 中找到更多AC数据，处理结束。\n', acFile);
            break; % 退出 while 循环
        end
        
        % 将表格转换为结构体数组
        % *要求Excel列标题与代码中的结构体字段名完全匹配*
        ACs = table2struct(acTable);
        
    catch ME
        fprintf('读取空调文件时出错或到达文件末尾: %s\n', ME.message);
        break; % 退出 while 循环
    end

    % --- 读取 EV 数据 ---
    EVs = []; % 清空上一个分块的数据
    try
        fprintf('正在从 %s 读取行 %s...\n', evFile, dataRange);
        evOpts = detectImportOptions(evFile, 'ReadVariableNames', true);
        evOpts.DataRange = dataRange; % 指定读取范围
        
        evTable = readtable(evFile, evOpts);
        
        if isempty(evTable)
            fprintf('未在 %s 中找到更多EV数据，处理结束。\n', evFile);
            break; % 退出 while 循环
        end
        
        % 将表格转换为结构体数组
        % *要求Excel列标题与代码中的结构体字段名完全匹配*
        EVs = table2struct(evTable);
        
    catch ME
        fprintf('读取EV文件时出错或到达文件末尾: %s\n', ME.message);
        break; % 退出 while 循环
    end

    % --- 检查数据 ---
    if isempty(ACs) || isempty(EVs)
         fprintf('空调或电动汽车数据为空 (可能一个文件比另一个短)，处理结束。\n');
         break;
    end
    
    num_AC = length(ACs);                          
    num_EV = length(EVs);
    fprintf('分块 %d: 成功加载 %d 台 AC 和 %d 台 EV。\n', chunkIndex, num_AC, num_EV);

    
    % =================================================================
    % == 以下是您原始脚本的核心逻辑 (Section 4, 5, 6) ==
    % == 这些代码现在在 while 循环内部，对当前分块的数据进行操作 ==
    % =================================================================

    %% 4. 激励响应模块
    %% 4.1 参数设定
    p_min = 15; p_max = 50;                         % 原始电价范围
    p_min_prime = 10; p_max_prime = 40; T_set_max = 3;                % 调整后电价范围

    %% AC温度设定调整 - 修改为使用临时结构体数组
    temp_ACs = ACs;  % 创建临时结构体数组
    parfor i = 1:num_AC
        participation = calculateParticipation(temp_ACs(i).p_incentive, base_price);
        [~, ~, deltaT] = incentiveTempAC(...        
            temp_ACs(i).p_incentive, p_min, p_max, p_min_prime, p_max_prime, T_set_max);
        
        temp_ACs(i).ptcp = (rand() < participation);     
        
        if temp_ACs(i).ptcp
            temp_ACs(i).Tmax = temp_ACs(i).Tset + deltaT;
            temp_ACs(i).Tmin = temp_ACs(i).Tset - deltaT;
        end
        
        base_temp = temp_ACs(i).Tset + 4*sin(2*pi*time_points/24); 
        
        % *** 关键修复：确保未参与的AC也有Tmax/Tmin，否则下一行出错 ***
        if ~temp_ACs(i).ptcp
             temp_ACs(i).Tmax = temp_ACs(i).Tset;
             temp_ACs(i).Tmin = temp_ACs(i).Tset;
        end
        
        temp_range = temp_ACs(i).Tmax - temp_ACs(i).Tmin;
        
        % *** 关键修复：如果未参与，temp_range可能为0，导致noise为0 ***
        if temp_range == 0 
            temp_range = 0.1; % 给予一个极小的范围，以防randn变为0
        end
        
        noise = 0.2 * temp_range * randn(size(time_points));
        
        % *** 关键修复：确保T_ja的边界对所有AC都有效 ***
        if temp_ACs(i).ptcp
            temp_ACs(i).T_ja = min(max(base_temp + noise, temp_ACs(i).Tmin), temp_ACs(i).Tmax);
        else
            % 未参与的AC，T_ja也应在基准附近波动，但不受Tmin/Tmax强约束
            % (或者您可以让它严格等于base_temp)
            temp_ACs(i).T_ja = base_temp + noise; 
        end
    end
    ACs = temp_ACs;  % 更新原始结构体数组

    %% 4.2 EV目标电量调整 - *** 更新为 new main file 逻辑 ***
    temp_EVs = EVs;  % 创建临时结构体数组
    parfor i = 1:num_EV
        temp_EVs(i).p_incentive = 11;                     % 设置激励电价 (保留原文件逻辑)
        participation = calculateParticipation(temp_EVs(i).p_incentive, base_price);
        ptcp_result = (rand() < participation);           % 生成参与决策
        temp_EVs(i).ptcp = ptcp_result;
        
        % 存储原始目标电量
        temp_EVs(i).E_tar_original = temp_EVs(i).E_tar;
        
        E_flex_max = 0.2 * temp_EVs(i).C_EV; % (来自 new main file)
        [deltaE_up, deltaE_down] = incentiveTempEV_updown(temp_EVs(i).p_incentive, p_min, p_max, p_min_prime, p_max_prime, E_flex_max); % (来自 new main file)
        
        if temp_EVs(i).ptcp
            temp_EVs(i).E_reg_min = temp_EVs(i).E_tar_original - deltaE_down; % (来自 new main file)
             if temp_EVs(i).E_reg_min <= temp_EVs(i).E_in
                 temp_EVs(i).E_reg_min = temp_EVs(i).E_in; % (来自 new main file)
             end
            temp_EVs(i).E_reg_max = temp_EVs(i).E_tar_original + deltaE_up; % (来自 new main file)
             if temp_EVs(i).E_reg_max >= temp_EVs(i).C_EV
                 temp_EVs(i).E_reg_max = temp_EVs(i).C_EV; % (来自 new main file)
             end
        else
            temp_EVs(i).E_reg_min = temp_EVs(i).E_tar_original; % (来自 new main file)
            temp_EVs(i).E_reg_max = temp_EVs(i).E_tar_original; % (来自 new main file)
        end
        
        % *** 注意：原文件修改 E_tar 的逻辑被替换为 E_reg_min 和 E_reg_max ***
    end
    EVs = temp_EVs;  % 更新原始结构体数组

    %% 5. 预计算模块
    %% 5.1 EV基线功率计算 - *** 更新为 new main file 逻辑 ***
    H = T_total;                                     % 预测时域（小时） (使用原T_total)
    H_steps = H / dt;                                % 转换为时间步数
    temp_EVs = EVs;  % 创建临时结构体数组
    num_time_points_for_baseline = length(time_points); % (来自 new main file)

    parfor i = 1:num_EV
        % (来自 new main file)
        temp_EVs(i).P_base_sequence = EVbaseP_ChargeUntilFull(...
            temp_EVs(i).C_EV, temp_EVs(i).eta,...
            temp_EVs(i).E_tar_original, temp_EVs(i).E_in,... % 使用 E_tar_original
            temp_EVs(i).t_dep, temp_EVs(i).t_in, dt, ... 
            temp_EVs(i).r, temp_EVs(i).p_on, temp_EVs(i).SOC, num_time_points_for_baseline, time_points); % 传入 time_points
    end
    EVs = temp_EVs;  % 更新原始结构体数组

    temp_ACs = ACs;  % 创建临时结构体数组
    parfor i = 1:num_AC
        % *** 修复：确保未参与的AC的Tmax/Tmin已在 4.1 中设置 ***
        [alpha, beta, gamma] = calculateACABC_single(...
            temp_ACs(i).R, temp_ACs(i).C, temp_ACs(i).eta,...
            temp_ACs(i).Tmax, temp_ACs(i).Tmin, temp_ACs(i).Tset, dt);
        temp_ACs(i).alpha = alpha;
        temp_ACs(i).beta = beta;
        temp_ACs(i).gamma = gamma;
    end
    ACs = temp_ACs;  % 更新原始结构体数组

    %% 6. 主时间循环
    %% 6.1 结果预分配
    AC_Up = zeros(length(time_points),1);
    AC_Down = zeros(length(time_points),1);
    EV_Up = zeros(length(time_points),1);
    EV_Down = zeros(length(time_points),1);
    SOC_AC = zeros(num_AC, length(time_points));
    SOC_EV = zeros(num_EV, length(time_points));
    m3 = zeros(num_EV,1);

    % 新增：存储每台设备的调节能力
    AC_Up_Individual = zeros(num_AC, length(time_points));
    AC_Down_Individual = zeros(num_AC, length(time_points));
    EV_Up_Individual = zeros(num_EV, length(time_points));
    EV_Down_Individual = zeros(num_EV, length(time_points));
    
    % *** 关键修复：为EVs结构体预先添加字段，避免parfor中出错 ***
    % 检查 'E_exp' 和 'E_current' 字段是否存在，如果不存在则初始化
    if ~isfield(EVs, 'E_exp')
        [EVs.E_exp] = deal(0); % 或者使用合适的初始值
    end
    if ~isfield(EVs, 'E_current')
         % E_current 似乎应该初始化为 E_in
         for i=1:num_EV
             EVs(i).E_current = EVs(i).E_in;
         end
    end
    if ~isfield(EVs, 'P_current')
        [EVs.P_current] = deal(0);
    end

    %% 6.2 时间步进循环
    fprintf('开始分块 %d 的时间步进仿真...\n', chunkIndex);
    for t_idx = 1:length(time_points)
        t = time_points(t_idx);
        
        % *** 修改：减少打印频率，避免刷屏 ***
        if mod(t_idx, round(length(time_points)/20)) == 1 || t_idx == length(time_points)
             fprintf('  分块 %d: 仿真时间 %.2f小时 (%.1f%%)\n', ...
                 chunkIndex, t, 100*t_idx/length(time_points));
        end
    
        %% 电动汽车集群处理 - *** 更新为 new main file 逻辑 ***
        temp_EVs = EVs;  % 创建临时结构体数组
        temp_EV_Up = 0;
        temp_EV_Down = 0;
        temp_SOC_EV = zeros(num_EV,1);
        temp_m3 = m3;
        
        parfor i = 1:num_EV
            online = (t >= temp_EVs(i).t_in) && (t < temp_EVs(i).t_dep); % (使用 < t_dep 匹配 new main file)
            
            if temp_EVs(i).ptcp % (只处理参与的EV, 来自 new main file)
                
                if online
                    %% 基线功率获取
                    temp_EVs(i).P_base = temp_EVs(i).P_base_sequence(t_idx);
                    
                    %% 状态更新
                    [~, ~, m3_val] = calculateEVABC_single(...
                        temp_EVs(i).C_EV, temp_EVs(i).eta,...
                        temp_EVs(i).E_tar_original, temp_EVs(i).E_in,... % 使用 E_tar_original
                        temp_EVs(i).t_dep, temp_EVs(i).t_in, dt, temp_EVs(i).r);
                    temp_m3(i) = m3_val;
                    
                    % (更新 calculateEVS_single 调用, 来自 new main file)
                    % *** E_exp, E_current, P_current 在parfor中被读取和写入 ***
                    % *** 它们必须在循环开始前存在于 temp_EVs(i) 中 ***
                    [temp_EVs(i).E_exp, temp_EVs(i).E_current,...
                        temp_EVs(i).P_current, temp_EVs(i).SOC] = ...
                        calculateEVS_single(m3_val, temp_EVs(i).E_exp, temp_EVs(i).E_tar_original,... % 使用 E_tar_original
                        temp_EVs(i).eta, temp_EVs(i).E_current,...
                        temp_EVs(i).P_current, temp_EVs(i).C_EV,...
                        temp_EVs(i).r, temp_EVs(i).p_on, dt, temp_EVs(i).t_dep, t); % 传入 t_dep 和 t
                    
                    %% 调节潜力计算
                    % (更新为 calculateEVAdjustmentPotentia_new, 来自 new main file)
                    [DeltaP_plus, DeltaP_minus] = calculateEVAdjustmentPotentia_new(...
                        temp_EVs(i).E_reg_min, temp_EVs(i).E_reg_max, ...
                        temp_EVs(i).E_current, temp_EVs(i).t_dep, t, ... % 传入 t
                        temp_EVs(i).p_on, temp_EVs(i).P_base, temp_EVs(i).eta, t_adj); % 传入 t_adj
                    
                    temp_EV_Up = temp_EV_Up + DeltaP_plus;
                    temp_EV_Down = temp_EV_Down + DeltaP_minus;
                    
                    % 新增：存储单台EV调节能力
                    EV_Up_Individual(i, t_idx) = DeltaP_plus;
                    EV_Down_Individual(i, t_idx) = DeltaP_minus;
                    
                    temp_SOC_EV(i) = temp_EVs(i).SOC;
                else
                    % (保持SOC, 来自 new main file)
                    if t_idx > 1 
                        temp_SOC_EV(i) = SOC_EV(i, t_idx-1); % 从上一时刻继承
                    else 
                        temp_SOC_EV(i) = temp_EVs(i).SOC; % 初始SOC
                    end
                    EV_Up_Individual(i, t_idx)   = 0; % 离线无潜力
                    EV_Down_Individual(i, t_idx) = 0; % 离线无潜力
                    
                    % *** 关键：确保离线时 E_current 和 P_current 不变 ***
                    % (parfor 要求所有分支都分配 temp_EVs(i).E_current 等)
                    if t_idx > 1
                        temp_EVs(i).E_current = EVs(i).E_current; % 保持上一时刻的值
                        temp_EVs(i).P_current = EVs(i).P_current; % 保持上一时刻的值
                        temp_EVs(i).E_exp     = EVs(i).E_exp;     % 保持上一时刻的值
                    end
                    % 初始SOC在循环开始前已设置
                    
                end
            
            else % 如果不参与 (ptcp == false)
                 if t_idx > 1 
                    temp_SOC_EV(i) = SOC_EV(i, t_idx-1);
                 else 
                    temp_SOC_EV(i) = temp_EVs(i).SOC; 
                 end
                 EV_Up_Individual(i, t_idx)   = 0; % 不参与无潜力
                 EV_Down_Individual(i, t_idx) = 0; % 不参与无潜力
                 
                 % *** 关键：确保不参与时 E_current 等也不变 ***
                 if t_idx > 1
                    temp_EVs(i).E_current = EVs(i).E_current; 
                    temp_EVs(i).P_current = EVs(i).P_current;
                    temp_EVs(i).E_exp     = EVs(i).E_exp;
                 end
            end
            
        end
        
        EVs = temp_EVs;  % 更新原始结构体数组 (*** 这是parfor后的关键更新 ***)
        m3 = temp_m3;
        EV_Up(t_idx) = temp_EV_Up;
        EV_Down(t_idx) = temp_EV_Down;
        SOC_EV(:, t_idx) = temp_SOC_EV;
        
        %% 空调分析 - (保持原文件不变)
        temp_ACs = ACs;  % 创建临时结构体数组
        temp_AC_Up = 0;
        temp_AC_Down = 0;
        temp_SOC_AC = zeros(num_AC,1);
        
        parfor i = 1:num_AC
            if temp_ACs(i).ptcp
                % 计算基线功率
                temp_ACs(i).P_base = ACbaseP_single(...
                    temp_ACs(i).T_ja(t_idx), temp_ACs(i).Tset, temp_ACs(i).R, temp_ACs(i).eta);
                
                % 更新SOC
                temp_ACs(i).SOC = calculateACS_single(...
                    temp_ACs(i).T_ja(t_idx), temp_ACs(i).Tmax, temp_ACs(i).Tmin);
                
                % 计算调节潜力
                [DeltaP_plus, DeltaP_minus] = calculateACAdjustmentPotentia(...
                    temp_ACs(i).P_base, 2*abs(temp_ACs(i).P_base), 0,...
                    temp_ACs(i).alpha, temp_ACs(i).beta, temp_ACs(i).gamma,...
                    temp_ACs(i).SOC, dt);
                
                temp_AC_Up = temp_AC_Up + DeltaP_plus;
                temp_AC_Down = temp_AC_Down + DeltaP_minus;
                
                % 新增：存储单台AC调节能力
                AC_Up_Individual(i, t_idx) = DeltaP_plus;
                AC_Down_Individual(i, t_idx) = DeltaP_minus;
                
                temp_SOC_AC(i) = temp_ACs(i).SOC;
            else
                % *** 关键修复：确保不参与的AC的SOC也被记录 ***
                if t_idx > 1
                    temp_SOC_AC(i) = SOC_AC(i, t_idx-1); % 假设SOC不变
                else
                    % 计算一个初始SOC (或者在 5.2 预计算)
                    temp_SOC_AC(i) = calculateACS_single(...
                         temp_ACs(i).T_ja(1), temp_ACs(i).Tmax, temp_ACs(i).Tmin);
                end
                AC_Up_Individual(i, t_idx) = 0;
                AC_Down_Individual(i, t_idx) = 0;
            end
        end
        
        ACs = temp_ACs;  % 更新原始结构体数组
        AC_Up(t_idx) = temp_AC_Up;
        AC_Down(t_idx) = temp_AC_Down;
        SOC_AC(:, t_idx) = temp_SOC_AC;
    end
    fprintf('分块 %d 仿真完成。\n', chunkIndex);

    %% 将结果组装回结构体 - (保持原文件不变)
    results = struct(...
        'AC_Up', AC_Up,...
        'AC_Down', AC_Down,...
        'EV_Up', EV_Up,...
        'EV_Down', EV_Down,...
        'SOC_AC', SOC_AC,...
        'SOC_EV', SOC_EV,...
        'm3', m3,...
        'AC_Up_Individual', AC_Up_Individual,...  % 新增：每台AC上调能力
        'AC_Down_Individual', AC_Down_Individual,... % 新增：每台AC下调能力
        'EV_Up_Individual', EV_Up_Individual,... % 新增：每台EV上调能力
        'EV_Down_Individual', EV_Down_Individual); % 新增：每台EV下调能力
    
    
    %% 7. *** 新增：保存当前分块的结果 ***
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
    
    %% 8. *** 新增：为下一次循环做准备 ***
    chunkIndex = chunkIndex + 1; % 移动到下一个分块
    
    % *** 清理本次循环的大型变量以释放内存 ***
    clear ACs EVs results num_AC num_EV AC_Up AC_Down EV_Up EV_Down ...
          SOC_AC SOC_EV m3 AC_Up_Individual AC_Down_Individual ...
          EV_Up_Individual EV_Down_Individual temp_ACs temp_EVs ...
          acTable evTable acOpts evOpts ...
          temp_AC_Up temp_AC_Down temp_EV_Up temp_EV_Down ...
          temp_SOC_AC temp_SOC_EV temp_m3;
      
end % *** 结束主 while 循环 ***

fprintf('\n====================================================\n');
fprintf('== 所有分块处理完毕。 ==\n');
fprintf('== 结果已保存于 %s 文件夹中。 ==\n', outputDir);
fprintf('====================================================\n');