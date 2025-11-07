
% 脚本功能:
% 1. 仅执行空调 (AC) 仿真，移除了所有电动汽车 (EV) 相关的代码。
% 2. 实现 AC 的 "状态化" (Stateful) 仿真，即 SOC 在时间步之间传递。
% 3. 调用 generate_hourly_regulation_signal 生成小时级变化的模拟电网指令。
% 4. 使用 updateACSOC_single 函数根据分解后的电网指令更新 SOC 状态。
% 5. 采用分块 (Chunk) 加载和保存机制，以处理大规模设备数据。

clear; close all; clc;

tic; % 启动一个总计时器

%% 1. 系统初始化 (全局参数)
rng(2023);                                      % 固定随机种子，保证结果可重复
T_total = 24;                                   % 总时长（小时）
dt = 0.05;                                      % 时间分辨率（小时）(例如 5/60)
time_points = 0:dt:T_total;                     % 生成时间序列
T_steps_total = length(time_points);            % 仿真总时间步数
steps_per_hour = round(1/dt);                 % 每小时的时间步数
num_hours = floor(T_steps_total / steps_per_hour); % 总小时数

base_price = 30;                                % 基础电价（元/kWh）
t_adj = 1;                                      % 调节时长（小时）(用于潜力计算)

%% 1.5. 模块运行控制
runAC = true;  % <--- 仅运行 AC 仿真
runEV = false; % <--- 禁用 EV 仿真 (在此脚本中固定为 false)

%% 1.6. 分块范围控制
% --- 在此处配置您想运行的分块范围 ---
startChunk = 1;    % <--- 设置起始块编号 (例如: 1)
endChunk = inf;    % <--- 设置结束块编号 (例如: inf 可运行到文件末尾)
% ---------------------------------

%% 2. 初始化参数 (全局参数)
acFile = 'AC_template1.xlsx';                   % AC 数据文件
outputDir = 'chunk_results_ac_only_responsive'; % 用于存储 AC 响应结果的文件夹
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
             
%% 3. 主分块循环
chunkSize = 10000;  % 每次处理的设备数量
chunkIndex = startChunk;   

while true
    fprintf('\n====================================================\n');
    fprintf('== 准备处理 AC 分块 %d (设备 %d 到 %d) ==\n', ...
             chunkIndex, (chunkIndex-1)*chunkSize + 1, chunkIndex*chunkSize);
    fprintf('====================================================\n');
    
    if chunkIndex > endChunk
        fprintf('已达到设定的结束块 %d，仿真停止。\n', endChunk);
        break; % 退出 while 循环
    end
    
    %% 3.1. 分块读取设备参数
    startRow_Excel = 1 + (chunkIndex - 1) * chunkSize + 1; % +1 是为了跳过标题行
    endRow_Excel = startRow_Excel + chunkSize - 1;
    dataRange = sprintf('%d:%d', startRow_Excel, endRow_Excel);
    
    ACs = []; % 初始化为空
    
    try
        % 设置导入选项以读取特定行
        acOpts = detectImportOptions(acFile, 'ReadVariableNames', true);
        acOpts.DataRange = dataRange; % 指定读取范围
        acTable = readtable(acFile, acOpts);
        
        if ~isempty(acTable)
            % 【修正】: 直接将读取的 table 转换为 struct
            % 这是 ac_ev_simulation_block.m 中的正确逻辑
            ACs = table2struct(acTable); 
        end
        
    catch ME
        fprintf('读取空调文件时出错或到达文件末尾: %s\n', ME.message);
        ACs = []; 
    end

    % 如果读取到的数据为空，则退出循环
    if isempty(ACs)
        fprintf('未在 %s 中找到更多AC数据 (已到达文件末尾或范围 %s 无数据)，处理结束。\n', acFile, dataRange);
        break; % 退出 while 循环
    end
    
    num_AC = length(ACs);
    fprintf('分块 %d: 成功加载 %d 台 AC。\n', chunkIndex, num_AC);

    
    %% 4. 激励响应模块 (仅 AC)
    fprintf('分块 %d: 正在计算 AC 激励响应...\n', chunkIndex);
    p_min = 15; p_max = 50;                         
    p_min_prime = 10; p_max_prime = 40; T_set_max = 3;                

    temp_ACs = ACs; 
    parfor i = 1:num_AC
        % 计算参与概率
        participation = calculateParticipation(temp_ACs(i).p_incentive, base_price);
        % 计算温度灵活性
        [~, ~, deltaT] = incentiveTempAC(...        
            temp_ACs(i).p_incentive, p_min, p_max, p_min_prime, p_max_prime, T_set_max);
        
        temp_ACs(i).ptcp = (rand() < participation);     
        
        % 根据参与意愿调整温度上下限
        if temp_ACs(i).ptcp
            temp_ACs(i).Tmax = temp_ACs(i).Tset + deltaT;
            temp_ACs(i).Tmin = temp_ACs(i).Tset - deltaT;
        else
             temp_ACs(i).Tmax = temp_ACs(i).Tset; 
             temp_ACs(i).Tmin = temp_ACs(i).Tset; 
        end
        
        % 生成环境温度曲线 (T_ja)
        base_temp = temp_ACs(i).Tset + 4*sin(2*pi*time_points/24); 
        temp_range = temp_ACs(i).Tmax - temp_ACs(i).Tmin;
        if abs(temp_range) < 1e-6; temp_range = 0.1; end % 避免范围为0
        
        noise = 0.2 * temp_range * randn(size(time_points));
        
        if temp_ACs(i).ptcp
            % 参与的AC，其T_ja被限制在新的[Tmin, Tmax]内
            temp_ACs(i).T_ja = min(max(base_temp + noise, temp_ACs(i).Tmin), temp_ACs(i).Tmax);
        else
            % 未参与的AC，T_ja仅基于室外温度波动
            temp_ACs(i).T_ja = base_temp + noise; 
        end
    end
    ACs = temp_ACs; % 将并行计算的结果写回主结构体

    %% 5. 预计算模块 (仅 AC)
    fprintf('分块 %d: 正在预计算 AC 参数 (Alpha, Beta, Gamma)...\n', chunkIndex);
    temp_ACs = ACs; 
    parfor i = 1:num_AC
        % 计算理论模型 (3-12) 中的状态转移系数
        [alpha, beta, gamma] = calculateACABC_single(...
            temp_ACs(i).R, temp_ACs(i).C, temp_ACs(i).eta,...
            temp_ACs(i).Tmax, temp_ACs(i).Tmin, temp_ACs(i).Tset, dt);
        temp_ACs(i).alpha = alpha;
        temp_ACs(i).beta = beta;
        temp_ACs(i).gamma = gamma;
    end
    ACs = temp_ACs;
    
    %% 5.5. 生成电网调节指令
    fprintf('分块 %d: 正在生成 %d 小时的电网调节指令...\n', chunkIndex, num_hours);
    % 调用新函数生成小时级指令
    P_grid_command_series = generate_hourly_regulation_signal(T_steps_total, steps_per_hour, num_hours, num_AC);
    
    % 计算参与的 AC 总数 (用于分解指令)
    num_AC_participating = sum([ACs.ptcp]);
    fprintf('  (共 %d / %d 台 AC 参与响应)\n', num_AC_participating, num_AC);


    %% 6. 主时间循环
    %% 6.1 结果预分配
    AC_Up = zeros(length(time_points),1);
    AC_Down = zeros(length(time_points),1);
    SOC_AC = zeros(num_AC, length(time_points));
    AC_Up_Individual = zeros(num_AC, length(time_points));
    AC_Down_Individual = zeros(num_AC, length(time_points));
    
    % 初始化当前状态向量 (即 SOC(t=0))
    % 我们使用从Excel加载的初始SOC值
    CURRENT_SOC_AC = [ACs.SOC]'; % [num_AC x 1] 向量

    %% 6.2 时间步进循环
    fprintf('开始分块 %d 的时间步进仿真 (状态化)...\n', chunkIndex);
    
    for t_idx = 1:T_steps_total
        t = time_points(t_idx);
        
        if mod(t_idx-1, steps_per_hour) == 0 % 每小时打印一次
             fprintf('  分块 %d: 仿真时间 %.2f小时 (%.1f%%)\n', ...
                 chunkIndex, t, 100*t_idx/T_steps_total);
        end
        
        % 获取当前时间步的电网总指令
        P_grid_command_t = P_grid_command_series(t_idx);
        
        % 分解指令到每个参与的设备 (平均分解)
        if num_AC_participating > 0
            Delta_Pj_command_t_step = P_grid_command_t / num_AC_participating;
        else
            Delta_Pj_command_t_step = 0;
        end
        
        % 临时存储当前步的聚合值
        temp_AC_Up_agg = 0;
        temp_AC_Down_agg = 0;
        
        % 临时存储下一时间步的状态
        temp_SOC_AC_for_next_step = zeros(num_AC, 1);
        % 临时存储当前时间步的个体结果 (用于保存)
        temp_AC_Up_Ind = zeros(num_AC, 1);
        temp_AC_Down_Ind = zeros(num_AC, 1);
        temp_SOC_AC_for_saving = zeros(num_AC, 1);
            
        parfor i = 1:num_AC
            
            % 1. 获取当前状态 SOC(t)
            SOC_current_i = CURRENT_SOC_AC(i);
            temp_SOC_AC_for_saving(i) = SOC_current_i; % 记录用于保存的 SOC(t)

            if ACs(i).ptcp
                % 2. 计算当前时间 t 的基线和潜力 (基于 SOC(t))
                P_base_i = ACbaseP_single(ACs(i).T_ja(t_idx), ACs(i).Tset, ACs(i).R, ACs(i).eta);
                
                % 计算潜力
                [DeltaP_plus, DeltaP_minus] = calculateACAdjustmentPotentia(...
                    P_base_i, 2*abs(P_base_i), 0,...
                    ACs(i).alpha, ACs(i).beta, ACs(i).gamma,...
                    SOC_current_i, dt);
                
                % 聚合潜力
                temp_AC_Up_agg = temp_AC_Up_agg + DeltaP_plus;
                temp_AC_Down_agg = temp_AC_Down_agg + DeltaP_minus;
                
                % 存储个体潜力
                temp_AC_Up_Ind(i) = DeltaP_plus;
                temp_AC_Down_Ind(i) = DeltaP_minus;
                
                % 3. 应用新函数和分解后的指令更新状态
                % 仅参与的 AC 接收指令
                Delta_Pj_command = Delta_Pj_command_t_step;
                
                % 检查指令是否超过潜力 (简单削减)
                if Delta_Pj_command > DeltaP_plus
                    Delta_Pj_command = DeltaP_plus;
                elseif Delta_Pj_command < DeltaP_minus
                    Delta_Pj_command = DeltaP_minus;
                end
                
                % 根据指令计算 SOC(t+1)
                SOC_next_i = updateACSOC_single(SOC_current_i, Delta_Pj_command, ...
                    ACs(i).alpha, ACs(i).beta, ACs(i).gamma);
                
                temp_SOC_AC_for_next_step(i) = SOC_next_i; % 存储 SOC(t+1)
            
            else % 如果不参与 (ptcp == false)
                
                % 不参与的设备按基线运行 (Delta_Pj = 0)
                SOC_next_i = updateACSOC_single(SOC_current_i, 0, ...
                    ACs(i).alpha, ACs(i).beta, ACs(i).gamma);
                
                temp_AC_Up_Ind(i) = 0;   % 无潜力
                temp_AC_Down_Ind(i) = 0; % 无潜力
                temp_SOC_AC_for_next_step(i) = SOC_next_i; % 状态仍然按基线更新
            end
        end % 结束 parfor AC
        
        % 4. 存储当前时间步 t 的聚合结果
        AC_Up(t_idx) = temp_AC_Up_agg;
        AC_Down(t_idx) = temp_AC_Down_agg;
        
        % 5. 存储当前时间步 t 的个体结果
        AC_Up_Individual(:, t_idx) = temp_AC_Up_Ind;
        AC_Down_Individual(:, t_idx) = temp_AC_Down_Ind;
        SOC_AC(:, t_idx) = temp_SOC_AC_for_saving; % 保存 SOC(t)
        
        % 6. 更新状态向量用于下一个循环
        CURRENT_SOC_AC = temp_SOC_AC_for_next_step; % CURRENT_SOC_AC 现在是 SOC(t+1)
        
    end % 结束 for t_idx
    fprintf('分块 %d 仿真完成。\n', chunkIndex);

    %% 7. 将结果组装回结构体
    results = struct(); % 创建一个空结构体
    
    results.AC_Up = AC_Up;
    results.AC_Down = AC_Down;
    results.SOC_AC = SOC_AC;
    results.AC_Up_Individual = AC_Up_Individual;
    results.AC_Down_Individual = AC_Down_Individual;
    results.P_grid_command_series = P_grid_command_series; % 保存电网指令
    
    %% 8. 保存当前分块的结果
    outputFileName = fullfile(outputDir, sprintf('results_chunk_ac_%d.mat', chunkIndex));
    try
        fprintf('正在保存分块 %d 的结果到 %s ...\n', chunkIndex, outputFileName);
        save(outputFileName, 'results', '-v7.3'); % 使用 -v7.3 以支持大文件
        fprintf('保存成功。\n');
    catch ME_save
        fprintf('*** 保存文件时出错: %s ***\n', ME_save.message);
        try
            save('LAST_CHUNK_AC_FAILED_SAVE.mat', 'results', '-v7.3');
            fprintf('已保存到 LAST_CHUNK_AC_FAILED_SAVE.mat\n');
        catch ME_save2
             fprintf('*** 彻底保存失败: %s ***\n', ME_save2.message);
        end
    end
    
    %% 9. 为下一次循环做准备 (动态清理)
    chunkIndex = chunkIndex + 1; % 移动到下一个分块
    
    % 清理本次循环的大型变量以释放内存
    clear results num_AC ...
          temp_ACs acTable acOpts ...
          AC_Up AC_Down SOC_AC AC_Up_Individual AC_Down_Individual ...
          CURRENT_SOC_AC P_grid_command_series;
      
end % *** 结束主 while 循环 ***

fprintf('\n====================================================\n');
fprintf('== 所有 AC 分块处理完毕。 ==\n');
fprintf('== 结果已保存于 %s 文件夹中。 ==\n', outputDir);
fprintf('====================================================\n');

toc; % 结束总计时