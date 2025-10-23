clear; close all; clc;
    
    tic; % 开始计时
    
    %% 1. 系统初始化
    rng(2023);                                      % 固定随机种子，保证结果可重复
    T_total = 24;                                   % 总时长（小时）
    dt = 5/60;                                      % 时间分辨率（小时）
    time_points = 0:dt:T_total;                     % 生成时间序列
    base_price = 30;                                % 基础电价（元/kWh）
    t_adj = 60/60;                                      % 调节时长（小时）
    num_time_points = length(time_points);
    
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
    
    %% 4. 激励响应模块
    fprintf('正在计算激励响应 (并行)...\n');
    %% 4.1 参数设定
    p_min = 15; p_max = 50;                         % 原始电价范围
    p_min_prime = 10; p_max_prime = 40; T_set_max = 3; % 调整后电价范围
    
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
        end
        base_temp = temp_ACs(i).Tset + 4*sin(2*pi*time_points/24); 
        temp_range = temp_ACs(i).Tmax - temp_ACs(i).Tmin;
        noise = 0.2 * temp_range * randn(size(time_points));
        temp_ACs(i).T_ja = min(max(base_temp + noise, temp_ACs(i).Tmin), temp_ACs(i).Tmax);
    end
    ACs = temp_ACs;
    clear temp_ACs; 
    
    %% 4.2 EV目标电量调整
    temp_EVs = EVs;
    parfor i = 1:num_EV
        temp_EVs(i).p_incentive = 11;
        participation = calculateParticipation(temp_EVs(i).p_incentive, base_price);
        temp_EVs(i).ptcp = (rand() < participation);
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
    clear temp_EVs; 
    
    %% 5. 预计算模块
    fprintf('正在执行预计算 (并行)...\n');
    %% 5.1 EV基线功率计算
    temp_EVs = EVs;
    parfor i = 1:num_EV
        % 假设 EVbaseP_ChargeUntilFull 返回 T x 1 (列向量)
        temp_EVs(i).P_base_sequence = EVbaseP_ChargeUntilFull(...
            temp_EVs(i).C_EV, temp_EVs(i).eta,...
            temp_EVs(i).E_tar_original, temp_EVs(i).E_in,...
            temp_EVs(i).t_dep, temp_EVs(i).t_in, dt, ... 
            temp_EVs(i).r, temp_EVs(i).p_on, temp_EVs(i).SOC, num_time_points, time_points);
    end
    EVs = temp_EVs;
    clear temp_EVs; 
    
    %% 5.2 AC 预计算
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
    % AC_T_ja 是 M x T 矩阵 (假设 T_ja 是 1xT 行向量)
    AC_T_ja     = cat(1, ACs.T_ja); 
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
    
    % --- (!!! 此处是关键修复 !!!) ---
    % 假设 EVs(i).P_base_sequence 是 T x 1 (列向量)
    % 我们使用 cat(2, ...) 将它们水平堆叠
    % EV_P_base_seq 变为 T x M 矩阵 (T行, M列)
    EV_P_base_seq = cat(2, EVs.P_base_sequence); 
    % --- (修复结束) ---
    
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
    EV_E_exp      = [EVs.E_exp]';
    EV_E_current  = [EVs.E_current]';
    EV_P_current  = [EVs.P_current]';
    
    m3 = zeros(num_EV,1);
    temp_m3 = m3;
    parfor i = 1:num_EV
        if EV_ptcp(i)
            [~, ~, m3_val] = calculateEVABC_single(EV_C_EV(i), EV_eta(i), EV_E_tar_orig(i), EV_E_in(i), EV_t_dep(i), EV_t_in(i), dt, EV_r(i));
            temp_m3(i) = m3_val;
        end
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
    
    % --- 个体结果：创建 matfile 对象 (硬盘映射) ---
    individual_results_file = 'individual_results.mat';
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
        
    catch E
        fprintf('创建 matfile 失败! 请检查硬盘权限或空间。\n');
        rethrow(E);
    end
    fprintf('硬盘预分配完成。\n');
    
    %% 6.2 时间步进循环 (!!! 硬盘存储修改 !!!)
    
    temp_EV_E_exp     = EV_E_exp;
    temp_EV_E_current = EV_E_current;
    temp_EV_P_current = EV_P_current;
            
    for t_idx = 1:num_time_points
        t = time_points(t_idx);
        if mod(t_idx-1, round(1/dt)) == 0 
            fprintf('== 正在处理时间 %.2f/%.2f 小时 (%.1f%%) ==\n', ...
                t, T_total, (t_idx/num_time_points)*100);
        end
        
        temp_EV_Up_Ind   = zeros(num_EV, 1, 'single');
        temp_EV_Down_Ind = zeros(num_EV, 1, 'single');
        temp_SOC_EV_Ind  = zeros(num_EV, 1, 'single');
        temp_AC_Up_Ind   = zeros(num_AC, 1, 'single');
        temp_AC_Down_Ind = zeros(num_AC, 1, 'single');
        temp_SOC_AC_Ind  = zeros(num_AC, 1, 'single');
    
        %% 电动汽车集群处理 - (使用解构后的数组)
        temp_EV_Up = 0;   
        temp_EV_Down = 0; 
        
        % --- (!!! 此处是关键修复 !!!) ---
        % EV_P_base_seq 是 T x M 矩阵
        % 我们获取第 t_idx 行，这是一个 1 x M 的行向量
        % 包含了所有 M 台设备在 t_idx 刻的基线功率
        current_P_base_row_vector = EV_P_base_seq(t_idx, :);
        % --- (修复结束) ---
    
        parfor i = 1:num_EV
            if EV_ptcp(i)
                online = (t >= EV_t_in(i)) && (t < EV_t_dep(i));
                if online
                    % --- (!!! 此处是关键修复 !!!) ---
                    % 从 1 x M 的行向量中获取第 i 个元素
                    P_base_i = current_P_base_row_vector(i);
                    % --- (修复结束) ---
                    
                    [E_exp_i, E_current_i, P_current_i, SOC_i] = ...
                        calculateEVS_single(m3(i), EV_E_exp(i), EV_E_tar_orig(i),...
                        EV_eta(i), EV_E_current(i), EV_P_current(i), EV_C_EV(i),...
                        EV_r(i), EV_p_on(i), dt, EV_t_dep(i), t);
                    
                    [DeltaP_plus, DeltaP_minus] = calculateEVAdjustmentPotentia_new(...
                        EV_E_reg_min(i), EV_E_reg_max(i), ...
                        E_current_i, EV_t_dep(i), t, ...
                        EV_p_on(i), P_base_i, EV_eta(i), t_adj);
                    
                    temp_EV_Up = temp_EV_Up + DeltaP_plus;
                    temp_EV_Down = temp_EV_Down + DeltaP_minus;
                    
                    temp_EV_Up_Ind(i)   = DeltaP_plus;
                    temp_EV_Down_Ind(i) = DeltaP_minus;
                    temp_SOC_EV_Ind(i)  = SOC_i;
    
                    temp_EV_E_exp(i)     = E_exp_i;
                    temp_EV_E_current(i) = E_current_i;
                    temp_EV_P_current(i) = P_current_i;
                else
                    temp_SOC_EV_Ind(i) = EV_E_current(i) / EV_C_EV(i);
                end
            else 
                temp_SOC_EV_Ind(i) = EV_E_current(i) / EV_C_EV(i);
            end
        end % 结束 parfor EV
        
        EV_Up(t_idx) = temp_EV_Up;
        EV_Down(t_idx) = temp_EV_Down;
        
        m.EV_Up_Individual(:, t_idx) = temp_EV_Up_Ind;
        m.EV_Down_Individual(:, t_idx) = temp_EV_Down_Ind;
        m.SOC_EV(:, t_idx) = temp_SOC_EV_Ind;
        
        EV_E_exp     = temp_EV_E_exp;
        EV_E_current = temp_EV_E_current;
        EV_P_current = temp_EV_P_current;
        
        %% 空调分析 - (AC 逻辑保持不变)
        temp_AC_Up = 0;   
        temp_AC_Down = 0; 
        
        % AC_T_ja 是 M x T 矩阵, 提取第 t_idx 列 (M x 1 向量)
        current_T_ja = AC_T_ja(:, t_idx);
        
        parfor i = 1:num_AC
            if AC_ptcp(i)
                % 从 M x 1 向量中提取第 i 个元素
                P_base_i = ACbaseP_single(current_T_ja(i), AC_Tset(i), AC_R(i), AC_eta(i));
                
                SOC_i = calculateACS_single(current_T_ja(i), AC_Tmax(i), AC_Tmin(i));
                
                [DeltaP_plus, DeltaP_minus] = calculateACAdjustmentPotentia(...
                    P_base_i, 2*abs(P_base_i), 0,...
                    AC_alpha(i), AC_beta(i), AC_gamma(i),...
                    SOC_i, dt);
                
                temp_AC_Up = temp_AC_Up + DeltaP_plus;
                temp_AC_Down = temp_AC_Down + DeltaP_minus;
                
                temp_AC_Up_Ind(i)   = DeltaP_plus;
                temp_AC_Down_Ind(i) = DeltaP_minus;
                temp_SOC_AC_Ind(i)  = SOC_i;
            else
                temp_SOC_AC_Ind(i) = 0.5; % 或者 NaN
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
        'm3', m3...
        );
    
    % 将聚合结果保存到单独的文件
    aggregate_results_file = 'aggregate_results.mat';
    save(aggregate_results_file, 'results');
    
    fprintf('仿真完成。\n');
    fprintf('======================================================\n');
    fprintf('聚合结果 (AC_Up, EV_Up...) 已保存到: %s\n', aggregate_results_file);
    fprintf('全量个体数据 (SOC_AC, EV_Up_Individual...) 已保存到: %s\n', individual_results_file);
    fprintf('======================================================\n');
    toc; % 结束计时