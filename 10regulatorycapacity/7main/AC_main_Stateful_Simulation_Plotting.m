function AC_main_Stateful_Simulation_Plotting()
    % 主函数：实现空调状态一致性聚合仿真与绘图
    %
    % 1. 基于 AC_main.m 和 ac_simulation_block.m 的逻辑重构。
    % 2. 实现大论文 2.4.1 节的聚合模型 (式 2-37) 和指令分解 (式 2-35)。
    % 3. 严格遵循流程图（image_912548.png）的仿真逻辑。
    % 4. [V2 绘图] 绘制用户要求的 SOC 对比图和功率跟踪对比图 (图例简化)。
    % 5. [V3 用户修改] 绘制所有单体SOC曲线为不同颜色。
    % 6. [V4 用户修改] 新增基于SOC历史反推的室内温度曲线图。
    
    clear; close all; clc;
    tic; 

    %% 1. 系统初始化
    rng(2023, 'Threefry');
    T_total = 24; dt = 5/60; % 5分钟步长
    time_points = 0:dt:T_total;
    T_steps_total = length(time_points);
    steps_per_hour = round(1/dt);
    num_hours = floor(T_steps_total / steps_per_hour);
    base_price = 30;

    %% 2. 初始化 AC 参数
    acFile = 'AC_template.xlsx'; 
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
        if ~isfield(ACs(i), 'p_incentive')
            ACs(i).p_incentive = round(60*rand(), 1);
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

        base_ambient_temp = temp_ACs(i).Tset_original + 4*sin(2*pi*time_points/24);
        actual_temp_range = temp_ACs(i).Tmax - temp_ACs(i).Tmin;
        if abs(actual_temp_range) < 1e-6; actual_temp_range = 0.1; end
        noise = 0.2 * actual_temp_range * randn(size(time_points));
        temp_ACs(i).T_ja = min(max(base_ambient_temp + noise, temp_ACs(i).Tmin), temp_ACs(i).Tmax);

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

    fprintf('  Step 5: 开始 %d 步的状态化仿真...\n', T_steps_total);
    
    for t_idx = 1:T_steps_total
        
        % 1. 获取当前聚合状态 SOC(t)
        SOC_agg_t = mean(CURRENT_SOC_AC, 'omitnan');
        Agg_SOC_History(t_idx) = SOC_agg_t;
        Individual_SOC_History(t_idx, :) = CURRENT_SOC_AC'; % 存储当前步的单体SOC
        
        % 2. 获取电网指令 ΔP_S (流程图 步骤4)
        Delta_P_S_command = P_grid_command_series(t_idx);
        Agg_P_Command_History(t_idx) = Delta_P_S_command;
        
        % 3. 预测目标聚合 SOC(t+1) (流程图 步骤5)
        SOC_target_next = AggParams.A * SOC_agg_t + AggParams.B * Delta_P_S_command + AggParams.C;
        SOC_target_next = max(0, min(1, SOC_target_next)); % 约束

        % 4. 临时变量 (用于 parfor)
        temp_AC_Up_agg = 0;
        temp_AC_Down_agg = 0;
        temp_SOC_for_next_step = zeros(num_AC_participating, 1);
        temp_P_achieved_this_step = zeros(num_AC_participating, 1);
        
        % 5. 【核心】状态转移 (parfor)
        parfor i = 1:num_AC_participating
            
            ac_i = ACs_participating(i);
            soc_current_i = CURRENT_SOC_AC(i); % 获取 SOC(t)
            
            % A. 计算当前物理潜力 (用于记录和约束)
            [P_plus, P_minus] = calculateACAdjustmentPotentia(...
                0, 1e6, -1e6, ... 
                ac_i.alpha, ac_i.beta, ac_i.gamma,...
                soc_current_i, dt);
            
            temp_AC_Up_agg = temp_AC_Up_agg + P_plus;
            temp_AC_Down_agg = temp_AC_Down_agg + P_minus;
            
            % B. 反解理论功率 ΔP_j (流程图 步骤6)
            delta_Pj_theory = 0;
            if abs(ac_i.beta) > 1e-9
                delta_Pj_theory = (SOC_target_next - ac_i.alpha * soc_current_i - ac_i.gamma) / ac_i.beta;
            end
            
            % C. 裁剪指令至物理可行域 (流程图 步骤4 "限幅")
            delta_Pj_clipped = max(P_minus, min(P_plus, delta_Pj_theory));
            
            % D. 更新状态 (实现 式 2-10)
            soc_next_i = updateACSOC_single(soc_current_i, delta_Pj_clipped, ...
                ac_i.alpha, ac_i.beta, ac_i.gamma);
                
            temp_SOC_for_next_step(i) = soc_next_i;
            temp_P_achieved_this_step(i) = delta_Pj_clipped;
        end
        
        % 6. 存储当前时间步 t 的聚合潜力
        Agg_P_Potential_Up_History(t_idx) = temp_AC_Up_agg;
        Agg_P_Potential_Down_History(t_idx) = temp_AC_Down_agg;

        % 7. 存储实际响应功率
        Agg_P_Achieved_History(t_idx) = sum(temp_P_achieved_this_step);
        
        % 8. 更新状态向量用于下一循环
        CURRENT_SOC_AC = temp_SOC_for_next_step; % CURRENT_SOC_AC 现在是 SOC(t+1)
        
    end % 结束 t_idx 循环
    
    fprintf('  Step 5: 仿真完成。\n');

    % --- [新增] 步骤 5.5: 反推室内温度 ---
    fprintf('  Step 5.5: 正在反推室内温度历史...\n');
    % 从参与的空调中获取 Tmax 和 Tmin 向量 (N_p x 1)
    Tmax_vec_p = [ACs_participating.Tmax]'; 
    Tmin_vec_p = [ACs_participating.Tmin]'; 
    TRange_vec_p = Tmax_vec_p - Tmin_vec_p;
    
    % 确保范围不为0，避免除零（尽管 SOC 应该始终在 0-1）
    TRange_vec_p(abs(TRange_vec_p) < 1e-6) = 1e-6; 
    
    % 将 Tmax 和 TRange 复制为 (T_steps_total x N_p) 矩阵
    Tmax_matrix_p = repmat(Tmax_vec_p', T_steps_total, 1); 
    Tmin_matrix_p = repmat(Tmin_vec_p', T_steps_total, 1); % (新增，用于绘图)
    TRange_matrix_p = repmat(TRange_vec_p', T_steps_total, 1);
    
    % Individual_SOC_History 是 (T_steps_total x N_p)
    % 反推公式: T_j(t) = Tmax_j - SOC_j(t) * (Tmax_j - Tmin_j)
    Individual_Temp_History = Tmax_matrix_p - Individual_SOC_History .* TRange_matrix_p;
    fprintf('  Step 5.5: 温度反推完成。\n');
    % --- [新增] 结束 ---

    
    %% 6. 绘图 (实现用户要求)
    fprintf('Step 6: 正在生成对比图...\n');
    
    % --- 图 1: 各空调响应功率之和 vs 电网调节指令 (您的要求 2) ---
    figure('Name', '功率跟踪对比 (理论分解)', 'Position', [100 100 1000 450]);
    ax1 = axes;
    hold(ax1, 'on');
    
    plot(ax1, time_points, Agg_P_Command_History, 'k:', 'LineWidth', 2.5, ...
        'DisplayName', '电网调节指令 (ΔP_s)');
    plot(ax1, time_points, Agg_P_Achieved_History, 'r-', 'LineWidth', 1.5, ...
        'DisplayName', '各空调的响应功率之和 (ΣΔP_j)');
    
    hold(ax1, 'off');
    xlabel(ax1, '时间 (小时)', 'FontSize', 12);
    ylabel(ax1, '聚合功率 (kW)', 'FontSize', 12);
    title(ax1, '图1：空调聚合响应功率 vs 电网指令', 'FontSize', 14);
    legend(ax1, 'show', 'Location', 'best');
    set(ax1, 'FontSize', 11);
    xticks(ax1, [0, 6, 12, 18, 24]);
    xticklabels(ax1, {'00:00', '06:00', '12:00', '18:00', '24:00'});
    xlim(ax1, [0, 24]);
    grid(ax1, 'on');

    
    % --- 图 2: 单体空调SOC vs 聚合模型SOC (*** 按您的要求修改 ***) ---
    figure('Name', 'SOC状态对比 (多曲线)', 'Position', [100 550 1000 450]);
    ax2 = axes;
    hold(ax2, 'on');
    
    % 1. 绘制所有单体SOC (自动分配不同颜色)
    h_individual = plot(ax2, time_points, Individual_SOC_History, 'LineWidth', 0.5);
    
    % 2. 绘制聚合平均SOC (黑色虚线)
    h_agg = plot(ax2, time_points, Agg_SOC_History, 'k--', 'LineWidth', 3, ...
        'DisplayName', '空调聚合模型的SOC (均值)');
    
    % 3. 图例管理
    if num_AC_participating > 0
        set(h_individual(1), 'DisplayName', '单体空调的SOC');
    end
    if num_AC_participating > 1
         set(h_individual(2:end), 'HandleVisibility', 'off');
    end
    
    hold(ax2, 'off');
    
    % 4. 格式化
    xlabel(ax2, '时间 (小时)', 'FontSize', 12);
    ylabel(ax2, '空调SOC', 'FontSize', 12);
    title(ax2, '图2：单体空调SOC 与 聚合模型SOC 对比 (多曲线)', 'FontSize', 14);
    
    % 5. 创建图例
    if num_AC_participating > 0
        legend(ax2, [h_agg, h_individual(1)], 'Location', 'best', 'FontSize', 11);
    else
        legend(ax2, h_agg, 'Location', 'best', 'FontSize', 11);
    end
    
    set(ax2, 'FontSize', 11);
    xticks(ax2, [0, 6, 12, 18, 24]);
    xticklabels(ax2, {'00:00', '06:00', '12:00', '18:00', '24:00'});
    xlim(ax2, [0, 24]);
    ylim(ax2, [-0.1, 1.1]);
    grid(ax2, 'on');
    % --- 修改结束 ---

    % --- [新增] 图 3: 单体空调室内温度 (反推) ---
    figure('Name', '室内温度变化 (反推)', 'Position', [100 300 1000 450]);
    ax3 = axes;
    hold(ax3, 'on');
    
    % 1. 绘制所有单体温度 (自动使用不同颜色)
    h_temp_individual = plot(ax3, time_points, Individual_Temp_History, 'LineWidth', 0.5);
    
    % 2. 绘制 Tmax 和 Tmin 的平均边界 (作为参考)
    h_tmax_avg = plot(ax3, time_points, mean(Tmax_matrix_p, 2), 'r--', 'LineWidth', 2, 'DisplayName', '平均 Tmax');
    h_tmin_avg = plot(ax3, time_points, mean(Tmin_matrix_p, 2), 'b--', 'LineWidth', 2, 'DisplayName', '平均 Tmin');

    % 3. 图例管理
    legend_handles = [h_tmax_avg, h_tmin_avg];
    if num_AC_participating > 0
        set(h_temp_individual(1), 'DisplayName', '单体空调温度 (反推)');
        legend_handles = [h_temp_individual(1), h_tmax_avg, h_tmin_avg];
    end
    if num_AC_participating > 1
         set(h_temp_individual(2:end), 'HandleVisibility', 'off');
    end
    
    hold(ax3, 'off');
    
    % 4. 格式化
    xlabel(ax3, '时间 (小时)', 'FontSize', 12);
    ylabel(ax3, '室内温度 (°C)', 'FontSize', 12);
    title(ax3, '图3：单体空调室内温度变化 (基于SOC反推)', 'FontSize', 14);
    
    % 5. 创建图例
    legend(ax3, legend_handles, 'Location', 'best', 'FontSize', 11);
    
    set(ax3, 'FontSize', 11);
    xticks(ax3, [0, 6, 12, 18, 24]);
    xticklabels(ax3, {'00:00', '06:00', '12:00', '18:00', '24:00'});
    xlim(ax3, [0, 24]);
    grid(ax3, 'on');
    % --- [新增] 结束 ---


    fprintf('所有仿真和绘图任务完成。\n');
    toc;
end