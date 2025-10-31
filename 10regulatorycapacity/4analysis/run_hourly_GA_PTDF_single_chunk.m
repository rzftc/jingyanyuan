% run_hourly_GA_PTDF_single_chunk.m
function run_hourly_GA_PTDF_single_chunk()
    %% 内存优化的分层调度脚本 (版本 8.3: 修复 optimoptions 错误和 matfile 警告)
    % 架构:
    % 1. [修改] 仅加载一个指定的分块 (chunk) .mat 文件。
    % 2. 上层(GA): 优化该分块的每小时参与数(n_ac, n_ev)，满足小时峰值需求+优化小时SDCI/Rho。
    % 3. 下层(MILP): 在GA给定的(n_ac, n_ev)数量内，
    %    在每个时间步t，求解该分块中成本最低的(u_i, u_j)组合，
    %    同时满足 P_req(t) 约束 和 线路潮流(PTDF)约束。
    % *** (新增) 4. MILP求解失败时，回退到 贪心算法 (Greedy) 确保功率满足 ***
    clc; close all; clear;

    % --- 1. [MODIFIED] 加载并准备 *单个* 分块数据 ---
    
    % *** [用户配置] ***
    % 请在此处指定要加载的单个分块文件:
    chunk_to_load = 'chunk_results/results_chunk_1.mat'; 
    % *****************
    
    aggregate_file_output = 'aggregate_results_slow_synced.mat'; % (基于该分块的聚合文件)
    individual_file_output = 'individual_results_slow_synced.mat'; % (基于该分块的个体文件)
    
    fprintf('正在从单个分块文件加载数据: %s\n', chunk_to_load);
    if ~exist(chunk_to_load, 'file')
        error('指定的分块文件未找到: %s', chunk_to_load);
    end

    % 加载分块数据
    try
        chunk_data = load(chunk_to_load);
        if ~isfield(chunk_data, 'results')
            error('文件 %s 中缺少 "results" 结构体。', chunk_to_load);
        end
        results_chunk = chunk_data.results;
    catch ME
        error('加载文件 %s 时出错: %s', chunk_to_load, ME.message);
    end

    % 1.1 提取聚合数据
    fprintf('  提取聚合数据...\n');
    results_agg = struct();
    T = 0; % 初始化时间步
    
    % 确定T (从AC或EV个体数据中)
    if isfield(results_chunk, 'AC_Up_Individual') && ~isempty(results_chunk.AC_Up_Individual)
        T = size(results_chunk.AC_Up_Individual, 2);
    elseif isfield(results_chunk, 'EV_Up_Individual') && ~isempty(results_chunk.EV_Up_Individual)
        T = size(results_chunk.EV_Up_Individual, 2);
    elseif isfield(results_chunk, 'AC_Up') % 后备
        T = length(results_chunk.AC_Up);
    elseif isfield(results_chunk, 'EV_Up') % 后备
        T = length(results_chunk.EV_Up);
    end
    
    if T == 0; error('未能从分块文件中确定时间步长 T。'); end
    
    % 提取AC聚合数据
    if isfield(results_chunk, 'AC_Up')
        results_agg.AC_Up = results_chunk.AC_Up;
        results_agg.AC_Down = results_chunk.AC_Down;
    else
        results_agg.AC_Up = zeros(T, 1);
        results_agg.AC_Down = zeros(T, 1);
    end
    
    % 提取EV聚合数据
    if isfield(results_chunk, 'EV_Up')
        results_agg.EV_Up = results_chunk.EV_Up;
        results_agg.EV_Down = results_chunk.EV_Down;
    else
        results_agg.EV_Up = zeros(T, 1);
        results_agg.EV_Down = zeros(T, 1);
    end
    
    % 保存这个分块的聚合数据 (使用'results'变量名以兼容)
    fprintf('  保存该分块的聚合数据到: %s\n', aggregate_file_output);
    results = results_agg; % 兼容 generate_demand 函数
    save(aggregate_file_output, 'results');

    % 1.2 提取并保存个体数据到 matfile
    fprintf('  正在创建并写入硬盘映射文件 (matfile): %s\n', individual_file_output);
    try
        if exist(individual_file_output, 'file'), delete(individual_file_output); end
        m = matfile(individual_file_output, 'Writable', true);
        
        % 写入 AC 数据
        if isfield(results_chunk, 'AC_Up_Individual') && ~isempty(results_chunk.AC_Up_Individual)
            m.AC_Up_Individual = results_chunk.AC_Up_Individual;
            m.AC_Down_Individual = results_chunk.AC_Down_Individual;
            num_ac_total = size(results_chunk.AC_Up_Individual, 1);
        else
            m.AC_Up_Individual = zeros(0, T); 
            m.AC_Down_Individual = zeros(0, T);
            num_ac_total = 0;
        end
        
        % 写入 EV 数据
        if isfield(results_chunk, 'EV_Up_Individual') && ~isempty(results_chunk.EV_Up_Individual)
            m.EV_Up_Individual = results_chunk.EV_Up_Individual;
            m.EV_Down_Individual = results_chunk.EV_Down_Individual;
            num_ev_total = size(results_chunk.EV_Up_Individual, 1);
        else
            m.EV_Up_Individual = zeros(0, T);
            m.EV_Down_Individual = zeros(0, T);
            num_ev_total = 0;
        end
        
        % *** [FIX 2] 确定并保存 dt ***
        if isfield(results_chunk, 'dt') && ~isempty(results_chunk.dt) % 假设 ac_ev_simulation_block.m 保存了 dt
            dt_from_chunk = results_chunk.dt;
        else
            dt_from_chunk = 5/60; % 默认 5 分钟
            warning('警告: 分块结果文件中未找到 dt, 使用默认值 %.3f 小时。', dt_from_chunk);
        end
        m.dt = dt_from_chunk; % 保存 dt 到 matfile
        
    catch E
        error('创建或写入 %s 失败: %s', individual_file_output, E.message);
    end
    
    % 1.3 *** [FIX 1] *** 创建 matfile 对象供下游使用
    % 必须在 'm' (写入句柄) 被清除或关闭 (隐式) 之后创建
    clear m; % 确保写入缓冲区已刷新并关闭文件
    fprintf('  完成写入, 创建 matfile 读取句柄...\n');
    m_individual = matfile(individual_file_output); % <-- 移动到此处
    
    % 1.4 准备原脚本的剩余变量
    AC_Up_raw = results_agg.AC_Up(:)'; 
    AC_Down_raw = abs(results_agg.AC_Down(:)');
    EV_Up_raw = results_agg.EV_Up(:)'; 
    EV_Down_raw = abs(results_agg.EV_Down(:)');
    
    fprintf('数据加载和预处理完成 (仅使用 %s 的数据)。\n', chunk_to_load);
    fprintf('总细分时间步数 T = %d\n', T);
    fprintf('该分块空调数 num_ac_total = %d\n', num_ac_total);
    fprintf('该分块电动汽车数 num_ev_total = %d\n', num_ev_total);
    
    clear chunk_data results_chunk; % 释放内存
    % --- [数据加载部分结束] ---


    % --- 2. 定义仿真参数 ---
    % *** [FIX 2 Cont.] ***
    if isprop(m_individual,'dt') && ~isempty(m_individual.dt) % Use isprop for matfile
        dt = m_individual.dt; % 从新创建的 matfile 中读取 dt
    else
        dt = 5/60; % 默认 5 分钟
        warning('matfile 中未找到 dt, 使用默认值 %.3f 小时。', dt);
    end
    steps_per_hour = round(1/dt);
    num_hours = floor(T / steps_per_hour);
    fprintf('仿真参数: dt=%.3f小时, 每小时步数=%d, 总小时数=%d\n', dt, steps_per_hour, num_hours);

    P_grid_up_demand = generate_demand(results_agg, 'P_grid_up_regulation_demand', T, (sum(AC_Up_raw) + sum(EV_Up_raw)) / T, [0.2, 0.5], 100);
    P_grid_down_demand = generate_demand(results_agg, 'P_grid_down_regulation_demand', T, (sum(AC_Down_raw) + sum(EV_Down_raw)) / T, [0.15, 0.4], 80);
    c_ac_up = ones(num_ac_total, 1) * 0.05; c_ev_up = ones(num_ev_total, 1) * 0.04;
    c_ac_down = ones(num_ac_total, 1) * 0.03; c_ev_down = ones(num_ev_total, 1) * 0.02;
    eps_val = 1e-6;

    % --- 2.5 定义网络拓扑参数 (示例数据) (保持不变) ---
    fprintf('正在定义网络拓扑参数 (示例)...\n');
    N_bus = 10;
    N_line = 8;
    Location_AC = randi([1, N_bus], num_ac_total, 1);
    Location_EV = randi([1, N_bus], num_ev_total, 1);
    PTDF_matrix = rand(N_line, N_bus) * 0.2 - 0.1; 
    P_Line_Base = rand(N_line, T) * 50 - 25; 
    P_Line_Max = rand(N_line, 1) * 50 + 50;
    % --- 结束 ---

    %% 3. 分层优化主循环 (保持不变)
    % 初始化最终结果存储
    U_ac_up_final = zeros(num_ac_total, T); U_ev_up_final = zeros(num_ev_total, T);
    U_ac_down_final = zeros(num_ac_total, T); U_ev_down_final = zeros(num_ev_total, T);
    cost_up_final = zeros(1, T); cost_down_final = zeros(1, T);
    P_ac_dispatched_opt_up_t = zeros(1,T); P_ev_dispatched_opt_up_t = zeros(1,T);
    P_ac_dispatched_opt_down_t = zeros(1,T); P_ev_dispatched_opt_down_t = zeros(1,T);

    fprintf('\n开始分层优化 (共 %d 小时)...\n', num_hours);
    tic_loop = tic;
    pool = gcp('nocreate');
  
    for h = 1:num_hours
        fprintf('--- 正在优化第 %d 小时 ---\n', h);

        % --- 3.1 定义当前小时的时间范围和数据 ---
        start_step = (h-1) * steps_per_hour + 1;
        end_step = h * steps_per_hour;
        if end_step > T; end_step = T; end
        current_steps_indices = start_step:end_step;
        current_steps_in_hour = length(current_steps_indices);

        fprintf('  加载第 %d 小时数据...\n', h);
        P_ac_up_hourly = double(m_individual.AC_Up_Individual(:, current_steps_indices));
        P_ev_up_hourly = double(m_individual.EV_Up_Individual(:, current_steps_indices));
        P_ac_down_hourly = double(abs(m_individual.AC_Down_Individual(:, current_steps_indices)));
        P_ev_down_hourly = double(abs(m_individual.EV_Down_Individual(:, current_steps_indices)));
        P_req_up_hourly = P_grid_up_demand(current_steps_indices);
        P_req_down_hourly = P_grid_down_demand(current_steps_indices);
        fprintf('  数据加载完毕.\n');

        % --- 3.2 上层GA ---
        fprintf('  上层GA: 优化设备参与数量 (满足峰值需求 + 优化小时级指标)...\n');
        
        % *** [FIX 3] ***
        % 传递 Name-Value 对 (Cell 数组)，而不是 options 对象
        ga_opts_up = {'Display', 'off', 'UseParallel', true};
        [n_ac_up_hourly, n_ev_up_hourly] = optimize_hourly_participation_GA_constrained_ptdf( ...
            num_ac_total, num_ev_total, P_ac_up_hourly, P_ev_up_hourly, P_req_up_hourly, current_steps_in_hour, ga_opts_up);
        fprintf('  GA结果(上调): n_ac=%d, n_ev=%d\n', n_ac_up_hourly, n_ev_up_hourly);

        % *** [FIX 3] ***
        ga_opts_down = {'Display', 'off', 'UseParallel', true};
        [n_ac_down_hourly, n_ev_down_hourly] = optimize_hourly_participation_GA_constrained_ptdf( ... % 使用 ptdf 版本的函数名
            num_ac_total, num_ev_total, P_ac_down_hourly, P_ev_down_hourly, P_req_down_hourly, current_steps_in_hour, ga_opts_down);
        fprintf('  GA结果(下调): n_ac=%d, n_ev=%d\n', n_ac_down_hourly, n_ev_down_hourly);

        % --- 3.3 下层MILP调度 (保持不变) ---
        fprintf('  下层MILP: 在 %d 个时间步内执行约束调度 (含PTDF)...\n', current_steps_in_hour);
        
        U_ac_up_h = zeros(num_ac_total, current_steps_in_hour);
        U_ev_up_h = zeros(num_ev_total, current_steps_in_hour);
        cost_up_h = zeros(1, current_steps_in_hour);
        p_ac_up_h = zeros(1, current_steps_in_hour);
        p_ev_up_h = zeros(1, current_steps_in_hour);

        U_ac_down_h = zeros(num_ac_total, current_steps_in_hour);
        U_ev_down_h = zeros(num_ev_total, current_steps_in_hour);
        cost_down_h = zeros(1, current_steps_in_hour);
        p_ac_down_h = zeros(1, current_steps_in_hour);
        p_ev_down_h = zeros(1, current_steps_in_hour);

        parfor t_local = 1:current_steps_in_hour
            t_global = current_steps_indices(t_local); 

            % --- 上调调度 ---
            p_ac_t_up = P_ac_up_hourly(:, t_local);
            p_ev_t_up = P_ev_up_hourly(:, t_local);
            
            [u_ac_up, u_ev_up, cost_up, flag_up] = solve_hourly_dispatch_ptdf( ...
                num_ac_total, num_ev_total, ...
                p_ac_t_up, p_ev_t_up, ...
                c_ac_up, c_ev_up, P_req_up_hourly(t_local), ...
                n_ac_up_hourly, n_ev_up_hourly, ... 
                Location_AC, Location_EV, PTDF_matrix, ... 
                P_Line_Base(:, t_global), P_Line_Max, N_bus, N_line); 
            
            %% --- MODIFICATION START (UP) ---
            if flag_up > 0 % MILP(PTDF) 成功
                U_ac_up_h(:, t_local) = u_ac_up;
                U_ev_up_h(:, t_local) = u_ev_up;
                cost_up_h(t_local) = cost_up;
                p_ac_up_h(t_local) = sum(u_ac_up .* p_ac_t_up);
                p_ev_up_h(t_local) = sum(u_ev_up .* p_ev_t_up);
            else % MILP(PTDF) 失败，回退到贪心算法 (不考虑PTDF)
                 % 在 parfor 中使用 fprintf 可能导致输出交错，但对于调试是可接受的
                 fprintf('  警告: 时段 %d (全局) MILP(PTDF) 上调求解失败 (Flag=%d)。回退到贪心算法。\n', t_global, flag_up);
                 
                 % 调用贪心算法 (需要 solve_hourly_dispatch_greedy_with_count.m 在路径中)
                 [u_ac_greedy, u_ev_greedy, cost_greedy] = solve_hourly_dispatch_greedy_with_count( ...
                    p_ac_t_up, p_ev_t_up, c_ac_up, c_ev_up, P_req_up_hourly(t_local), n_ac_up_hourly, n_ev_up_hourly);
            
                U_ac_up_h(:, t_local) = u_ac_greedy;
                U_ev_up_h(:, t_local) = u_ev_greedy;
                cost_up_h(t_local) = cost_greedy; % 使用贪心算法的成本
                p_ac_up_h(t_local) = sum(u_ac_greedy .* p_ac_t_up);
                p_ev_up_h(t_local) = sum(u_ev_greedy .* p_ev_t_up);
            end
            %% --- MODIFICATION END (UP) ---


            % --- 下调调度 ---
            p_ac_t_down = P_ac_down_hourly(:, t_local);
            p_ev_t_down = P_ev_down_hourly(:, t_local);

            [u_ac_down, u_ev_down, cost_down, flag_down] = solve_hourly_dispatch_ptdf( ...
                num_ac_total, num_ev_total, ...
                p_ac_t_down, p_ev_t_down, ...
                c_ac_down, c_ev_down, P_req_down_hourly(t_local), ...
                n_ac_down_hourly, n_ev_down_hourly, ... 
                Location_AC, Location_EV, PTDF_matrix, ... 
                P_Line_Base(:, t_global), P_Line_Max, N_bus, N_line); 

            %% --- MODIFICATION START (DOWN) ---
            if flag_down > 0 % MILP(PTDF) 成功
                U_ac_down_h(:, t_local) = u_ac_down;
                U_ev_down_h(:, t_local) = u_ev_down;
                cost_down_h(t_local) = cost_down;
                p_ac_down_h(t_local) = sum(u_ac_down .* p_ac_t_down);
                p_ev_down_h(t_local) = sum(u_ev_down .* p_ev_t_down);
            else % MILP(PTDF) 失败，回退到贪心算法 (不考虑PTDF)
                fprintf('  警告: 时段 %d (全局) MILP(PTDF) 下调求解失败 (Flag=%d)。回退到贪心算法。\n', t_global, flag_down);
                
                [u_ac_greedy_d, u_ev_greedy_d, cost_greedy_d] = solve_hourly_dispatch_greedy_with_count( ...
                    p_ac_t_down, p_ev_t_down, c_ac_down, c_ev_down, P_req_down_hourly(t_local), n_ac_down_hourly, n_ev_down_hourly);
                    
                U_ac_down_h(:, t_local) = u_ac_greedy_d;
                U_ev_down_h(:, t_local) = u_ev_greedy_d;
                cost_down_h(t_local) = cost_greedy_d;
                p_ac_down_h(t_local) = sum(u_ac_greedy_d .* p_ac_t_down);
                p_ev_down_h(t_local) = sum(u_ev_greedy_d .* p_ev_t_down);
            end
            %% --- MODIFICATION END (DOWN) ---
            
        end % 结束 parfor t_local

        U_ac_up_final(:, current_steps_indices) = U_ac_up_h;
        U_ev_up_final(:, current_steps_indices) = U_ev_up_h;
        cost_up_final(current_steps_indices) = cost_up_h;
        P_ac_dispatched_opt_up_t(current_steps_indices) = p_ac_up_h;
        P_ev_dispatched_opt_up_t(current_steps_indices) = p_ev_up_h;

        U_ac_down_final(:, current_steps_indices) = U_ac_down_h;
        U_ev_down_final(:, current_steps_indices) = U_ev_down_h;
        cost_down_final(current_steps_indices) = cost_down_h;
        P_ac_dispatched_opt_down_t(current_steps_indices) = p_ac_down_h;
        P_ev_dispatched_opt_down_t(current_steps_indices) = p_ev_down_h;
        
        clear P_ac_up_hourly P_ev_up_hourly P_ac_down_hourly P_ev_down_hourly ...
              U_ac_up_h U_ev_up_h cost_up_h p_ac_up_h p_ev_up_h ...
              U_ac_down_h U_ev_down_h cost_down_h p_ac_down_h p_ev_down_h;
    end % 结束 for h (小时循环)

    toc(tic_loop);
    fprintf('分层优化(含PTDF)完成。\n');

    %% 4. 计算优化前后的SDCI和Rho指标 (保持不变)
    fprintf('计算指标 (基于分块数据)...\n');
    n_ac_raw_t = ones(T,1) * num_ac_total; n_ev_raw_t = ones(T,1) * num_ev_total;
    if num_ac_total == 0; n_ac_raw_t = zeros(T,1); end; if num_ev_total == 0; n_ev_raw_t = zeros(T,1); end

    AC_Up_Aggregated_Raw = AC_Up_raw'; EV_Up_Aggregated_Raw = EV_Up_raw';
    AC_Down_Aggregated_Raw = AC_Down_raw'; EV_Down_Aggregated_Raw = EV_Down_raw';
    Total_Up_Aggregated_Raw = AC_Up_Aggregated_Raw + EV_Up_Aggregated_Raw;
    Total_Down_Aggregated_Raw = AC_Down_Aggregated_Raw + EV_Down_Aggregated_Raw;

    avg_P_ac_raw_up_t = AC_Up_Aggregated_Raw ./ (n_ac_raw_t + eps_val);
    avg_P_ev_raw_up_t = EV_Up_Aggregated_Raw ./ (n_ev_raw_t + eps_val);
    avg_P_ac_raw_down_t = AC_Down_Aggregated_Raw ./ (n_ac_raw_t + eps_val);
    avg_P_ev_raw_down_t = EV_Down_Aggregated_Raw ./ (n_ev_raw_t + eps_val);

    SDCI_up_raw = ensureScalar(calculateSDCI(n_ac_raw_t, n_ev_raw_t, avg_P_ac_raw_up_t, avg_P_ev_raw_up_t));
    rho_up_raw  = ensureScalar(calculateSpearmanRho(n_ac_raw_t, avg_P_ac_raw_up_t, n_ev_raw_t, avg_P_ev_raw_up_t));
    SDCI_down_raw = ensureScalar(calculateSDCI(n_ac_raw_t, n_ev_raw_t, avg_P_ac_raw_down_t, avg_P_ev_raw_down_t));
    rho_down_raw  = ensureScalar(calculateSpearmanRho(n_ac_raw_t, avg_P_ac_raw_down_t, n_ev_raw_t, avg_P_ev_raw_down_t));

    n_ac_opt_up_count_t = sum(U_ac_up_final, 1)'; n_ev_opt_up_count_t = sum(U_ev_up_final, 1)';
    n_ac_opt_down_count_t = sum(U_ac_down_final, 1)'; n_ev_opt_down_count_t = sum(U_ev_down_final, 1)';
    avg_P_ac_opt_up = P_ac_dispatched_opt_up_t' ./ (n_ac_opt_up_count_t + eps_val);
    avg_P_ev_opt_up = P_ev_dispatched_opt_up_t' ./ (n_ev_opt_up_count_t + eps_val);
    avg_P_ac_opt_down = P_ac_dispatched_opt_down_t' ./ (n_ac_opt_down_count_t + eps_val);
    avg_P_ev_opt_down = P_ev_dispatched_opt_down_t' ./ (n_ev_opt_down_count_t + eps_val);
    SDCI_up_opt = ensureScalar(calculateSDCI(n_ac_opt_up_count_t, n_ev_opt_up_count_t, avg_P_ac_opt_up, avg_P_ev_opt_up));
    rho_up_opt  = ensureScalar(calculateSpearmanRho(n_ac_opt_up_count_t, avg_P_ac_opt_up, n_ev_opt_up_count_t, avg_P_ev_opt_up));
    SDCI_down_opt = ensureScalar(calculateSDCI(n_ac_opt_down_count_t, n_ev_opt_down_count_t, avg_P_ac_opt_down, avg_P_ev_opt_down));
    rho_down_opt  = ensureScalar(calculateSpearmanRho(n_ac_opt_down_count_t, avg_P_ac_opt_down, n_ev_opt_down_count_t, avg_P_ev_opt_down));

    %% 5. 结果分析与可视化 (保持不变)
    fprintf('生成可视化图表...\n');
    time_axis = (1:T)';
    Total_Up_Optimal_Agg = P_ac_dispatched_opt_up_t' + P_ev_dispatched_opt_up_t';
    Total_Down_Optimal_Agg = P_ac_dispatched_opt_down_t' + P_ev_dispatched_opt_down_t';

    figure('Position', [100, 100, 1200, 800]);
    sgtitle(sprintf('VPP 逐时优化结果 (分层GA+MILP, N_{chunk}=%d)', num_ac_total + num_ev_total), 'FontSize', 16);
    subplot(2,2,1); plot(time_axis, Total_Up_Aggregated_Raw, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Raw Aggregated (Chunk)'); hold on; plot(time_axis, Total_Up_Optimal_Agg, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Optimized (Layered+PTDF)'); plot(time_axis, P_grid_up_demand, 'k:', 'LineWidth', 1, 'DisplayName', 'Up-Regulation Demand'); hold off; title('Up-Regulation Power Comparison'); xlabel('Time Step'); ylabel('Power (kW)'); legend('show', 'Location', 'best'); grid on;
    subplot(2,2,2); plot(time_axis, Total_Down_Aggregated_Raw, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Raw Aggregated (Chunk)'); hold on; plot(time_axis, Total_Down_Optimal_Agg, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Optimized (Layered+PTDF)'); plot(time_axis, P_grid_down_demand, 'k:', 'LineWidth', 1, 'DisplayName', 'Down-Regulation Demand'); hold off; title('Down-Regulation Power Comparison'); xlabel('Time Step'); ylabel('Power (kW)'); legend('show', 'Location', 'best'); grid on;
    subplot(2,2,3); indicator_names_up = {'SDCI⁺', 'Spearman ρ⁺'}; values_raw_up = [SDCI_up_raw; rho_up_raw]; values_opt_up = [SDCI_up_opt; rho_up_opt]; bar_data_up = [values_raw_up, values_opt_up]; b_up = bar(bar_data_up); set(gca, 'XTickLabel', indicator_names_up); ylabel('Indicator Value'); ylim([-1.1, 1.1]); legend([b_up(1) b_up(2)], {'Raw Aggregated', 'Optimized (Layered+PTDF)'}, 'Location', 'northoutside', 'Orientation','horizontal'); title('Up-Regulation Complementarity & Correlation'); grid on; add_bar_labels(b_up, bar_data_up);
    subplot(2,2,4); indicator_names_down = {'SDCI⁻', 'Spearman ρ⁻'}; values_raw_down = [SDCI_down_raw; rho_down_raw]; values_opt_down = [SDCI_down_opt; rho_down_opt]; bar_data_down = [values_raw_down, values_opt_down]; b_down = bar(bar_data_down); set(gca, 'XTickLabel', indicator_names_down); ylabel('Indicator Value'); ylim([-1.1, 1.1]); legend([b_down(1) b_down(2)],{'Raw Aggregated', 'Optimized (Layered+PTDF)'}, 'Location', 'northoutside', 'Orientation','horizontal'); title('Down-Regulation Complementarity & Correlation'); grid on; add_bar_labels(b_down, bar_data_down);

    %% 6. 命令行输出汇总 (保持不变)
    disp(' ');
    disp('=== 互补性与相关性指标对比 (分层GA + MILP(含PTDF) 优化 - 单个分块) ===');
    disp('【上调】'); fprintf('优化前 (Raw Aggregated): SDCI⁺ = %.4f, ρ⁺ = %.4f\n', SDCI_up_raw, rho_up_raw); fprintf('优化后 (Layered Opt): SDCI⁺ = %.4f, ρ⁺ = %.4f, 总成本 = %.2f\n', SDCI_up_opt, rho_up_opt, nansum(cost_up_final));
    disp(' ');
    disp('【下调】'); fprintf('优化前 (Raw Aggregated): SDCI⁻ = %.4f, ρ⁻ = %.4f\n', SDCI_down_raw, rho_down_raw); fprintf('优化后 (Layered Opt): SDCI⁻ = %.4f, ρ⁻ = %.4f, 总成本 = %.2f\n', SDCI_down_opt, rho_down_opt, nansum(cost_down_final));
    disp(' ');

end % 结束主函数


%% 局部函数: generate_demand 
function P_demand = generate_demand(results_agg, field_name, T, avg_raw_power, factors, default_avg)
    if isfield(results_agg, field_name) && ~isempty(results_agg.(field_name)) && length(results_agg.(field_name)) == T
        P_demand = results_agg.(field_name)(:);
    else
        if isnan(avg_raw_power) || avg_raw_power <= 0; avg_raw_power = default_avg; end
        P_demand = avg_raw_power * (factors(1) + (factors(2)-factors(1)) * rand(T,1));
        fprintf('未在聚合文件中找到 "%s"，已生成示例需求。\n', field_name);
    end
end

%% 局部函数: ensureScalar
function val = ensureScalar(inputVal)
    if isscalar(inputVal); val = inputVal;
    elseif isempty(inputVal); val = NaN;
    else; val = mean(inputVal(:), 'omitnan'); end
end

%% 局部函数: add_bar_labels
function add_bar_labels(bar_handles, bar_data)
    for k_bar = 1:size(bar_data, 1)
        for j_bar = 1:size(bar_data, 2)
             try
                 if ~ishandle(bar_handles(j_bar)); continue; end
                 if k_bar > length(bar_handles(j_bar).XData); continue; end
                 if isnan(bar_data(k_bar, j_bar)); continue; end
                 
                text(bar_handles(j_bar).XData(k_bar) + bar_handles(j_bar).XOffset, ...
                     bar_data(k_bar, j_bar), sprintf('%.3f', bar_data(k_bar, j_bar)), ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                     'FontSize', 8, 'Color', 'k');
             catch
             end
        end
    end
end