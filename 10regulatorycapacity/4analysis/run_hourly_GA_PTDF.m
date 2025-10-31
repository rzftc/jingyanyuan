% run_hourly_GA_PTDF.m (已修正)
function run_hourly_GA_PTDF()
    %% 内存优化的分层调度脚本 (版本 8: 小时级GA + 分钟级MILP(含PTDF)调度)
    % 架构:
    % 1. 上层(GA): 优化每小时参与数(n_ac, n_ev)，满足小时峰值需求+优化小时SDCI/Rho。
    % 2. 下层(MILP): 在GA给定的(n_ac, n_ev)数量内，
    %    在每个时间步t，求解成本最低的(u_i, u_j)组合，
    %    同时满足 P_req(t) 约束 和 线路潮流(PTDF)约束。
    clc; close all; clear;

    % --- 1. 加载数据 (*** 已按您的要求修改 ***) ---
    chunk_file_to_load = 'chunk_results/results_chunk_1.mat'; % <--- 加载指定的分块文件
    
    fprintf('正在加载分块数据: %s\n', chunk_file_to_load);
    try
        data_struct = load(chunk_file_to_load);
        if isfield(data_struct, 'results')
            sim_results = data_struct.results;
        else
            error('在 %s 文件中未找到 "results" 结构体。', chunk_file_to_load);
        end
    catch ME
        error('加载数据文件 %s 失败。请确保文件路径正确。\n错误信息: %s', chunk_file_to_load, ME.message);
    end

    % --- 从 results 结构体中提取数据 (假设是上调场景) ---
    if isfield(sim_results, 'AC_Up_Individual')
        P_ac_potential_full = sim_results.AC_Up_Individual;
    else
        P_ac_potential_full = zeros(0, 481); % 假设 T=481
    end
    
    if isfield(sim_results, 'EV_Up_Individual')
        P_ev_potential_full = sim_results.EV_Up_Individual;
    else
        P_ev_potential_full = zeros(0, 481);
    end
    
    % (假设下调也需要)
    if isfield(sim_results, 'AC_Down_Individual')
        P_ac_potential_down_full = abs(sim_results.AC_Down_Individual);
    else
        P_ac_potential_down_full = zeros(0, 481);
    end
    
    if isfield(sim_results, 'EV_Down_Individual')
        P_ev_potential_down_full = abs(sim_results.EV_Down_Individual);
    else
        P_ev_potential_down_full = zeros(0, 481);
    end

    [num_ac_total, T] = size(P_ac_potential_full);
    [num_ev_total, ~] = size(P_ev_potential_full);

    % --- 将数据加载到硬盘 matfile 以模拟原始脚本的内存优化 ---
    individual_file = 'temp_individual_data_for_opt.mat';
    if exist(individual_file, 'file'), delete(individual_file); end
    m_individual = matfile(individual_file, 'Writable', true);
    m_individual.AC_Up_Individual = P_ac_potential_full;
    m_individual.EV_Up_Individual = P_ev_potential_full;
    m_individual.AC_Down_Individual = -P_ac_potential_down_full; % 存为负值
    m_individual.EV_Down_Individual = -P_ev_potential_down_full;
    clear P_ac_potential_full P_ev_potential_full P_ac_potential_down_full P_ev_potential_down_full data_struct sim_results;

    % 计算原始聚合数据
    AC_Up_raw = sum(m_individual.AC_Up_Individual, 1);
    EV_Up_raw = sum(m_individual.EV_Up_Individual, 1);
    AC_Down_raw = sum(abs(m_individual.AC_Down_Individual), 1);
    EV_Down_raw = sum(abs(m_individual.EV_Down_Individual), 1);
    
    fprintf('数据维度加载完成。\nT=%d, N_AC=%d, N_EV=%d\n', T, num_ac_total, num_ev_total);
    
    
    % --- 2. 定义仿真参数 ---
    dt = 5/60; % 假设 5 分钟
    steps_per_hour = round(1/dt);
    num_hours = floor(T / steps_per_hour);
    fprintf('仿真参数: dt=%.3f小时, 每小时步数=%d, 总小时数=%d\n', dt, steps_per_hour, num_hours);

    P_grid_up_demand = (sum(AC_Up_raw) + sum(EV_Up_raw)) / T * (0.2 + 0.3*rand(T,1));
    P_grid_down_demand = (sum(AC_Down_raw) + sum(EV_Down_raw)) / T * (0.15 + 0.25*rand(T,1));
    c_ac_up = ones(num_ac_total, 1) * 0.05; c_ev_up = ones(num_ev_total, 1) * 0.04;
    c_ac_down = ones(num_ac_total, 1) * 0.03; c_ev_down = ones(num_ev_total, 1) * 0.02;
    eps_val = 1e-6;

    % --- 2.5 新增：定义网络拓扑参数 (示例数据) ---
    fprintf('正在定义网络拓扑参数 (示例)...\n');
    N_bus = 10;  % 假设有10个节点
    N_line = 8; % 假设有8条线路
    
    Location_AC = randi([1, N_bus], num_ac_total, 1);
    Location_EV = randi([1, N_bus], num_ev_total, 1);
    PTDF_matrix = rand(N_line, N_bus) * 0.2 - 0.1; % 示例 PTDF 矩阵 (N_line x N_bus)
    P_Line_Base = rand(N_line, T) * 50 - 25; % 示例 线路基础潮流 (N_line x T)
    P_Line_Max = rand(N_line, 1) * 50 + 50; % 示例 线路容量 (N_line x 1), 50-100 kW

    %% 3. 分层优化主循环
    U_ac_up_final = zeros(num_ac_total, T); U_ev_up_final = zeros(num_ev_total, T);
    U_ac_down_final = zeros(num_ac_total, T); U_ev_down_final = zeros(num_ev_total, T);
    cost_up_final = zeros(1, T); cost_down_final = zeros(1, T);
    P_ac_dispatched_opt_up_t = zeros(1,T); P_ev_dispatched_opt_up_t = zeros(1,T);
    P_ac_dispatched_opt_down_t = zeros(1,T); P_ev_dispatched_opt_down_t = zeros(1,T);

    fprintf('\n开始分层优化 (共 %d 小时)...\n', num_hours);
    tic_loop = tic;
    pool = gcp('nocreate');
    if isempty(pool)
        parpool(); % 启动并行池
    end
  
    % --- 主循环改为串行 ---
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

        % --- 3.2 上层GA (保持不变) ---
        fprintf('  上层GA: 优化设备参与数量 (满足峰值需求 + 优化小时级指标)...\n');
        
        % ***【修正 1/2】***
        % 传递一个包含名称-值对的 *单元格数组*
        ga_opts_up = {'Display', 'off', 'UseParallel', true};
        [n_ac_up_hourly, n_ev_up_hourly] = optimize_hourly_participation_GA_constrained_ptdf( ...
            num_ac_total, num_ev_total, P_ac_up_hourly, P_ev_up_hourly, P_req_up_hourly, current_steps_in_hour, ga_opts_up);
        fprintf('  GA结果(上调): n_ac=%d, n_ev=%d\n', n_ac_up_hourly, n_ev_up_hourly);
        
        % ***【修正 2/2】***
        ga_opts_down = {'Display', 'off', 'UseParallel', true};
        [n_ac_down_hourly, n_ev_down_hourly] = optimize_hourly_participation_GA_constrained_ptdf( ...
            num_ac_total, num_ev_total, P_ac_down_hourly, P_ev_down_hourly, P_req_down_hourly, current_steps_in_hour, ga_opts_down);
        fprintf('  GA结果(下调): n_ac=%d, n_ev=%d\n', n_ac_down_hourly, n_ev_down_hourly);

        % --- 3.3 下层MILP调度 ---
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

        % --- 使用 parfor 并行处理时间步 ---
        parfor t_local = 1:current_steps_in_hour
            t_global = current_steps_indices(t_local); % 全局时间索引

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
            
            if flag_up > 0
                U_ac_up_h(:, t_local) = u_ac_up;
                U_ev_up_h(:, t_local) = u_ev_up;
                cost_up_h(t_local) = cost_up;
                p_ac_up_h(t_local) = sum(u_ac_up .* p_ac_t_up);
                p_ev_up_h(t_local) = sum(u_ev_up .* p_ev_t_up);
            else
                cost_up_h(t_local) = NaN;
            end

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

            if flag_down > 0
                U_ac_down_h(:, t_local) = u_ac_down;
                U_ev_down_h(:, t_local) = u_ev_down;
                cost_down_h(t_local) = cost_down;
                p_ac_down_h(t_local) = sum(u_ac_down .* p_ac_t_down);
                p_ev_down_h(t_local) = sum(u_ev_down .* p_ev_t_down);
            else
                cost_down_h(t_local) = NaN;
            end
        end % 结束 parfor t_local

        % --- 将小时结果存入最终结果矩阵 ---
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

    % 清理临时 mat 文件
    delete(individual_file);

    %% 4. 计算优化前后的SDCI和Rho指标
    fprintf('计算指标 (基于全量数据)...\n');
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

    % --- "优化后"场景的SDCI和Rho计算 ---
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

    %% 5. 结果分析与可视化
    fprintf('生成可视化图表...\n');
    time_axis = (1:T)';
    Total_Up_Optimal_Agg = P_ac_dispatched_opt_up_t' + P_ev_dispatched_opt_up_t';
    Total_Down_Optimal_Agg = P_ac_dispatched_opt_down_t' + P_ev_dispatched_opt_down_t';

    figure('Position', [100, 100, 1200, 800]);
    sgtitle(sprintf('VPP 逐时优化结果 (分层GA+MILP N_{total}=%d)', num_ac_total + num_ev_total), 'FontSize', 16);
    subplot(2,2,1); plot(time_axis, Total_Up_Aggregated_Raw, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Raw Aggregated'); hold on; plot(time_axis, Total_Up_Optimal_Agg, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Optimized (Layered+PTDF)'); plot(time_axis, P_grid_up_demand, 'k:', 'LineWidth', 1, 'DisplayName', 'Up-Regulation Demand'); hold off; title('Up-Regulation Power Comparison'); xlabel('Time Step'); ylabel('Power (kW)'); legend('show', 'Location', 'best'); grid on;
    subplot(2,2,2); plot(time_axis, Total_Down_Aggregated_Raw, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Raw Aggregated'); hold on; plot(time_axis, Total_Down_Optimal_Agg, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Optimized (Layered+PTDF)'); plot(time_axis, P_grid_down_demand, 'k:', 'LineWidth', 1, 'DisplayName', 'Down-Regulation Demand'); hold off; title('Down-Regulation Power Comparison'); xlabel('Time Step'); ylabel('Power (kW)'); legend('show', 'Location', 'best'); grid on;
    subplot(2,2,3); indicator_names_up = {'SDCI⁺', 'Spearman ρ⁺'}; values_raw_up = [SDCI_up_raw; rho_up_raw]; values_opt_up = [SDCI_up_opt; rho_up_opt]; bar_data_up = [values_raw_up, values_opt_up]; b_up = bar(bar_data_up); set(gca, 'XTickLabel', indicator_names_up); ylabel('Indicator Value'); ylim([-1.1, 1.1]); legend([b_up(1) b_up(2)], {'Raw Aggregated', 'Optimized (Layered+PTDF)'}, 'Location', 'northoutside', 'Orientation','horizontal'); title('Up-Regulation Complementarity & Correlation'); grid on; add_bar_labels(b_up, bar_data_up);
    subplot(2,2,4); indicator_names_down = {'SDCI⁻', 'Spearman ρ⁻'}; values_raw_down = [SDCI_down_raw; rho_down_raw]; values_opt_down = [SDCI_down_opt; rho_down_opt]; bar_data_down = [values_raw_down, values_opt_down]; b_down = bar(bar_data_down); set(gca, 'XTickLabel', indicator_names_down); ylabel('Indicator Value'); ylim([-1.1, 1.1]); legend([b_down(1) b_down(2)],{'Raw Aggregated', 'Optimized (Layered+PTDF)'}, 'Location', 'northoutside', 'Orientation','horizontal'); title('Down-Regulation Complementarity & Correlation'); grid on; add_bar_labels(b_down, bar_data_down);

    %% 6. 命令行输出汇总
    disp(' ');
    disp('=== 互补性与相关性指标对比 (分层GA + MILP(含PTDF) 优化) ===');
    disp('【上调】'); fprintf('优化前 (Raw Aggregated): SDCI⁺ = %.4f, ρ⁺ = %.4f\n', SDCI_up_raw, rho_up_raw); fprintf('优化后 (Layered Opt): SDCI⁺ = %.4f, ρ⁺ = %.4f, 总成本 = %.2f\n', SDCI_up_opt, rho_up_opt, nansum(cost_up_final));
    disp(' ');
    disp('【下调】'); fprintf('优化前 (Raw Aggregated): SDCI⁻ = %.4f, ρ⁻ = %.4f\n', SDCI_down_raw, rho_down_raw); fprintf('优化后 (Layered Opt): SDCI⁻ = %.4f, ρ⁻ = %.4f, 总成本 = %.2f\n', SDCI_down_opt, rho_down_opt, nansum(cost_down_final));
    disp(' ');

end % 结束主函数

%% -----------------------------------------------------------------------
%                       局部辅助函数
% ------------------------------------------------------------------------

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