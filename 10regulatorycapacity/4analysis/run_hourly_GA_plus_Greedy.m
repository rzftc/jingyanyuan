function run_hourly_GA_plus_Greedy()
    %% 内存优化的分层调度脚本 (版本 7: 强化功率约束 + 小时级GA + 分钟级贪心)
    % 上层GA优化每小时参与数，首要满足小时峰值需求，然后优化小时级SDCI/Rho
    % 下层贪心根据GA数量选择成本最低设备满足实时需求
    clc; close all; clear;

    % --- 1. 加载数据 (同前) ---
    aggregate_file = 'aggregate_results_slow_synced.mat';
    individual_file = 'individual_results_slow_synced.mat';
    fprintf('正在加载聚合数据: %s\n', aggregate_file);
    if exist(aggregate_file, 'file'); agg_data = load(aggregate_file);
    else; error('聚合数据文件 %s 未找到。', aggregate_file); end
    if ~isfield(agg_data, 'results'); error('文件 %s 中未找到 "results" 结构体。', aggregate_file); end
    results_agg = agg_data.results;
    fprintf('创建个体数据文件对象 (matfile): %s\n', individual_file);
    if exist(individual_file, 'file'); m_individual = matfile(individual_file);
    else; error('个体数据文件 %s 未找到。', individual_file); end
    fprintf('正在读取维度信息...\n');
    try % 尝试读取维度
        s = whos(m_individual); ac_up_info = s(strcmp({s.name}, 'AC_Up_Individual')); ev_up_info = s(strcmp({s.name}, 'EV_Up_Individual'));
        num_ac_total = ac_up_info.size(1); num_ev_total = ev_up_info.size(1); T = ac_up_info.size(2);
        if T ~= ev_up_info.size(2) || T ~= length(results_agg.AC_Up); warning('个体数据和聚合数据的时间步长不匹配。'); end
    catch ME; error('从 %s 读取维度信息时出错: %s', individual_file, ME.message); end

    AC_Up_raw = results_agg.AC_Up(:)'; AC_Down_raw = abs(results_agg.AC_Down(:)');
    EV_Up_raw = results_agg.EV_Up(:)'; EV_Down_raw = abs(results_agg.EV_Down(:)');
    fprintf('数据维度加载完成。\nT=%d, N_AC=%d, N_EV=%d\n', T, num_ac_total, num_ev_total);

    % --- 2. 定义仿真参数 ---
    % 检查 dt 是否存在于聚合结果中，否则设定默认值
    if isfield(results_agg,'dt') && ~isempty(results_agg.dt)
        dt = results_agg.dt;
    else
        dt = 5/60; % 默认 5 分钟
        warning('聚合结果文件中未找到 dt, 使用默认值 %.3f 小时。', dt);
    end
    steps_per_hour = round(1/dt);
    num_hours = floor(T / steps_per_hour);
    fprintf('仿真参数: dt=%.3f小时, 每小时步数=%d, 总小时数=%d\n', dt, steps_per_hour, num_hours);

    P_grid_up_demand = generate_demand(results_agg, 'P_grid_up_regulation_demand', T, (sum(AC_Up_raw) + sum(EV_Up_raw)) / T, [0.2, 0.5], 100);
    P_grid_down_demand = generate_demand(results_agg, 'P_grid_down_regulation_demand', T, (sum(AC_Down_raw) + sum(EV_Down_raw)) / T, [0.15, 0.4], 80);
    c_ac_up = ones(num_ac_total, 1) * 0.05; c_ev_up = ones(num_ev_total, 1) * 0.04;
    c_ac_down = ones(num_ac_total, 1) * 0.03; c_ev_down = ones(num_ev_total, 1) * 0.02;
    eps_val = 1e-6;

    %% 3. 分层优化主循环
    % 初始化最终结果存储
    U_ac_up_final = zeros(num_ac_total, T); U_ev_up_final = zeros(num_ev_total, T);
    U_ac_down_final = zeros(num_ac_total, T); U_ev_down_final = zeros(num_ev_total, T);
    cost_up_final = zeros(1, T); cost_down_final = zeros(1, T);

    fprintf('\n开始分层优化 (共 %d 小时)...\n', num_hours);
    tic_loop = tic;

     % 启动并行池 (主循环非并行，GA内部和贪心可以并行)
    pool = gcp('nocreate');
  
    % --- 主循环改为串行，GA内部可以并行 ---
    for h = 1:num_hours
        fprintf('--- 正在优化第 %d 小时 ---\n', h);

        % --- 3.1 定义当前小时的时间范围和数据 ---
        start_step = (h-1) * steps_per_hour + 1;
        end_step = h * steps_per_hour;
        if end_step > T; end_step = T; end % 处理最后不足一小时的情况
        current_steps_indices = start_step:end_step;
        current_steps_in_hour = length(current_steps_indices); % 当前小时的实际步数

        % 加载当前小时所需的全量数据
        fprintf('  加载第 %d 小时数据...\n', h);
        P_ac_up_hourly = double(m_individual.AC_Up_Individual(:, current_steps_indices));
        P_ev_up_hourly = double(m_individual.EV_Up_Individual(:, current_steps_indices));
        P_ac_down_hourly = double(abs(m_individual.AC_Down_Individual(:, current_steps_indices)));
        P_ev_down_hourly = double(abs(m_individual.EV_Down_Individual(:, current_steps_indices)));
        P_req_up_hourly = P_grid_up_demand(current_steps_indices);
        P_req_down_hourly = P_grid_down_demand(current_steps_indices);
        fprintf('  数据加载完毕.\n');

        % --- 3.2 上层GA：优化本小时的AC/EV参与数量 (加入功率约束) ---
        fprintf('  上层GA: 优化设备参与数量 (满足峰值需求 + 优化小时级指标)...\n');

        % 上调
        ga_opts_up = optimoptions('ga', 'Display', 'off', 'UseParallel', true); % GA内部并行
        [n_ac_up_hourly, n_ev_up_hourly] = optimize_hourly_participation_GA_constrained( ...
            num_ac_total, num_ev_total, P_ac_up_hourly, P_ev_up_hourly, P_req_up_hourly, current_steps_in_hour, ga_opts_up);
        fprintf('  GA结果(上调): n_ac=%d, n_ev=%d\n', n_ac_up_hourly, n_ev_up_hourly);

        % 下调
        ga_opts_down = optimoptions('ga', 'Display', 'off', 'UseParallel', true); % GA内部并行
        [n_ac_down_hourly, n_ev_down_hourly] = optimize_hourly_participation_GA_constrained( ...
            num_ac_total, num_ev_total, P_ac_down_hourly, P_ev_down_hourly, P_req_down_hourly, current_steps_in_hour, ga_opts_down);
        fprintf('  GA结果(下调): n_ac=%d, n_ev=%d\n', n_ac_down_hourly, n_ev_down_hourly);


        % --- 3.3 下层贪心：根据GA决定的数量，在每个时间步执行调度 ---
        fprintf('  下层贪心: 在 %d 个时间步内执行调度...\n', current_steps_in_hour);
        
        % --- 为下层循环预分配结果 ---
        U_ac_up_h = zeros(num_ac_total, current_steps_in_hour);
        U_ev_up_h = zeros(num_ev_total, current_steps_in_hour);
        cost_up_h = zeros(1, current_steps_in_hour);
        U_ac_down_h = zeros(num_ac_total, current_steps_in_hour);
        U_ev_down_h = zeros(num_ev_total, current_steps_in_hour);
        cost_down_h = zeros(1, current_steps_in_hour);

        % --- 使用 parfor 并行处理时间步 ---
        parfor t_local = 1:current_steps_in_hour
            t_global = current_steps_indices(t_local); % 全局时间索引

            % 上调调度
            p_ac_t_up = P_ac_up_hourly(:, t_local);
            p_ev_t_up = P_ev_up_hourly(:, t_local);
            [u_ac_up, u_ev_up, cost_up] = solve_hourly_dispatch_greedy_with_count( ...
                p_ac_t_up, p_ev_t_up, c_ac_up, c_ev_up, P_req_up_hourly(t_local), n_ac_up_hourly, n_ev_up_hourly);
            U_ac_up_h(:, t_local) = u_ac_up; % 存储到临时矩阵
            U_ev_up_h(:, t_local) = u_ev_up;
            cost_up_h(t_local) = cost_up;

            % 下调调度
            p_ac_t_down = P_ac_down_hourly(:, t_local);
            p_ev_t_down = P_ev_down_hourly(:, t_local);
            [u_ac_down, u_ev_down, cost_down] = solve_hourly_dispatch_greedy_with_count( ...
                p_ac_t_down, p_ev_t_down, c_ac_down, c_ev_down, P_req_down_hourly(t_local), n_ac_down_hourly, n_ev_down_hourly);
            U_ac_down_h(:, t_local) = u_ac_down;
            U_ev_down_h(:, t_local) = u_ev_down;
            cost_down_h(t_local) = cost_down;
        end % 结束 parfor t_local

        % --- 将小时结果存入最终结果矩阵 ---
        U_ac_up_final(:, current_steps_indices) = U_ac_up_h;
        U_ev_up_final(:, current_steps_indices) = U_ev_up_h;
        cost_up_final(current_steps_indices) = cost_up_h;
        U_ac_down_final(:, current_steps_indices) = U_ac_down_h;
        U_ev_down_final(:, current_steps_indices) = U_ev_down_h;
        cost_down_final(current_steps_indices) = cost_down_h;

        % --- 清理内存 (可选) ---
        clear P_ac_up_hourly P_ev_up_hourly P_ac_down_hourly P_ev_down_hourly P_req_up_hourly P_req_down_hourly ...
              U_ac_up_h U_ev_up_h cost_up_h U_ac_down_h U_ev_down_h cost_down_h;

    end % 结束 for h (小时循环)

    toc(tic_loop);
    fprintf('分层优化完成。\n');

    %% 4. 计算最终的全局SDCI和Rho指标
    fprintf('计算最终全局指标...\n');
    % --- 优化后指标 ---
    % 优化读取：一次性读取所有需要的最终决策变量对应的功率
    fprintf('  读取优化后的AC上调功率...\n');
    P_ac_opt_up_all = double(m_individual.AC_Up_Individual); % 读取全量数据
    P_ac_opt_up = sum(U_ac_up_final .* P_ac_opt_up_all, 1)';   % 计算选中的功率
    clear P_ac_opt_up_all; % 及时释放内存
    fprintf('  读取优化后的EV上调功率...\n');
    P_ev_opt_up_all = double(m_individual.EV_Up_Individual);
    P_ev_opt_up = sum(U_ev_up_final .* P_ev_opt_up_all, 1)';
    clear P_ev_opt_up_all;
    fprintf('  读取优化后的AC下调功率...\n');
    P_ac_opt_down_all = double(abs(m_individual.AC_Down_Individual));
    P_ac_opt_down = sum(U_ac_down_final .* P_ac_opt_down_all, 1)';
    clear P_ac_opt_down_all;
    fprintf('  读取优化后的EV下调功率...\n');
    P_ev_opt_down_all = double(abs(m_individual.EV_Down_Individual));
    P_ev_opt_down = sum(U_ev_down_final .* P_ev_opt_down_all, 1)';
    clear P_ev_opt_down_all;

    n_ac_opt_up_t = sum(U_ac_up_final, 1)'; n_ev_opt_up_t = sum(U_ev_up_final, 1)';
    n_ac_opt_down_t = sum(U_ac_down_final, 1)'; n_ev_opt_down_t = sum(U_ev_down_final, 1)';

    avg_P_ac_opt_up = P_ac_opt_up ./ (n_ac_opt_up_t + eps_val);
    avg_P_ev_opt_up = P_ev_opt_up ./ (n_ev_opt_up_t + eps_val);
    avg_P_ac_opt_down = P_ac_opt_down ./ (n_ac_opt_down_t + eps_val);
    avg_P_ev_opt_down = P_ev_opt_down ./ (n_ev_opt_down_t + eps_val);

    SDCI_up_opt = ensureScalar(calculateSDCI(n_ac_opt_up_t, n_ev_opt_up_t, avg_P_ac_opt_up, avg_P_ev_opt_up));
    rho_up_opt  = ensureScalar(calculateSpearmanRho(n_ac_opt_up_t, avg_P_ac_opt_up, n_ev_opt_up_t, avg_P_ev_opt_up));
    SDCI_down_opt = ensureScalar(calculateSDCI(n_ac_opt_down_t, n_ev_opt_down_t, avg_P_ac_opt_down, avg_P_ev_opt_down));
    rho_down_opt  = ensureScalar(calculateSpearmanRho(n_ac_opt_down_t, avg_P_ac_opt_down, n_ev_opt_down_t, avg_P_ev_opt_down));

    % --- 优化前指标 (Raw) (同前) ---
    n_ac_raw_t = ones(T,1) * num_ac_total; n_ev_raw_t = ones(T,1) * num_ev_total;
    AC_Up_Aggregated_Raw = AC_Up_raw'; EV_Up_Aggregated_Raw = EV_Up_raw';
    AC_Down_Aggregated_Raw = AC_Down_raw'; EV_Down_Aggregated_Raw = EV_Down_raw';
    avg_P_ac_raw_up_t = AC_Up_Aggregated_Raw ./ (n_ac_raw_t + eps_val); avg_P_ev_raw_up_t = EV_Up_Aggregated_Raw ./ (n_ev_raw_t + eps_val);
    avg_P_ac_raw_down_t = AC_Down_Aggregated_Raw ./ (n_ac_raw_t + eps_val); avg_P_ev_raw_down_t = EV_Down_Aggregated_Raw ./ (n_ev_raw_t + eps_val);
    SDCI_up_raw = ensureScalar(calculateSDCI(n_ac_raw_t, n_ev_raw_t, avg_P_ac_raw_up_t, avg_P_ev_raw_up_t)); rho_up_raw  = ensureScalar(calculateSpearmanRho(n_ac_raw_t, avg_P_ac_raw_up_t, n_ev_raw_t, avg_P_ev_raw_up_t));
    SDCI_down_raw = ensureScalar(calculateSDCI(n_ac_raw_t, n_ev_raw_t, avg_P_ac_raw_down_t, avg_P_ev_raw_down_t)); rho_down_raw  = ensureScalar(calculateSpearmanRho(n_ac_raw_t, avg_P_ac_raw_down_t, n_ev_raw_t, avg_P_ev_raw_down_t));


    %% 5. 结果可视化
    fprintf('生成可视化图表...\n');
    Total_Up_Optimal_Agg = P_ac_opt_up + P_ev_opt_up;
    Total_Down_Optimal_Agg = P_ac_opt_down + P_ev_opt_down;
    % 修正：使用全量的 Raw 数据进行比较
    Total_Up_Aggregated_Raw_Full = results_agg.AC_Up(:) + results_agg.EV_Up(:);
    Total_Down_Aggregated_Raw_Full = abs(results_agg.AC_Down(:)) + abs(results_agg.EV_Down(:));

    time_axis = (1:T)';
    figure('Position', [100, 100, 1200, 800]);
    sgtitle(sprintf('VPP 分层优化结果 (全局 N=%d)', num_ac_total + num_ev_total), 'FontSize', 16);
    subplot(2,2,1); plot(time_axis, Total_Up_Aggregated_Raw_Full, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Raw Aggregated (Full)'); hold on; plot(time_axis, Total_Up_Optimal_Agg, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Optimized (Layered)'); plot(time_axis, P_grid_up_demand, 'k:', 'LineWidth', 1, 'DisplayName', 'Up Demand'); hold off; title('Up-Regulation Power'); xlabel('Time Step'); ylabel('Power (kW)'); legend('show', 'Location', 'best'); grid on;
    subplot(2,2,2); plot(time_axis, Total_Down_Aggregated_Raw_Full, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Raw Aggregated (Full)'); hold on; plot(time_axis, Total_Down_Optimal_Agg, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Optimized (Layered)'); plot(time_axis, P_grid_down_demand, 'k:', 'LineWidth', 1, 'DisplayName', 'Down Demand'); hold off; title('Down-Regulation Power'); xlabel('Time Step'); ylabel('Power (kW)'); legend('show', 'Location', 'best'); grid on;
    subplot(2,2,3); indicator_names_up = {'SDCI⁺', 'Spearman ρ⁺'}; values_raw_up = [SDCI_up_raw; rho_up_raw]; values_opt_up = [SDCI_up_opt; rho_up_opt]; bar_data_up = [values_raw_up, values_opt_up]; b_up = bar(bar_data_up); set(gca, 'XTickLabel', indicator_names_up); ylabel('Value'); ylim([-1.1, 1.1]); legend([b_up(1) b_up(2)], {'Raw', 'Optimized (Layered)'}, 'Location', 'northoutside', 'Orientation','horizontal'); title('Up Complementarity & Correlation'); grid on; add_bar_labels(b_up, bar_data_up);
    subplot(2,2,4); indicator_names_down = {'SDCI⁻', 'Spearman ρ⁻'}; values_raw_down = [SDCI_down_raw; rho_down_raw]; values_opt_down = [SDCI_down_opt; rho_down_opt]; bar_data_down = [values_raw_down, values_opt_down]; b_down = bar(bar_data_down); set(gca, 'XTickLabel', indicator_names_down); ylabel('Value'); ylim([-1.1, 1.1]); legend([b_down(1) b_down(2)],{'Raw', 'Optimized (Layered)'}, 'Location', 'northoutside', 'Orientation','horizontal'); title('Down Complementarity & Correlation'); grid on; add_bar_labels(b_down, bar_data_down);

    %% 6. 命令行输出汇总
    disp(' ');
    disp('=== 互补性与相关性指标对比 (分层优化 - 全量数据) ===');
    disp('【上调】'); fprintf('优化前 (Raw Aggregated): SDCI⁺ = %.4f, ρ⁺ = %.4f\n', SDCI_up_raw, rho_up_raw); fprintf('优化后 (Layered Opt): SDCI⁺ = %.4f, ρ⁺ = %.4f, 总成本 = %.2f\n', SDCI_up_opt, rho_up_opt, nansum(cost_up_final));
    disp(' ');
    disp('【下调】'); fprintf('优化前 (Raw Aggregated): SDCI⁻ = %.4f, ρ⁻ = %.4f\n', SDCI_down_raw, rho_down_raw); fprintf('优化后 (Layered Opt): SDCI⁻ = %.4f, ρ⁻ = %.4f, 总成本 = %.2f\n', SDCI_down_opt, rho_down_opt, nansum(cost_down_final));
    disp(' ');

end

%% -----------------------------------------------------------------------
%                       局部辅助函数
% ------------------------------------------------------------------------

%% 局部函数: 上层GA，优化每小时的参与设备数 (带功率约束)
function [n_ac_opt, n_ev_opt] = optimize_hourly_participation_GA_constrained(num_ac_total, num_ev_total, P_ac_hourly, P_ev_hourly, P_req_hourly, steps_per_hour, ga_options)

    nvars = 2; % [n_ac, n_ev]
    lb = [0, 0];
    ub = [num_ac_total, num_ev_total];
    intcon = 1:2;

    % 目标函数句柄 (包含功率约束检查)
    fitness_fcn = @(x) hourly_ga_fitness_constrained(x, P_ac_hourly, P_ev_hourly, P_req_hourly, steps_per_hour);

    % --- 设置 GA 选项 ---
    % 允许并行计算 (如果主循环是串行)
    ga_opts_default = optimoptions('ga', ...
        'PopulationSize', 50, ... % 可以适当增大
        'MaxGenerations', 80, ... % 可以适当增大
        'Display', 'off', ...
        'EliteCount', 3, ...
        'UseParallel', true, ... % 允许 GA 内部并行
        'ConstraintTolerance', 1e-3, ... % 放宽一点约束容忍度
        'FunctionTolerance', 1e-4);

    % 合并传入的选项和默认选项
    if nargin > 6 && ~isempty(ga_options)
        ga_opts_final = optimoptions(ga_opts_default, ga_options); % User options override defaults
    else
        ga_opts_final = ga_opts_default;
    end
    
    % --- 运行 GA ---
    try
        [x_opt, fval_opt, exitflag] = ga(fitness_fcn, nvars, [], [], [], [], lb, ub, [], intcon, ga_opts_final);
        if exitflag <= 0
             warning('小时级 GA 未找到最优解或提前终止 (exitflag=%d)。可能无法满足峰值需求。使用最大可用数量作为备选。', exitflag);
             % 备选方案：如果GA失败，返回最大数量以保证功率
             x_opt = ub;
        end
    catch ME_ga
        warning('小时级 GA 运行失败: %s。使用最大可用数量作为备选。', ME_ga.message);
        x_opt = ub; % 备选方案
    end


    n_ac_opt = round(x_opt(1));
    n_ev_opt = round(x_opt(2));
end

%% 局部函数: 上层GA的目标函数 (带功率约束惩罚)
function fitness = hourly_ga_fitness_constrained(x, P_ac_hourly, P_ev_hourly, P_req_hourly, steps_per_hour)
    n_ac = round(x(1));
    n_ev = round(x(2));
    
    % --- 1. 检查功率约束 ---
    power_penalty = 0;
    large_penalty_factor = 1e9; % 巨大的惩罚因子
    max_demand_hourly = max(P_req_hourly); % 本小时最大需求
    min_potential_power_hourly = inf; % 初始化为无穷大

    num_avail_ac = size(P_ac_hourly, 1);
    num_avail_ev = size(P_ev_hourly, 1);
    
    actual_n_ac = min(n_ac, num_avail_ac); % 防止索引越界
    actual_n_ev = min(n_ev, num_avail_ev);

    for t = 1:steps_per_hour
        current_potential_power = 0;
        if actual_n_ac > 0 && num_avail_ac > 0
            sorted_p_ac = sort(P_ac_hourly(:, t), 'descend');
            current_potential_power = current_potential_power + sum(sorted_p_ac(1:actual_n_ac));
        end
        if actual_n_ev > 0 && num_avail_ev > 0
            sorted_p_ev = sort(P_ev_hourly(:, t), 'descend');
            current_potential_power = current_potential_power + sum(sorted_p_ev(1:actual_n_ev));
        end
        min_potential_power_hourly = min(min_potential_power_hourly, current_potential_power);
    end

    power_shortage = max(0, max_demand_hourly - min_potential_power_hourly);
    if power_shortage > 1e-3 % 如果功率不足 (允许一点点误差)
        power_penalty = power_shortage * large_penalty_factor;
        fitness = power_penalty; % 直接返回大惩罚，无需计算后续指标
        return;
    end

    % --- 2. 如果功率约束满足，计算 SDCI 和 Rho 指标 ---
    p_ac_agg = zeros(steps_per_hour, 1);
    p_ev_agg = zeros(steps_per_hour, 1);
    if actual_n_ac > 0 && num_avail_ac > 0
        for t = 1:steps_per_hour
            sorted_p_ac = sort(P_ac_hourly(:, t), 'descend');
            p_ac_agg(t) = sum(sorted_p_ac(1:actual_n_ac));
        end
    end
     if actual_n_ev > 0 && num_avail_ev > 0
        for t = 1:steps_per_hour
            sorted_p_ev = sort(P_ev_hourly(:, t), 'descend');
            p_ev_agg(t) = sum(sorted_p_ev(1:actual_n_ev));
        end
    end
    
    n_ac_vec = ones(steps_per_hour, 1) * actual_n_ac; % 使用实际选取的数量
    n_ev_vec = ones(steps_per_hour, 1) * actual_n_ev;
    avg_p_ac = p_ac_agg ./ (actual_n_ac + 1e-9);
    avg_p_ev = p_ev_agg ./ (actual_n_ev + 1e-9);

    sdci_hourly = calculateSDCI(n_ac_vec, n_ev_vec, avg_p_ac, avg_p_ev);
    rho_hourly = calculateSpearmanRho(n_ac_vec, avg_p_ac, n_ev_vec, avg_p_ev);
    
    % 目标：最小化 SDCI 和 |Rho| (如果希望相关性尽量小，无论正负)
    w_sdci = 0.5;
    w_rho = 0.5;
    
    fitness = w_sdci * sdci_hourly + w_rho * abs(rho_hourly); % 或者 w_rho * (rho_hourly + 1) 如果希望负相关

    % 增加一个惩罚项，如果总功率为0 (虽然理论上功率约束已经处理了)
    if sum(p_ac_agg) + sum(p_ev_agg) < 1e-6
        fitness = fitness + 1e6;
    end
end


%% 局部函数: 下层贪心算法 (带数量约束) - solve_hourly_dispatch_greedy_with_count (同前)
function [u_ac_optimal, u_ev_optimal, total_cost] = solve_hourly_dispatch_greedy_with_count(p_ac, p_ev, c_ac, c_ev, P_req_t, n_ac_max, n_ev_max)
    % ... (函数体与上一版本完全相同) ...
    u_ac_optimal = zeros(length(p_ac), 1); u_ev_optimal = zeros(length(p_ev), 1);
    total_cost = 0;
    if P_req_t <= 0; return; end

    num_ac_avail = sum(p_ac > 1e-6); num_ev_avail = sum(p_ev > 1e-6);
    num_devices = num_ac_avail + num_ev_avail;
    if num_devices == 0; return; end

    device_data = struct('id', cell(num_devices, 1), 'type', cell(num_devices, 1), ...
                         'power', zeros(num_devices, 1), 'cost_per_kw', zeros(num_devices, 1), ...
                         'total_cost', zeros(num_devices, 1));
    idx = 1;
    for i = 1:length(p_ac); if p_ac(i) > 1e-6; device_data(idx).id=i; device_data(idx).type='AC'; device_data(idx).power=p_ac(i); device_data(idx).cost_per_kw=c_ac(i); device_data(idx).total_cost=p_ac(i)*c_ac(i); idx=idx+1; end; end
    for i = 1:length(p_ev); if p_ev(i) > 1e-6; device_data(idx).id=i; device_data(idx).type='EV'; device_data(idx).power=p_ev(i); device_data(idx).cost_per_kw=c_ev(i); device_data(idx).total_cost=p_ev(i)*c_ev(i); idx=idx+1; end; end

    [~, sorted_indices] = sort([device_data.cost_per_kw]);
    sorted_devices = device_data(sorted_indices);

    power_accumulated = 0; count_ac = 0; count_ev = 0;
    for i = 1:length(sorted_devices)
        if power_accumulated >= P_req_t; break; end
        device = sorted_devices(i);
        if strcmp(device.type, 'AC')
            if count_ac >= n_ac_max; continue; end
            count_ac = count_ac + 1; u_ac_optimal(device.id) = 1;
        else % EV
            if count_ev >= n_ev_max; continue; end
            count_ev = count_ev + 1; u_ev_optimal(device.id) = 1;
        end
        power_accumulated = power_accumulated + device.power;
        total_cost = total_cost + device.total_cost;
    end
    % 可以在这里加一个警告，如果最终 power_accumulated < P_req_t
     if power_accumulated < P_req_t - 1e-3 % 允许一点误差
          %fprintf('警告: 时间步 %d 贪心未能满足需求 %.2f kW (调度 %.2f kW) n_ac=%d/%d, n_ev=%d/%d\n', t_global, P_req_t, power_accumulated, count_ac, n_ac_max, count_ev, n_ev_max);
          % 在 parfor 中无法直接访问 t_global，可以在主循环结束后统一检查
     end
end


%% 局部函数: 其他 (generate_demand, ensureScalar, add_bar_labels - 同前)
function P_demand = generate_demand(results_agg, field_name, T, avg_raw_power, factors, default_avg)
    if isfield(results_agg, field_name) && ~isempty(results_agg.(field_name)) && length(results_agg.(field_name)) == T
        P_demand = results_agg.(field_name)(:);
    else
        if isnan(avg_raw_power) || avg_raw_power <= 0; avg_raw_power = default_avg; end
        P_demand = avg_raw_power * (factors(1) + (factors(2)-factors(1)) * rand(T,1));
        fprintf('未在聚合文件中找到 "%s"，已生成示例需求。\n', field_name);
    end
end

function val = ensureScalar(inputVal)
    if isscalar(inputVal); val = inputVal;
    elseif isempty(inputVal); val = NaN;
    else; val = mean(inputVal(:), 'omitnan'); end
end

function add_bar_labels(bar_handles, bar_data)
    for k_bar = 1:size(bar_data, 1)
        for j_bar = 1:size(bar_data, 2)
             try % 增加错误处理
                 if ~ishandle(bar_handles(j_bar)); continue; end % 检查句柄是否有效
                 if k_bar > length(bar_handles(j_bar).XData); continue; end % 检查索引
                 if isnan(bar_data(k_bar, j_bar)); continue; end % 跳过 NaN
                 
                text(bar_handles(j_bar).XData(k_bar) + bar_handles(j_bar).XOffset, ...
                     bar_data(k_bar, j_bar), sprintf('%.3f', bar_data(k_bar, j_bar)), ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                     'FontSize', 8, 'Color', 'k');
             catch ME_label
                 % fprintf('添加标签时出错: %s\n', ME_label.message); % 调试时可以取消注释
             end
        end
    end
end