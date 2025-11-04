function run_hourly_optimization_greedy_full()
    %% 内存优化的逐时段调度脚本 (版本 3: 贪心算法 + 全量数据)
    clc; close all; clear;

    % --- 1. 加载数据 (优化) ---
    aggregate_file = 'aggregate_results_slow_synced.mat';
    individual_file = 'individual_results_slow_synced.mat';

    fprintf('正在加载聚合数据: %s\n', aggregate_file);
    if exist(aggregate_file, 'file')
        agg_data = load(aggregate_file);
        if ~isfield(agg_data, 'results'); error('文件 %s 中未找到 "results" 结构体。', aggregate_file); end
        results_agg = agg_data.results;
    else; error('聚合数据文件 %s 未找到。', aggregate_file); end

    fprintf('创建个体数据文件对象 (matfile): %s\n', individual_file);
    if exist(individual_file, 'file')
        m_individual = matfile(individual_file);
    else; error('个体数据文件 %s 未找到。', individual_file); end

    % --- 从 matfile 对象中读取维度信息 ---
    fprintf('正在读取维度信息...\n');
    try
        s = whos(m_individual);
        ac_up_info = s(strcmp({s.name}, 'AC_Up_Individual'));
        ev_up_info = s(strcmp({s.name}, 'EV_Up_Individual'));
        num_ac_total = ac_up_info.size(1);
        num_ev_total = ev_up_info.size(1);
        T = ac_up_info.size(2);
        if T ~= ev_up_info.size(2) || T ~= length(results_agg.AC_Up); warning('个体数据和聚合数据的时间步长不匹配。'); end
    catch ME; error('从 %s 读取维度信息时出错: %s', individual_file, ME.message); end

    AC_Up_raw = results_agg.AC_Up(:)'; AC_Down_raw = abs(results_agg.AC_Down(:)');
    EV_Up_raw = results_agg.EV_Up(:)'; EV_Down_raw = abs(results_agg.EV_Down(:)');

    fprintf('数据维度加载完成。\n');
    fprintf('总细分时间步数 T = %d\n', T);
    fprintf('总空调数 num_ac_total = %d\n', num_ac_total);
    fprintf('总电动汽车数 num_ev_total = %d\n', num_ev_total);

    % --- 2. 定义仿真参数 ---
    P_grid_up_demand = generate_demand(results_agg, 'P_grid_up_regulation_demand', T, (sum(AC_Up_raw) + sum(EV_Up_raw)) / T, [0.2, 0.5], 100);
    P_grid_down_demand = generate_demand(results_agg, 'P_grid_down_regulation_demand', T, (sum(AC_Down_raw) + sum(EV_Down_raw)) / T, [0.15, 0.4], 80);
    c_ac_up = ones(num_ac_total, 1) * 0.05; c_ev_up = ones(num_ev_total, 1) * 0.04;
    c_ac_down = ones(num_ac_total, 1) * 0.03; c_ev_down = ones(num_ev_total, 1) * 0.02;
    eps_val = 1e-6;

    %% 3. 逐时段优化 (核心修改：使用贪心算法处理全量数据)

    fprintf('\n--- 正在准备对所有 %d 台设备进行逐时段贪心优化 ---\n', num_ac_total + num_ev_total);
    
    % 使用全量索引
    ac_indices = 1:num_ac_total;
    ev_indices = 1:num_ev_total;
    num_ac_opt = num_ac_total; % 优化全部 AC
    num_ev_opt = num_ev_total; % 优化全部 EV

    % 成本向量 (全量)
    c_ac_up_opt = c_ac_up; c_ev_up_opt = c_ev_up;
    c_ac_down_opt = c_ac_down; c_ev_down_opt = c_ev_down;

    % 初始化结果存储
    u_ac_opt_up_t = zeros(num_ac_opt, T); u_ev_opt_up_t = zeros(num_ev_opt, T);
    optimal_cost_up_t = zeros(1, T);
    u_ac_opt_down_t = zeros(num_ac_opt, T); u_ev_opt_down_t = zeros(num_ev_opt, T);
    optimal_cost_down_t = zeros(1, T);
    P_ac_dispatched_opt_up_t = zeros(1,T); P_ev_dispatched_opt_up_t = zeros(1,T);
    P_ac_dispatched_opt_down_t = zeros(1,T); P_ev_dispatched_opt_down_t = zeros(1,T);

    fprintf('\n开始逐时段贪心优化 (共 %d 步)...\n', T);
    tic_loop = tic;

    % 启动并行池
    pool = gcp('nocreate');
    

    % --- 临时存储单元，用于 parfor ---
    tmp_u_ac_up = cell(T, 1); tmp_u_ev_up = cell(T, 1); tmp_cost_up = zeros(T, 1);
    tmp_p_ac_up = zeros(T, 1); tmp_p_ev_up = zeros(T, 1);
    tmp_u_ac_down = cell(T, 1); tmp_u_ev_down = cell(T, 1); tmp_cost_down = zeros(T, 1);
    tmp_p_ac_down = zeros(T, 1); tmp_p_ev_down = zeros(T, 1);

    % 使用 parfor
    parfor t = 1:T
        if mod(t, 10) == 1 % 每 10 步打印一次
            fprintf('  并行处理时间步 %d 附近...\n', t);
        end

        % --- 上调优化 ---
        % 读取当前时间步 t 的全量数据
        p_ac_t_up = double(m_individual.AC_Up_Individual(:, t));
        p_ev_t_up = double(m_individual.EV_Up_Individual(:, t));

        % !!! 调用贪心算法 !!!
        [u_ac_up, u_ev_up, cost_val_up] = solve_hourly_dispatch_greedy(num_ac_opt, num_ev_opt, p_ac_t_up, p_ev_t_up, c_ac_up_opt, c_ev_up_opt, P_grid_up_demand(t));

        % 存储结果 (贪心算法总是能找到一个解，除非总容量不足)
        tmp_u_ac_up{t} = u_ac_up; tmp_u_ev_up{t} = u_ev_up;
        tmp_cost_up(t) = cost_val_up;
        tmp_p_ac_up(t) = sum(u_ac_up .* p_ac_t_up);
        tmp_p_ev_up(t) = sum(u_ev_up .* p_ev_t_up);


        % --- 下调优化 ---
        % 读取当前时间步 t 的全量数据
        p_ac_t_down = double(abs(m_individual.AC_Down_Individual(:, t)));
        p_ev_t_down = double(abs(m_individual.EV_Down_Individual(:, t)));

        % !!! 调用贪心算法 !!!
        [u_ac_down, u_ev_down, cost_val_down] = solve_hourly_dispatch_greedy(num_ac_opt, num_ev_opt, p_ac_t_down, p_ev_t_down, c_ac_down_opt, c_ev_down_opt, P_grid_down_demand(t));

        % 存储结果
        tmp_u_ac_down{t} = u_ac_down; tmp_u_ev_down{t} = u_ev_down;
        tmp_cost_down(t) = cost_val_down;
        tmp_p_ac_down(t) = sum(u_ac_down .* p_ac_t_down);
        tmp_p_ev_down(t) = sum(u_ev_down .* p_ev_t_down);
    end

    % 将 parfor 的 cell 结果复制回矩阵
    for t = 1:T
        if ~isempty(tmp_u_ac_up{t}); u_ac_opt_up_t(:, t) = tmp_u_ac_up{t}; end
        if ~isempty(tmp_u_ev_up{t}); u_ev_opt_up_t(:, t) = tmp_u_ev_up{t}; end
        optimal_cost_up_t(t) = tmp_cost_up(t);
        P_ac_dispatched_opt_up_t(t) = tmp_p_ac_up(t);
        P_ev_dispatched_opt_up_t(t) = tmp_p_ev_up(t);

        if ~isempty(tmp_u_ac_down{t}); u_ac_opt_down_t(:, t) = tmp_u_ac_down{t}; end
        if ~isempty(tmp_u_ev_down{t}); u_ev_opt_down_t(:, t) = tmp_u_ev_down{t}; end
        optimal_cost_down_t(t) = tmp_cost_down(t);
        P_ac_dispatched_opt_down_t(t) = tmp_p_ac_down(t);
        P_ev_dispatched_opt_down_t(t) = tmp_p_ev_down(t);
    end

    toc(tic_loop);
    fprintf('逐时段贪心优化完成。\n');


    %% 4. 计算优化前后的SDCI和Rho指标 (基于全量数据)
    fprintf('计算指标 (基于全量数据)...\n');
    n_ac_raw_t = ones(T,1) * num_ac_total; n_ev_raw_t = ones(T,1) * num_ev_total;
    if num_ac_total == 0; n_ac_raw_t = zeros(T,1); end; if num_ev_total == 0; n_ev_raw_t = zeros(T,1); end

    % 直接使用加载的聚合数据作为 "Raw"
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
    n_ac_opt_up_count_t = sum(u_ac_opt_up_t, 1)'; n_ev_opt_up_count_t = sum(u_ev_opt_up_t, 1)';
    n_ac_opt_down_count_t = sum(u_ac_opt_down_t, 1)'; n_ev_opt_down_count_t = sum(u_ev_opt_down_t, 1)';
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
    sgtitle(sprintf('VPP 逐时优化结果 (贪心算法 N_{total}=%d)', num_ac_total + num_ev_total), 'FontSize', 16);
    subplot(2,2,1); plot(time_axis, Total_Up_Aggregated_Raw, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Raw Aggregated'); hold on; plot(time_axis, Total_Up_Optimal_Agg, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Optimized (Greedy)'); plot(time_axis, P_grid_up_demand, 'k:', 'LineWidth', 1, 'DisplayName', 'Up-Regulation Demand'); hold off; title('Up-Regulation Power Comparison'); xlabel('Time Step'); ylabel('Power (kW)'); legend('show', 'Location', 'best'); grid on;
    subplot(2,2,2); plot(time_axis, Total_Down_Aggregated_Raw, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Raw Aggregated'); hold on; plot(time_axis, Total_Down_Optimal_Agg, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Optimized (Greedy)'); plot(time_axis, P_grid_down_demand, 'k:', 'LineWidth', 1, 'DisplayName', 'Down-Regulation Demand'); hold off; title('Down-Regulation Power Comparison'); xlabel('Time Step'); ylabel('Power (kW)'); legend('show', 'Location', 'best'); grid on;
    subplot(2,2,3); indicator_names_up = {'SDCI⁺', 'Spearman ρ⁺'}; values_raw_up = [SDCI_up_raw; rho_up_raw]; values_opt_up = [SDCI_up_opt; rho_up_opt]; bar_data_up = [values_raw_up, values_opt_up]; b_up = bar(bar_data_up); set(gca, 'XTickLabel', indicator_names_up); ylabel('Indicator Value'); ylim([-1.1, 1.1]); legend([b_up(1) b_up(2)], {'Raw Aggregated', 'Optimized (Greedy)'}, 'Location', 'northoutside', 'Orientation','horizontal'); title('Up-Regulation Complementarity & Correlation'); grid on; add_bar_labels(b_up, bar_data_up);
    subplot(2,2,4); indicator_names_down = {'SDCI⁻', 'Spearman ρ⁻'}; values_raw_down = [SDCI_down_raw; rho_down_raw]; values_opt_down = [SDCI_down_opt; rho_down_opt]; bar_data_down = [values_raw_down, values_opt_down]; b_down = bar(bar_data_down); set(gca, 'XTickLabel', indicator_names_down); ylabel('Indicator Value'); ylim([-1.1, 1.1]); legend([b_down(1) b_down(2)],{'Raw Aggregated', 'Optimized (Greedy)'}, 'Location', 'northoutside', 'Orientation','horizontal'); title('Down-Regulation Complementarity & Correlation'); grid on; add_bar_labels(b_down, bar_data_down);

    %% 6. 命令行输出汇总
    disp(' ');
    disp('=== 互补性与相关性指标对比 (Greedy Optimization - Full Data) ===');
    disp('【上调】'); fprintf('优化前 (Raw Aggregated): SDCI⁺ = %.4f, ρ⁺ = %.4f\n', SDCI_up_raw, rho_up_raw); fprintf('优化后 (Greedy Dispatch): SDCI⁺ = %.4f, ρ⁺ = %.4f, 总成本 = %.2f\n', SDCI_up_opt, rho_up_opt, nansum(optimal_cost_up_t));
    disp(' ');
    disp('【下调】'); fprintf('优化前 (Raw Aggregated): SDCI⁻ = %.4f, ρ⁻ = %.4f\n', SDCI_down_raw, rho_down_raw); fprintf('优化后 (Greedy Dispatch): SDCI⁻ = %.4f, ρ⁻ = %.4f, 总成本 = %.2f\n', SDCI_down_opt, rho_down_opt, nansum(optimal_cost_down_t));
    disp(' ');

end % 结束主函数

%% -----------------------------------------------------------------------
%                       局部辅助函数
% ------------------------------------------------------------------------

%% 局部函数: solve_hourly_dispatch_greedy (!!! 新增 !!!)
function [u_ac_optimal, u_ev_optimal, total_cost] = solve_hourly_dispatch_greedy(num_ac, num_ev, p_ac, p_ev, c_ac, c_ev, P_req_t)
    % 使用贪心算法解决单时段调度问题 (最小化成本)
    
    u_ac_optimal = zeros(num_ac, 1);
    u_ev_optimal = zeros(num_ev, 1);
    total_cost = 0;
    
    if P_req_t <= 0 || (num_ac == 0 && num_ev == 0)
        return; % 无需求或无设备
    end
    
    % --- 1. 准备设备列表 ---
    num_devices = num_ac + num_ev;
    device_data = struct('id', cell(num_devices, 1), ...
                         'type', cell(num_devices, 1), ...
                         'power', zeros(num_devices, 1), ...
                         'cost_per_kw', zeros(num_devices, 1), ... % 单位功率成本
                         'total_cost', zeros(num_devices, 1)); % 参与的总成本
                         
    idx = 1;
    % 添加 AC 设备信息
    for i = 1:num_ac
        if p_ac(i) > 1e-6 % 只考虑有实际能力的设备
            device_data(idx).id = i;
            device_data(idx).type = 'AC';
            device_data(idx).power = p_ac(i);
            device_data(idx).cost_per_kw = c_ac(i); % 使用单位成本 c_i 作为排序依据
            device_data(idx).total_cost = p_ac(i) * c_ac(i);
            idx = idx + 1;
        end
    end
    % 添加 EV 设备信息
    for i = 1:num_ev
        if p_ev(i) > 1e-6
            device_data(idx).id = i;
            device_data(idx).type = 'EV';
            device_data(idx).power = p_ev(i);
            device_data(idx).cost_per_kw = c_ev(i);
            device_data(idx).total_cost = p_ev(i) * c_ev(i);
            idx = idx + 1;
        end
    end
    
    % 移除未使用的结构体元素
    device_data(idx:end) = [];
    
    if isempty(device_data)
        return; % 没有可用的设备
    end
    
    % --- 2. 按成本效益排序 ---
    % 按 cost_per_kw (即 c_i) 升序排序
    [~, sorted_indices] = sort([device_data.cost_per_kw]);
    sorted_devices = device_data(sorted_indices);
    
    % --- 3. 贪心选择 ---
    power_accumulated = 0;
    
    for i = 1:length(sorted_devices)
        if power_accumulated >= P_req_t
            break; % 需求已满足
        end
        
        device = sorted_devices(i);
        
        % 选择该设备
        if strcmp(device.type, 'AC')
            u_ac_optimal(device.id) = 1;
        else % EV
            u_ev_optimal(device.id) = 1;
        end
        
        power_accumulated = power_accumulated + device.power;
        total_cost = total_cost + device.total_cost;
    end
    
     % (可选) 检查是否满足需求
     if power_accumulated < P_req_t
         % warning('时间步 %d: 贪心算法未能完全满足需求 %.2f kW (已调度 %.2f kW)', t, P_req_t, power_accumulated);
         % 在 parfor 中不方便打印 t，可以在外部检查
     end

end

%% 局部函数: generate_demand (保持不变)
function P_demand = generate_demand(results_agg, field_name, T, avg_raw_power, factors, default_avg)
    if isfield(results_agg, field_name) && ~isempty(results_agg.(field_name)) && length(results_agg.(field_name)) == T
        P_demand = results_agg.(field_name)(:);
        % fprintf('已从聚合文件加载 "%s"。\n', field_name); % 减少打印
    else
        if isnan(avg_raw_power) || avg_raw_power <= 0; avg_raw_power = default_avg; end
        P_demand = avg_raw_power * (factors(1) + (factors(2)-factors(1)) * rand(T,1));
        fprintf('未在聚合文件中找到 "%s"，已生成示例需求。\n', field_name);
    end
end

%% 局部函数: ensureScalar (保持不变)
function val = ensureScalar(inputVal)
    if isscalar(inputVal); val = inputVal;
    elseif isempty(inputVal); val = NaN;
    else; val = mean(inputVal(:), 'omitnan'); end
end

%% 局部函数: add_bar_labels (保持不变)
function add_bar_labels(bar_handles, bar_data)
    for k_bar = 1:size(bar_data, 1)
        for j_bar = 1:size(bar_data, 2)
             try
                text(bar_handles(j_bar).XData(k_bar) + bar_handles(j_bar).XOffset, ...
                     bar_data(k_bar, j_bar), ...
                     sprintf('%.3f', bar_data(k_bar, j_bar)), ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                     'FontSize', 8, 'Color', 'k');
             catch; end % Ignore errors
        end
    end
end