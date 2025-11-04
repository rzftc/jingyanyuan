% main_vpp_optimizer.m
% 解决VPP大规模调度优化问题的主脚本 (v4 - 支持上调和下调)
clear; close all; clc;
tic;

%% 1. 定义全局输入
% --- 用户必须修改的参数 ---
chunk_results_dir = 'chunk_results'; % ac_ev_simulation_block.m 的输出目录

% *** 在此处手动指定要加载的MAT文件列表 ***
selected_mat_files = {
    'results_chunk_1.mat',
    'results_chunk_2.mat',
    'results_chunk_3.mat' 
    % ... (添加任意数量您想要加载的文件)
};

% *** 新增：定义要运行的调节方向 ***
directions_to_run = {'Up', 'Down'};
% -------------------------

fprintf('开始执行VPP大规模优化...\n');
fprintf('将依次处理 %d 个调节方向: %s\n', ...
        length(directions_to_run), strjoin(directions_to_run, ', '));

%% 2. 调节方向主循环
for dir_idx = 1:length(directions_to_run)
    
    output_direction = directions_to_run{dir_idx};
    
    fprintf('\n====================================================\n');
    fprintf('== [正在处理方向: %s 调节] ==\n', output_direction);
    fprintf('====================================================\n');
    
    %% 2.1 加载仿真数据 (Up 或 Down)
    fprintf('正在从 %s 加载并聚合 %s 潜力数据...\n', ...
            chunk_results_dir, output_direction);

    SimData = load_simulation_data_from_chunks(...
        chunk_results_dir, ...
        selected_mat_files, ...
        output_direction ...
    );

    fprintf('加载完成: %d 台AC, %d 台EV, %d 个时间步。\n', ...
            SimData.nAC, SimData.nEV, SimData.T);

    %% 2.2 加载优化输入 (动态生成)
    fprintf('正在为 %s 调节创建动态虚拟电网/成本输入...\n', output_direction);
    OptInputs = create_dummy_optimization_inputs(SimData);
    fprintf('虚拟数据创建完成。\n');

    %% 2.3 计算基准互补性指标 (Ori)
    fprintf('正在计算 %s 基准 (Ori) 互补性指标...\n', output_direction);
    u_baseline_ac = SimData.p_AC > 0; % nAC x T
    u_baseline_ev = SimData.p_EV > 0; % nEV x T
    u_baseline = [u_baseline_ac; u_baseline_ev];
    Metrics_ori = calculate_solution_metrics(u_baseline, SimData, 'ori');

    %% 2.4 执行逐时贪心优化
    fprintf('开始为 %s 调节进行逐时贪心优化 (共 %d 步)...\n', ...
            output_direction, SimData.T);
    u_optimal = zeros(SimData.N_total, SimData.T, 'logical'); 
    total_cost = 0;
    p_total_optimal = zeros(1, SimData.T); % (1 x T)

    p_all_devices = [SimData.p_AC; SimData.p_EV];         
    c_all_devices = [OptInputs.DeviceCosts.c_AC; ...
                     OptInputs.DeviceCosts.c_EV];         
    device_nodes = OptInputs.Network.DeviceNodes;         

    parfor t = 1:SimData.T
        % (在 parfor 循环中减少打印)
        % if mod(t, 24) == 0
        %     fprintf('  [%s] 正在处理时间步 %d / %d ...\n', output_direction, t, SimData.T);
        % end
        
        p_t = p_all_devices(:, t);
        P_req_t = OptInputs.P_req(t);
        Network_t = struct(...
            'PTDF', OptInputs.Network.PTDF, ...
            'LineLimits', OptInputs.Network.LineLimits, ...
            'BaseFlow', OptInputs.Network.BaseFlow(:, t), ...
            'DeviceNodes', device_nodes ...
        );
        
        [u_opt_t, P_total_t, cost_t, is_feasible] = ...
            solve_times_step_greedy(p_t, c_all_devices, P_req_t, Network_t);
            
        if ~is_feasible
            fprintf('警告 [%s]: 时间步 %d 无法满足需求 P_req=%.2f. 仅提供了 %.2f\n', ...
                    output_direction, t, P_req_t, P_total_t);
        end
        
        u_optimal(:, t) = u_opt_t;
        p_total_optimal(t) = P_total_t;
        total_cost = total_cost + cost_t;
    end

    fprintf('[%s] 贪心优化完成。总调节成本: %.2f\n', output_direction, total_cost);

    %% 2.5 计算优化后互补性指标 (Opt)
    fprintf('正在计算 %s 优化后 (Opt) 互补性指标...\n', output_direction);
    Metrics_opt = calculate_solution_metrics(u_optimal, SimData, 'opt');

    %% 2.6 结果对比与报告
    fprintf('\n================ [%s] 优化结果报告 ================\n', output_direction);
    fprintf('目标: 最小化总成本: %.2f (元)\n', total_cost);
    fprintf('约束 1 (满足需求): 见逐时警告。\n');
    fprintf('约束 4 (网络潮流): 在贪心算法中已严格执行。\n');
    fprintf('\n--- 互补性指标对比 (目标: Opt < Ori) ---\n');
    fprintf('约束 2 (SDCI): \n');
    fprintf('  - SDCI_ori: %.4f (基准)\n', Metrics_ori.SDCI);
    fprintf('  - SDCI_opt: %.4f (优化后)\n', Metrics_opt.SDCI);
    if Metrics_opt.SDCI < Metrics_ori.SDCI
        fprintf('  -> 结果: 成功改善 (%.4f < %.4f)\n', Metrics_opt.SDCI, Metrics_ori.SDCI);
    else
        fprintf('  -> 结果: 未改善 (%.4f >= %.4f)\n', Metrics_opt.SDCI, Metrics_ori.SDCI);
    end

    fprintf('\n约束 3 (Spearman Rho): \n');
    fprintf('  - Rho_ori: %.4f (基准)\n', Metrics_ori.Rho);
    fprintf('  - Rho_opt: %.4f (优化后)\n', Metrics_opt.Rho);
    if abs(Metrics_opt.Rho) < abs(Metrics_ori.Rho)
        fprintf('  -> 结果: 成功改善 (|%.4f| < |%.4f|)\n', Metrics_opt.Rho, Metrics_ori.Rho);
    else
        fprintf('  -> 结果: 未改善 (|%.4f| >= |%.4f|)\n', Metrics_opt.Rho, Metrics_ori.Rho);
    end
    fprintf('=======================================================\n');

    %% 2.7 绘制结果对比图
    fprintf('正在为 %s 调节生成优化结果对比图...\n', output_direction);

    try
        % 1. 准备数据
        time_axis = 1:SimData.T;
        
        % 基准潜力 (T x 1)
        P_ori_total = Metrics_ori.P_AC_agg + Metrics_ori.P_EV_agg; 
        
        % 优化后调度 (T x 1) - p_total_optimal 是 (1 x T)，需转置
        P_opt_total = p_total_optimal'; 
        
        % 电网需求 (T x 1) - OptInputs.P_req 是 (1 x T)，需转置
        P_demand = OptInputs.P_req'; 

        % 2. 创建图形
        figure('Name', sprintf('VPP Optimization Results (%s)', output_direction), ...
               'Position', [100 100 1200 600]);
        hold on;
        
        % 绘制曲线
        plot(time_axis, P_ori_total, 'r--', 'LineWidth', 2, 'DisplayName', '基准总潜力 (Ori)');
        plot(time_axis, P_opt_total, 'b-', 'LineWidth', 2, 'DisplayName', '优化后调度 (Opt)');
        plot(time_axis, P_demand, 'k:', 'LineWidth', 2.5, 'DisplayName', '电网需求 (P_{req})');
        
        % 3. 添加标签和图例
        xlabel('时间步 (Time Step)');
        ylabel(sprintf('VPP %s 调节功率 (kW)', output_direction));
        title(sprintf('VPP 优化调度 (%s) vs 基准潜力 vs 电网需求', output_direction));
        legend('show', 'Location', 'best');
        grid on;
        hold off;
        
        % 4. 保存图形
        output_plot_filename = sprintf('vpp_optimization_results_%s.png', output_direction);
        saveas(gcf, output_plot_filename);
        fprintf('结果图已保存为: %s\n', output_plot_filename);
        
    catch ME_plot
        fprintf('*** 绘制 [%s] 结果图时出错: %s ***\n', output_direction, ME_plot.message);
    end

end % 结束调节方向主循环

fprintf('\n所有调节方向处理完毕。\n');
toc;