% main_vpp_optimizer_multi_objective.m
% 解决VPP大规模调度优化问题的主脚本 (v5 - 多目标: 成本 + 互补性[SDCI & Rho])
clear; close all; clc;
tic;

%% 1. 定义全局输入
% --- 用户必须修改的参数 ---
chunk_results_dir = 'chunk_results'; 

selected_mat_files = {
    'results_chunk_1.mat',
    'results_chunk_2.mat',
    'results_chunk_3.mat' 
};

directions_to_run = {'Up', 'Down'};

% *** 多目标优化权重 ***
% W_cost + W_complementarity 应该等于 1.0
% W_cost = 1.0 -> 纯成本优化 (忽略 SDCI 和 Rho)
% W_cost = 0.0 -> 纯互补性优化 (忽略成本)
W_cost = 0.7; 
W_complementarity = 0.3; % 这个权重将同时优化 SDCI 和 Rho
% -------------------------

fprintf('开始执行VPP大规模优化 (多目标启发式)...\n');
fprintf('优化权重: Cost=%.2f, Complementarity(SDCI+Rho)=%.2f\n', ...
        W_cost, W_complementarity);

%% 2. 调节方向主循环
for dir_idx = 1:length(directions_to_run)
    
    output_direction = directions_to_run{dir_idx};
    
    fprintf('\n====================================================\n');
    fprintf('== [正在处理方向: %s 调节] ==\n', output_direction);
    fprintf('====================================================\n');
    
    %% 2.1 加载仿真数据
    fprintf('正在从 %s 加载并聚合 %s 潜力数据...\n', ...
            chunk_results_dir, output_direction);
    SimData = load_simulation_data_from_chunks(...
        chunk_results_dir, selected_mat_files, output_direction);

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
    p_total_optimal = zeros(1, SimData.T); 

    p_all_devices = [SimData.p_AC; SimData.p_EV];         
    c_all_devices = [OptInputs.DeviceCosts.c_AC; ...
                     OptInputs.DeviceCosts.c_EV];         
    device_nodes = OptInputs.Network.DeviceNodes;         

    % *** 为多目标启发式准备额外参数 ***
    % 1. 设备类型标识 (1=AC, 2=EV)
    device_types = [ones(SimData.nAC, 1); 2*ones(SimData.nEV, 1)];
    % 2. 聚合潜力曲线 (1 x T)
    P_AC_agg_total_t = sum(SimData.p_AC, 1);
    P_EV_agg_total_t = sum(SimData.p_EV, 1);
    % 3. 避免除零
    P_Total_agg_t = P_AC_agg_total_t + P_EV_agg_total_t;
    P_Total_agg_t(P_Total_agg_t == 0) = 1; % 避免0/0
    % 4. 计算AC资源占比 (1 x T)
    Ratio_AC_t = P_AC_agg_total_t ./ P_Total_agg_t;
    % *** 结束准备 ***
    
    parfor t = 1:SimData.T
        p_t = p_all_devices(:, t);
        P_req_t = OptInputs.P_req(t);
        Network_t = struct(...
            'PTDF', OptInputs.Network.PTDF, ...
            'LineLimits', OptInputs.Network.LineLimits, ...
            'BaseFlow', OptInputs.Network.BaseFlow(:, t), ...
            'DeviceNodes', device_nodes ...
        );
        
        % *** 修改：调用新的多目标求解器 ***
        [u_opt_t, P_total_t, cost_t, is_feasible] = ...
            solve_times_step_greedy_multi_obj(...
                p_t, c_all_devices, P_req_t, Network_t, ...
                W_cost, W_complementarity, Ratio_AC_t(t), device_types ...
            );
        % *** 结束修改 ***
            
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
    fprintf('优化权重: Cost=%.2f, Complementarity=%.2f\n', W_cost, W_complementarity);
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
        time_axis = 1:SimData.T;
        P_ori_total = Metrics_ori.P_AC_agg + Metrics_ori.P_EV_agg; 
        P_opt_total = p_total_optimal'; 
        P_demand = OptInputs.P_req'; 

        figure('Name', sprintf('VPP Optimization Results (%s, W_comp=%.1f)', output_direction, W_complementarity), ...
               'Position', [100 100 1200 600]);
        hold on;
        
        plot(time_axis, P_ori_total, 'r--', 'LineWidth', 2, 'DisplayName', '基准总潜力 (Ori)');
        plot(time_axis, P_opt_total, 'b-', 'LineWidth', 2, 'DisplayName', '优化后调度 (Opt)');
        plot(time_axis, P_demand, 'k:', 'LineWidth', 2.5, 'DisplayName', '电网需求 (P_{req})');
        
        xlabel('时间步 (Time Step)');
        ylabel(sprintf('VPP %s 调节功率 (kW)', output_direction));
        title_str = sprintf('VPP 优化调度 (%s, W_{cost}=%.1f, W_{comp}=%.1f) vs 需求', ...
                            output_direction, W_cost, W_complementarity);
        title(title_str);
        legend('show', 'Location', 'best');
        grid on;
        hold off;
        
        output_plot_filename = sprintf('vpp_optimization_results_%s_Wcomp_%.1f.png', ...
                                       output_direction, W_complementarity);
        saveas(gcf, output_plot_filename);
        fprintf('结果图已保存为: %s\n', output_plot_filename);
        
    catch ME_plot
        fprintf('*** 绘制 [%s] 结果图时出错: %s ***\n', output_direction, ME_plot.message);
    end

end % 结束调节方向主循环

fprintf('\n所有调节方向处理完毕。\n');
toc;