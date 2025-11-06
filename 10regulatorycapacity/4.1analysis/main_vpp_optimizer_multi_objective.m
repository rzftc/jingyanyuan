% main_vpp_optimizer_multi_objective.m
% 解决VPP大规模调度优化问题的主脚本 (v6 - 多目标 + 聚合数据读取 + 保留完整分析)
clear; close all; clc;
tic;

%% 1. 全局配置
chunk_results_dir = 'chunk_results'; 
selected_files = {
    'results_chunk_1.mat',
    'results_chunk_2.mat',
    'results_chunk_3.mat'
};
run_directions = {'Up', 'Down'}; 

W_cost = 0.6;              % 成本权重
W_complementarity = 0.4;   % 互补性权重 (SDCI & Rho)

fprintf('启动 VPP 多目标优化 (Cost=%.1f, Comp=%.1f)...\n', W_cost, W_complementarity);

%% 2. 主循环：按方向依次执行
for d_idx = 1:length(run_directions)
    this_dir = run_directions{d_idx};
    fprintf('\n>>> 开始处理 %s (调节) 方向 >>>\n', upper(this_dir));

    % --- 2.1 加载数据 ---
    SimData = load_simulation_data_from_chunks(chunk_results_dir, selected_files, this_dir);
    fprintf('数据加载完毕: AC=%d, EV=%d, T=%d\n', SimData.nAC, SimData.nEV, SimData.T);

    % --- 2.2 生成配套输入 ---
    OptIn = create_dummy_optimization_inputs(SimData);

    % --- 2.3 准备优化变量 ---
    % [修改] 使用预加载聚合数据的【幅度】来计算资源占比
    P_AC_mag = abs(SimData.P_AC_agg_loaded);
    P_EV_mag = abs(SimData.P_EV_agg_loaded);
    P_Total_mag = P_AC_mag + P_EV_mag;
    
    Ratio_AC_t = zeros(1, SimData.T);
    valid_t = P_Total_mag > 1e-9;
    Ratio_AC_t(valid_t) = P_AC_mag(valid_t) ./ P_Total_mag(valid_t);

    dev_types = [ones(SimData.nAC,1); 2*ones(SimData.nEV,1)]; % 1=AC, 2=EV
    p_all = [SimData.p_AC; SimData.p_EV];
    c_all = [OptIn.DeviceCosts.c_AC; OptIn.DeviceCosts.c_EV];

    % --- [保留] 计算基准互补性指标 (Ori) ---
    fprintf('正在计算 %s 基准 (Ori) 互补性指标...\n', this_dir);
    % 使用绝对值判断基准参与状态 (兼容下调的负值)
    u_baseline = [abs(SimData.p_AC) > 1e-4; abs(SimData.p_EV) > 1e-4];
    % Metrics_ori = calculate_solution_metrics(u_baseline, SimData, 'ori');
    % u_baseline = [abs(SimData.p_AC) ~= 1e-4; abs(SimData.p_EV) ~= 1e-4];
    Metrics_ori = calculate_solution_metrics(u_baseline, SimData, 'ori');

    % --- 2.4 逐时贪心优化 ---
    fprintf('正在执行 %d 步逐时优化...\n', SimData.T);
    u_opt = false(SimData.N_total, SimData.T);
    p_opt_t = zeros(1, SimData.T);
    cost_total = 0;

    parfor t = 1:SimData.T
        Network_t = OptIn.Network; 
        Network_t.BaseFlow = OptIn.Network.BaseFlow(:, t);
        
        [u_t, p_t, c_t, ~] = solve_times_step_greedy_multi_obj(...
            p_all(:,t), c_all, OptIn.P_req(t), Network_t, ...
            W_cost, W_complementarity, Ratio_AC_t(t), dev_types);
            
        u_opt(:,t) = u_t;
        p_opt_t(t) = p_t;
        cost_total = cost_total + c_t;
    end
    fprintf('优化完成。总成本: %.2f\n', cost_total);

    % --- [保留] 2.5 计算优化后互补性指标 (Opt) ---
    fprintf('正在计算 %s 优化后 (Opt) 互补性指标...\n', this_dir);
    Metrics_opt = calculate_solution_metrics(u_opt, SimData, 'opt');

    % --- [保留] 2.6 结果对比报告 ---
    fprintf('\n--- [%s] 互补性指标对比 (目标: Opt < Ori) ---\n', this_dir);
    fprintf('SDCI: Ori=%.4f -> Opt=%.4f\n', Metrics_ori.SDCI, Metrics_opt.SDCI);
    fprintf('Rho:  Ori=%.4f -> Opt=%.4f\n', Metrics_ori.Rho, Metrics_opt.Rho);

    % --- [保留] 2.7 结果可视化 ---
    if ~isdeployed
        figure('Name', ['Optimization Results - ', this_dir], 'Color', 'w', 'Position', [100,100,800,400]);
        plot(1:SimData.T, OptIn.P_req, 'k--', 'LineWidth', 2, 'DisplayName', '电网需求');
        hold on;
        plot(1:SimData.T, p_opt_t, 'r-', 'LineWidth', 1.5, 'DisplayName', '优化调度出力');
        % 绘制总潜力作为参考 (使用直接加载的聚合数据)
        plot(1:SimData.T, SimData.P_AC_agg_loaded + SimData.P_EV_agg_loaded, ...
             'Color', [0.7 0.7 0.7], 'DisplayName', '系统总潜力');
        
        legend('Location', 'bestoutside'); grid on;
        xlabel('时间步'); ylabel('功率 (kW)');
        title(sprintf('%s 方向调节优化结果 (W_{cost}=%.1f, W_{comp}=%.1f)', this_dir, W_cost, W_complementarity));
    end
end

toc;
fprintf('\n全部完成。\n');