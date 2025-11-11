
clear; close all; clc;
tic;

plot_data_filename = 'multi_obj_plot_data3.mat'; 

%% 1. 全局配置
chunk_results_dir = 'chunk_results_abs_hour'; 

% [!!! 修改点：在这里配置要加载的文件和数量 !!!]
% 您可以按需添加或修改这个结构体数组
% 使用 'inf' 表示加载该文件中的所有设备
chunk_load_config = [
    struct('file', 'results_chunk_1.mat', 'ac_rows', 0, 'ev_rows', 100),
 struct('file', 'results_chunk_1.mat', 'ac_rows', 0, 'ev_rows', 100),
  struct('file', 'results_chunk_1.mat', 'ac_rows', 0, 'ev_rows', 100),
];
% [!!! 修改结束 !!!]


run_directions = {'Up', 'Down'}; 
% run_directions = { 'Down'}; 

W_cost = 0.6;              % 成本权重
W_complementarity = 0.4;   % 互补性权重 (SDCI & Rho)

plotData = struct(); 

fprintf('启动 VPP 多目标优化 (Cost=%.1f, Comp=%.1f)...\n', W_cost, W_complementarity);

%% 2. 主循环：按方向依次执行
for d_idx = 1:length(run_directions)
    this_dir = run_directions{d_idx};
    fprintf('\n>>> 开始处理 %s (调节) 方向 >>>\n', upper(this_dir));

    % --- 2.1 加载数据 ---
    % [!!! 修改点：传递新的配置结构体 !!!]
    SimData = load_simulation_data_from_chunks(chunk_results_dir, chunk_load_config, this_dir);
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
    u_baseline = [abs(SimData.p_AC) > 1e-4; abs(SimData.p_EV) > 1e-4];
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

    % --- [新增] 2.7 将绘图数据暂存到结构体 ---
    fprintf('正在为 %s 方向暂存绘图数据...\n', this_dir);
    plotData.(this_dir).time_axis = (1:SimData.T);
    plotData.(this_dir).P_req = OptIn.P_req; % (1 x T)
    plotData.(this_dir).p_opt_t = p_opt_t; % (1 x T)
    plotData.(this_dir).P_potential_total = SimData.P_AC_agg_loaded + SimData.P_EV_agg_loaded; % (1 x T)
    plotData.(this_dir).W_cost = W_cost;
    plotData.(this_dir).W_complementarity = W_complementarity;
    
end % 结束主循环

toc;
fprintf('\n全部完成。\n');

%% 3. [新增] 保存绘图数据到 MAT 文件
try
    fprintf('\n正在将所有绘图数据保存到 %s ...\n', plot_data_filename);
    save(plot_data_filename, 'plotData');
    fprintf('数据保存成功。\n');
catch ME_save
    fprintf('*** 保存绘图 MAT 文件失败: %s ***\n', ME_save.message);
end

% %% 4. [新增] 从 MAT 文件加载数据并执行绘图
% if ~isdeployed && exist(plot_data_filename, 'file')
%     fprintf('\n正在从 %s 加载数据并开始绘图...\n', plot_data_filename);
%     loadedPlotData = load(plot_data_filename);
% 
%     if isfield(loadedPlotData, 'plotData')
%         plotData = loadedPlotData.plotData;
%         saved_directions = fieldnames(plotData);
% 
%         for i = 1:length(saved_directions)
%             this_dir = saved_directions{i};
%             data = plotData.(this_dir);
% 
%             fprintf('  正在绘制 %s 方向的图像...\n', this_dir);
% 
%             figure('Name', ['Optimization Results - ', this_dir, ' (from MAT)'], 'Color', 'w', 'Position', [100+i*50, 100+i*50, 800, 400]);
%             plot(data.time_axis, data.P_req, 'k--', 'LineWidth', 2, 'DisplayName', '电网需求');
%             hold on;
%             plot(data.time_axis, data.p_opt_t, 'r-', 'LineWidth', 1.5, 'DisplayName', '优化调度出力');
%             plot(data.time_axis, data.P_potential_total, ...
%                  'Color', [0.7 0.7 0.7], 'DisplayName', '系统总潜力');
% 
%             legend('Location', 'bestoutside'); grid on;
%             xlabel('时间步'); ylabel('功率 (kW)');
%             title(sprintf('%s 方向调节优化结果 (W_{cost}=%.1f, W_{comp}=%.1f)', this_dir, data.W_cost, data.W_complementarity));
%         end
%         fprintf('绘图完成。\n');
%     else
%         fprintf('错误: %s 中未找到 "plotData" 结构体。\n', plot_data_filename);
%     end
% elseif isdeployed
%      fprintf('\n处于部署环境 (isdeployed=true)，跳过绘图。\n');
% else
%      fprintf('\n未找到绘图数据文件 %s，跳过绘图。\n', plot_data_filename);
% end