% main_stochastic_dispatch_full.m
% 包含网络潮流约束(PTDF)的 VPP 协同调度主程序

clear; close all; clc;

%% 1. 加载场景集与可靠调节域数据
data_file = 'reliable_regulation_domain_1000.mat';
if ~exist(data_file, 'file')
    error('未找到数据文件 %s，请先运行 main_scenario_generation_diff.m 生成场景数据。', data_file);
end
load(data_file); 

% 确保维度正确 (列向量)
Reliable_AC_Up = Reliable_AC_Up(:);
Reliable_EV_Up = Reliable_EV_Up(:);
[T_steps, N_scenarios] = size(Scenarios_AC_Up);
dt = time_points(2) - time_points(1); 
fprintf('成功加载数据：T=%d, 场景数=%d, dt=%.2f小时\n', T_steps, N_scenarios, dt);

%% 2. 调度参数设置
% --- 成本系数 ---
cost_params.c1_ac = 0.5;  cost_params.c2_ac = 0.01;
cost_params.c1_ev = 0.4;  cost_params.c2_ev = 0.01;
cost_params.c_slack = 10000; % 弃风/切负荷惩罚

% --- 风险偏好参数 (CVaR) ---
risk_params.beta = 10;      
risk_params.confidence = 0.95; 
risk_params.rho_pen = 100;   

% --- 生成物理可行的电网需求 ---
Max_System_Capacity = Reliable_AC_Up + Reliable_EV_Up;
rng(100); 
safety_margin = 0.90; 
demand_ratio = 0.4 + 0.4 * rand(T_steps, 1); 
P_grid_demand = Max_System_Capacity .* demand_ratio * safety_margin;

fprintf('已生成电网需求，最大需求占可靠容量的 %.1f%%\n', max(P_grid_demand ./ (Max_System_Capacity+1e-5))*100);

%% 3. [新增] 网络拓扑参数初始化
fprintf('正在初始化网络参数 (PTDF)... \n');

% 假设一个简单的 5节点 6线路 系统
N_bus = 5;
N_line = 6;

% 3.1 生成 PTDF 矩阵 (示例随机生成，实际应来自 matpower 或 grid 模型)
% 行: 线路, 列: 节点
rng(42); % 固定种子保证复现
net_params.PTDF = (rand(N_line, N_bus) - 0.5) * 0.5; 

% 3.2 线路容量限制 (kW)
% 设定一个稍微宽裕的限制，避免一开始就无解，例如总需求的 60% 分流
avg_total_power = mean(P_grid_demand);
net_params.LineLimit = ones(N_line, 1) * (avg_total_power * 0.6); 

% 3.3 线路基础潮流 (Base Flow)
% 假设线路上已有一定的背景潮流 (正弦波动)
base_flow_profile = 0.2 * net_params.LineLimit(1) * sin(linspace(0, 2*pi, T_steps))';
net_params.BaseFlow = repmat(base_flow_profile', N_line, 1); % [N_line x T]

% 3.4 资源节点分布 (关键步骤)
% 因为优化变量是聚合的 P_AC 和 P_EV，我们需要定义它们在各个节点的分布比例
% 假设 AC 均匀分布，EV 集中在节点 2 和 4
net_params.AcDist = ones(N_bus, 1) / N_bus; % [0.2, 0.2, 0.2, 0.2, 0.2]
net_params.EvDist = zeros(N_bus, 1);
net_params.EvDist(2) = 0.4;
net_params.EvDist(4) = 0.6; % EV 只分布在节点 2,4

% 检查分布和为 1
assert(abs(sum(net_params.AcDist) - 1) < 1e-5, 'AC 分布比例和不为1');
assert(abs(sum(net_params.EvDist) - 1) < 1e-5, 'EV 分布比例和不为1');

fprintf('网络参数设置完成: %d 节点, %d 线路\n', N_bus, N_line);

%% 4. 构建风险规避型协同调度模型 (含网络约束)
fprintf('正在构建包含网络约束的大规模 QP 模型...\n');

% 调用新的构建函数
[H, f, A, b, Aeq, beq, lb, ub, mapping_info] = construct_risk_constrained_qp_robust_network(...
    P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
    Reliable_AC_Up, Reliable_EV_Up, ... 
    cost_params, risk_params, net_params); % 传入 net_params

%% 5. 求解优化问题
fprintf('开始求解 QP (变量维度: %d)...\n', length(f));
options = optimoptions('quadprog', 'Display', 'iter-detailed', ...
    'Algorithm', 'interior-point-convex', 'ConstraintTolerance', 1e-5);

tic;
[x_opt, fval, exitflag, output] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
solve_time = toc;

%% 6. 结果解析与绘图
if exitflag > 0
    fprintf('\n====== 优化求解成功 (耗时: %.4f 秒) ======\n', solve_time);
    
    P_AC_opt = x_opt(mapping_info.idx_P_AC);
    P_EV_opt = x_opt(mapping_info.idx_P_EV);
    P_Slack  = x_opt(mapping_info.idx_Slack);
    
    % 检查潮流是否越限
    fprintf('检查线路潮流约束...\n');
    Sens_AC = net_params.PTDF * net_params.AcDist;
    Sens_EV = net_params.PTDF * net_params.EvDist;
    
    max_line_loading = 0;
    
    % 计算实际潮流
    Real_Flows = zeros(N_line, T_steps);
    for t = 1:T_steps
        % Flow = BaseFlow + Sens * P_opt
        flow_t = net_params.BaseFlow(:, t) + Sens_AC * P_AC_opt(t) + Sens_EV * P_EV_opt(t);
        Real_Flows(:, t) = flow_t;
        
        loading = max(abs(flow_t) ./ net_params.LineLimit);
        max_line_loading = max(max_line_loading, loading);
    end
    fprintf('最大线路负载率: %.2f%%\n', max_line_loading * 100);

    % --- 绘图 1: 调度结果 ---
    figure('Position', [100, 100, 1200, 800], 'Color', 'w');
    subplot(2, 2, 1);
    area(time_points, [P_AC_opt, P_EV_opt]);
    hold on;
    plot(time_points, P_grid_demand, 'k--', 'LineWidth', 2);
    if max(abs(P_Slack)) > 1e-3
        bar(time_points, P_Slack, 'FaceColor', 'r', 'EdgeColor', 'none');
        legend('AC出力', 'EV出力', '电网需求', 'Slack');
    else
        legend('AC出力', 'EV出力', '电网需求');
    end
    xlabel('时间 (h)'); ylabel('功率 (kW)');
    title('协同调度结果'); grid on;

    subplot(2, 2, 2);
    plot(time_points, Reliable_AC_Up, 'b--', 'LineWidth', 1.5); hold on;
    plot(time_points, P_AC_opt, 'b-', 'LineWidth', 2);
    plot(time_points, Reliable_EV_Up, 'r--', 'LineWidth', 1.5); 
    plot(time_points, P_EV_opt, 'r-', 'LineWidth', 2);
    legend('AC容量', 'AC调度', 'EV容量', 'EV调度');
    title('容量约束校验'); grid on;
    
    % --- 绘图 2: 线路潮流热力图 ---
    subplot(2, 1, 2);
    % 计算负载率矩阵
    Load_Percent = abs(Real_Flows) ./ repmat(net_params.LineLimit, 1, T_steps) * 100;
    imagesc(time_points, 1:N_line, Load_Percent);
    colorbar;
    clim([0, 100]); % 固定色标 0-100%
    xlabel('时间 (h)'); ylabel('线路编号');
    title('线路负载率 (%) (PTDF约束效果)');
    
else
    error('优化求解失败！Exitflag: %d. 可能由于网络约束过紧导致不可行，请尝试增大 LineLimit。', exitflag);
end