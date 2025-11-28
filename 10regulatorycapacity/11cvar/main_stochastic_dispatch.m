clear; close all; clc;

%% 1. 加载场景集与可靠调节域数据
data_file = 'reliable_regulation_domain.mat';
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
% 弃风/切负荷惩罚 (用于松弛变量)
cost_params.c_slack = 10000; 

% --- 风险偏好参数 (CVaR) ---
risk_params.beta = 10;      
risk_params.confidence = 0.95; 
risk_params.rho_pen = 100;   

% --- 【关键修改】生成物理可行的电网需求 ---
% 计算系统在每一时刻的“最大可靠调节能力”
Max_System_Capacity = Reliable_AC_Up + Reliable_EV_Up;

% 生成需求：设定为最大能力的 40% ~ 80% 之间，确保一定小于 100%
% 留出 10% 的裕度防止数值误差导致无解
rng(100); 
safety_margin = 0.90; % 安全系数
demand_ratio = 0.4 + 0.4 * rand(T_steps, 1); % 随机比例
% 对正弦波进行包络限制，确保 P_req(t) <= 0.9 * Capacity(t)
P_grid_demand = Max_System_Capacity .* demand_ratio * safety_margin;

fprintf('已生成电网需求，最大需求占可靠容量的 %.1f%%\n', max(P_grid_demand ./ (Max_System_Capacity+1e-5))*100);

%% 3. 构建风险规避型协同调度模型 (QP转化)
fprintf('正在构建大规模二次规划(QP)模型 (包含松弛变量)...\n');

% 调用核心构建函数 (函数已更新，见下文)
[H, f, A, b, Aeq, beq, lb, ub, mapping_info] = construct_risk_constrained_qp_robust(...
    P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
    Reliable_AC_Up, Reliable_EV_Up, ... 
    cost_params, risk_params);

%% 4. 求解优化问题
fprintf('开始求解 QP (变量维度: %d)...\n', length(f));
% 增大容差以提高求解成功率
options = optimoptions('quadprog', 'Display', 'iter-detailed', ...
    'Algorithm', 'interior-point-convex', 'ConstraintTolerance', 1e-5);

tic;
[x_opt, fval, exitflag, output] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
solve_time = toc;

%% 5. 结果解析与保护
if exitflag > 0
    fprintf('\n====== 优化求解成功 (耗时: %.4f 秒) ======\n', solve_time);
    
    % 提取决策变量
    P_AC_opt = x_opt(mapping_info.idx_P_AC);
    P_EV_opt = x_opt(mapping_info.idx_P_EV);
    P_Slack  = x_opt(mapping_info.idx_Slack); % 松弛变量
    
    % 检查是否有切负荷
    max_slack = max(abs(P_Slack));
    if max_slack > 1e-3
        warning('注意：存在功率缺额 (Slack > 0)，最大缺额 %.2f kW。这说明物理容量不足。', max_slack);
    end

    % 计算后评估指标
    [SDCI, Rho] = calculate_metrics_post(P_AC_opt, P_EV_opt);
    fprintf('总成本: %.2f\n', fval);
    fprintf('协同互补性 SDCI: %.4f\n', SDCI);
    fprintf('时间相关性 Rho:  %.4f\n', Rho);

    % --- 绘图 ---
    figure('Position', [100, 100, 1200, 600], 'Color', 'w');
    
    % 子图1: 功率堆叠
    subplot(1, 2, 1);
    area(time_points, [P_AC_opt, P_EV_opt]);
    hold on;
    plot(time_points, P_grid_demand, 'k--', 'LineWidth', 2);
    if max_slack > 1e-3
        bar(time_points, P_Slack, 'FaceColor', 'r', 'EdgeColor', 'none');
        legend('AC出力', 'EV出力', '电网需求', '功率缺额(Slack)');
    else
        legend('AC出力', 'EV出力', '电网需求');
    end
    xlabel('时间 (h)'); ylabel('功率 (kW)');
    title('协同调度结果'); grid on;

    % 子图2: 可靠性校验
    subplot(1, 2, 2);
    plot(time_points, Reliable_AC_Up, 'b--', 'LineWidth', 1.5); hold on;
    plot(time_points, P_AC_opt, 'b-', 'LineWidth', 2);
    plot(time_points, Reliable_EV_Up, 'r--', 'LineWidth', 1.5); 
    plot(time_points, P_EV_opt, 'r-', 'LineWidth', 2);
    legend('AC容量上界', 'AC调度', 'EV容量上界', 'EV调度');
    xlabel('时间 (h)'); ylabel('功率 (kW)');
    title('调度指令可行性校验'); grid on;
    
else
    % 错误处理
    error('优化求解失败！Exitflag: %d. \n可能原因：即使加了松弛变量，问题仍存在严重的数值病态或约束冲突。', exitflag);
end