%% test_suite_comprehensive_dispatch.m (修正版)
% 第四章算例综合测试套件：风险灵敏度(B)、网络阻塞(C)、鲁棒性(D)

clear; close all; clc;

%% ================= 1. 全局初始化与数据加载 =================
fprintf('正在加载场景数据...\n');
data_file = 'reliable_regulation_domain_1000.mat';
if ~exist(data_file, 'file')
    error('数据文件缺失！请先运行 main_scenario_generation_diff.m');
end
load(data_file);

% 数据维度对齐
Reliable_AC_Up = Reliable_AC_Up(:);
Reliable_EV_Up = Reliable_EV_Up(:);
[T_steps, N_scenarios] = size(Scenarios_AC_Up);

% --- 基础参数设置 ---
cost_params.c1_ac = 0.5;  cost_params.c2_ac = 0.01;
cost_params.c1_ev = 0.4;  cost_params.c2_ev = 0.01;
cost_params.c_slack = 5000; % 切负荷惩罚

% [关键修改] 生成更高水平的电网需求，制造“困难”场景
% 设定为最大容量的 85% - 115%，确保部分时刻必须动用风险储备或切负荷
Max_Cap = Reliable_AC_Up + Reliable_EV_Up;
rng(105); % 使用新的随机种子
P_grid_demand = Max_Cap .* (0.85 + 0.30 * rand(T_steps, 1));

% 初始化网络参数
N_bus = 5; N_line = 6;
net_params.PTDF = (rand(N_line, N_bus) - 0.5) * 0.5; 
net_params.AcDist = [0.2; 0.2; 0.2; 0.2; 0.2];
net_params.EvDist = [0; 0.4; 0; 0.6; 0];
net_params.BaseFlow = zeros(N_line, T_steps);
net_params.LineLimit = ones(N_line, 1) * 1e5; % 默认无阻塞

fprintf('初始化完成。开始执行测试序列...\n');

%% ================= 测试场景 B: 风险偏好灵敏度分析 =================
fprintf('\n>>> 正在执行场景 B: 风险偏好灵敏度分析 (Risk Sensitivity) <<<\n');

beta_values = [0, 10, 100];
% 预分配结果数组，避免绘图长度错误
b_run_cost = nan(1, length(beta_values));
b_risk_cost = nan(1, length(beta_values));
strategies = cell(1, length(beta_values)); % 存储策略供场景D使用

for i = 1:length(beta_values)
    beta = beta_values(i);
    fprintf('  - 测试工况 %d: Beta = %d ... ', i, beta);
    
    risk_p.beta = beta;
    risk_p.confidence = 0.95;
    risk_p.rho_pen = 100;
    
    % 构建并求解
    [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_robust_network(...
        P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, Reliable_AC_Up, Reliable_EV_Up, ...
        cost_params, risk_p, net_params);
    
    options = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'interior-point-convex');
    [x_opt, fval, exitflag] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
    
    if exitflag > 0
        % 提取结果
        P_AC = x_opt(info.idx_P_AC);
        P_EV = x_opt(info.idx_P_EV);
        P_Slack = x_opt(info.idx_Slack);
        
        % 计算经济成本 (不含风险项)
        cost_running = sum(cost_params.c1_ac*P_AC + cost_params.c2_ac*P_AC.^2 + ...
                           cost_params.c1_ev*P_EV + cost_params.c2_ev*P_EV.^2 + ...
                           cost_params.c_slack * abs(P_Slack)); % Slack 算作运行惩罚
        
        % 计算风险成本 (总目标 - 运行成本)
        cost_risk = fval - cost_running; 
        
        % 存储数据
        b_run_cost(i) = cost_running;
        b_risk_cost(i) = cost_risk;
        strategies{i}.P_AC = P_AC;
        strategies{i}.P_EV = P_EV;
        
        fprintf('成功。运行成本: %.2f, 风险部分: %.2f\n', cost_running, cost_risk);
    else
        fprintf('失败 (Exitflag %d)\n', exitflag);
    end
end

% --- 绘图 B ---
if any(~isnan(b_run_cost))
    figure('Name', '场景B_风险灵敏度分析', 'Color', 'w', 'Position', [100, 100, 800, 400]);
    yyaxis left; 
    bar(1:3, b_run_cost, 0.5, 'FaceColor', [0.2 0.6 0.8]); 
    ylabel('运行经济成本 (元)');
    set(gca, 'XTick', 1:3, 'XTickLabel', beta_values);
    
    yyaxis right; 
    plot(1:3, b_risk_cost, 'r-o', 'LineWidth', 2); 
    ylabel('CVaR风险项数值');
    
    xlabel('风险厌恶系数 \beta');
    title('风险偏好对调度策略的影响 (Pareto权衡)');
    grid on;
end

%% ================= 测试场景 C: 网络阻塞管理 =================
fprintf('\n>>> 正在执行场景 C: 网络阻塞管理测试 (Network Congestion) <<<\n');

net_params_C = net_params;
% 制造更明显的阻塞：取平均需求的 40%
limit_val = mean(P_grid_demand) * 0.4; 
net_params_C.LineLimit(1) = limit_val; 
fprintf('  - 设置线路1容量限制为: %.2f kW\n', limit_val);

risk_p.beta = 10; 
[H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_robust_network_improve(...
    P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, Reliable_AC_Up, Reliable_EV_Up, ...
    cost_params, risk_p, net_params_C);

[x_opt_C, ~, exitflag_C] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);

if exitflag_C > 0
    P_AC_C = x_opt_C(info.idx_P_AC);
    P_EV_C = x_opt_C(info.idx_P_EV);
    
    % 计算实际潮流
    Sens_AC = net_params_C.PTDF * net_params_C.AcDist;
    Sens_EV = net_params_C.PTDF * net_params_C.EvDist;
    Flow_L1 = Sens_AC(1)*P_AC_C + Sens_EV(1)*P_EV_C;
    
    figure('Name', '场景C_网络阻塞管理', 'Color', 'w', 'Position', [100, 550, 800, 400]);
    plot(1:T_steps, Flow_L1, 'b-', 'LineWidth', 1.5); hold on;
    yline(limit_val, 'r--', 'LineWidth', 2, 'Label', '线路上限');
    yline(-limit_val, 'r--', 'LineWidth', 2);
    ylabel('线路1 有功潮流 (kW)'); xlabel('时间步');
    title('网络约束下的潮流控制效果');
    grid on;
    fprintf('  - 优化求解成功。最大潮流: %.2f (限额 %.2f)\n', max(abs(Flow_L1)), limit_val);
else
    fprintf('  - 求解失败。\n');
end

%% ================= 测试场景 D: 鲁棒性测试 =================
fprintf('\n>>> 正在执行场景 D: 极端场景鲁棒性测试 (Robustness) <<<\n');

if isempty(strategies{1}) || isempty(strategies{3})
    fprintf('  - 警告：场景B中 Beta=0 或 Beta=100 未成功求解，跳过鲁棒性测试。\n');
else
    % 1. 挑选“黑天鹅”场景（总容量最小的场景）
    Total_Cap = sum(Scenarios_AC_Up + Scenarios_EV_Up, 1);
    [~, worst_idx] = min(Total_Cap);
    fprintf('  - 最恶劣场景索引: #%d\n', worst_idx);

    Real_Cap_AC = Scenarios_AC_Up(:, worst_idx);
    Real_Cap_EV = Scenarios_EV_Up(:, worst_idx);

    % 2. 策略对比
    Strategy_Neutral = strategies{1}; % Beta=0
    Strategy_Robust  = strategies{3}; % Beta=100

    % 计算物理违约量（调度指令超过该场景实际能力的部分）
    calc_viol = @(P_ac, P_ev) sum(max(0, P_ac - Real_Cap_AC) + max(0, P_ev - Real_Cap_EV));
    
    viol_neutral = calc_viol(Strategy_Neutral.P_AC, Strategy_Neutral.P_EV);
    viol_robust  = calc_viol(Strategy_Robust.P_AC,  Strategy_Robust.P_EV);
    
    fprintf('  - 风险中性策略(Beta=0) 违约量: %.2f kW\n', viol_neutral);
    fprintf('  - 风险规避策略(Beta=100) 违约量: %.2f kW\n', viol_robust);
    
    figure('Name', '场景D_鲁棒性对比', 'Color', 'w', 'Position', [600, 300, 500, 400]);
    bar([viol_neutral, viol_robust], 0.4, 'FaceColor', [0.8 0.3 0.3]);
    set(gca, 'XTickLabel', {'Beta=0 (中性)', 'Beta=100 (规避)'});
    ylabel('极端场景违约电量 (kWh)');
    title('黑天鹅场景下的策略鲁棒性对比');
    grid on;
end

fprintf('\n所有测试场景执行完毕。\n');