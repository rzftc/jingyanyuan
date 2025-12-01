%% test_suite_comprehensive_dispatch.m (逻辑完全修正版)
% 修复了之前因物理边界锁死导致Beta失效的问题
% 重新设计了需求生成逻辑，确保触发风险权衡

clear; close all; clc;

%% ================= 1. 全局初始化 =================
fprintf('正在加载场景数据...\n');
data_file = 'reliable_regulation_domain_1000.mat';
if ~exist(data_file, 'file')
    error('数据文件缺失！请先运行 main_scenario_generation_diff.m');
end
load(data_file);

% 数据维度对齐
[T_steps, N_scenarios] = size(Scenarios_AC_Up);

% --- 关键修正 1: 计算物理最大边界 ---
% 之前错误地使用了 Reliable (95%分位数) 作为硬约束 UB，导致无法探索风险区域
% 现在使用场景集中的最大值作为物理硬约束
Physical_AC_Up = max(Scenarios_AC_Up, [], 2);
Physical_EV_Up = max(Scenarios_EV_Up, [], 2);
Reliable_AC_Up = Reliable_AC_Up(:); % 仅作参考和绘图
Reliable_EV_Up = Reliable_EV_Up(:);

% --- 基础参数设置 ---
cost_params.c1_ac = 0.5;  cost_params.c2_ac = 0.01;
cost_params.c1_ev = 0.4;  cost_params.c2_ev = 0.01;
% 切负荷成本 (设为中等偏高，允许在极端风险下切负荷，但平时尽量满足)
cost_params.c_slack = 200; 

% --- 关键修正 2: 构造"风险权衡"型需求 ---
% 需求设定在 [可靠容量, 物理最大容量] 之间
% 这样，满足需求就必须动用"不可靠"的容量部分，从而触发 CVaR 惩罚
rng(105); 
Total_Reliable = Reliable_AC_Up + Reliable_EV_Up;
Total_Physical = Physical_AC_Up + Physical_EV_Up;
% 需求 = 可靠容量 + 40% * (风险容量区间)
P_grid_demand = Total_Reliable + 0.40 * (Total_Physical - Total_Reliable) .* rand(T_steps, 1);

% 网络参数
N_bus = 5; N_line = 6;
net_params.PTDF = (rand(N_line, N_bus) - 0.5) * 0.5; 
net_params.AcDist = [0.2; 0.2; 0.2; 0.2; 0.2];
net_params.EvDist = [0; 0.4; 0; 0.6; 0];
net_params.BaseFlow = zeros(N_line, T_steps);
net_params.LineLimit = ones(N_line, 1) * 1e5; % 默认无阻塞

fprintf('初始化完成。最大物理容量: %.2f kW, 最大需求: %.2f kW\n', max(Total_Physical), max(P_grid_demand));

%% ================= 场景 B: 风险偏好灵敏度分析 =================
fprintf('\n>>> 场景 B: 风险偏好灵敏度分析 <<<\n');

beta_values = [0, 10, 100];
% 初始化结果容器
b_run_cost = nan(1, length(beta_values));
b_risk_val = nan(1, length(beta_values));
b_slack_sum = nan(1, length(beta_values)); % 记录切负荷总量
strategies = cell(1, length(beta_values));

for i = 1:length(beta_values)
    beta = beta_values(i);
    fprintf('  工况 %d (Beta=%d): ', i, beta);
    
    risk_p.beta = beta;
    risk_p.confidence = 0.95;
    % 违约惩罚单价 (应小于 Slack 成本，否则永远优先切负荷)
    risk_p.rho_pen = 50; 
    
    % [关键] 传入 Physical_AC_Up 作为 UB (R_AC 参数)
    [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_robust_network_improve(...
        P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
        Physical_AC_Up, Physical_EV_Up, ... % <--- 修正此处
        cost_params, risk_p, net_params);
    
    options = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'interior-point-convex');
    [x_opt, fval, exitflag] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
    
    if exitflag > 0
        P_AC = x_opt(info.idx_P_AC);
        P_EV = x_opt(info.idx_P_EV);
        P_Slack = x_opt(info.idx_Slack);
        eta_val = x_opt(info.idx_eta);
        z_val = x_opt(info.idx_z);
        
        % 1. 发电成本 (二次项 + 一次项)
        cost_gen = sum(cost_params.c1_ac*P_AC + cost_params.c2_ac*P_AC.^2 + ...
                       cost_params.c1_ev*P_EV + cost_params.c2_ev*P_EV.^2);
                   
        % 2. 切负荷成本
        cost_slack = sum(cost_params.c_slack * abs(P_Slack));
        
        % 3. CVaR 风险值 (物理含义：预期在最坏5%情况下的平均违约量)
        cvar_val = eta_val + (1 / (N_scenarios * (1 - risk_p.confidence))) * sum(z_val);
        
        % 记录
        b_run_cost(i) = cost_gen;      
        b_slack_sum(i) = sum(abs(P_Slack)); % 记录切负荷总量
        b_risk_val(i) = cvar_val;
        
        strategies{i}.P_AC = P_AC;
        strategies{i}.P_EV = P_EV;
        
        fprintf('发电成本: %.2f, 切负荷量: %.2f kW, CVaR风险: %.2f\n', ...
            cost_gen, b_slack_sum(i), cvar_val);
    else
        fprintf('失败 (Exitflag %d)\n', exitflag);
    end
end

% --- 绘图 B ---
if any(~isnan(b_run_cost))
    figure('Name', '场景B_风险灵敏度', 'Color', 'w', 'Position', [100, 100, 900, 400]);
    
    % 左轴：切负荷量 (反映保守程度)
    yyaxis left; 
    bar(1:3, b_slack_sum, 0.5, 'FaceColor', [0.8 0.3 0.3]); 
    ylabel('总切负荷量 (kW) [越保守越高]');
    set(gca, 'XTick', 1:3, 'XTickLabel', beta_values);
    
    % 右轴：CVaR 风险值 (反映潜在违约)
    yyaxis right; 
    plot(1:3, b_risk_val, 'b-o', 'LineWidth', 2, 'MarkerSize', 8); 
    ylabel('CVaR 潜在违约风险 (kW) [越激进越高]');
    
    xlabel('风险厌恶系数 \beta');
    % title('风险偏好对调度策略的影响: 随着Beta增加，宁愿切负荷也不冒违约风险'); % [已删除标题]
    legend('切负荷量 (安全性)', '潜在违约风险 (经济性代价)', 'Location', 'best');
    grid on;
    
    % 保存为高DPI PNG
    print(gcf, '风险偏好灵敏度分析.png', '-dpng', '-r300');
end

%% ================= 场景 C: 网络阻塞管理 =================
fprintf('\n>>> 场景 C: 网络阻塞管理测试 <<<\n');

net_params_C = net_params;
% 制造明显阻塞: 限制为平均需求的 20% (非常严格，迫使调整)
limit_val = mean(P_grid_demand) * 0.2; 
net_params_C.LineLimit(1) = limit_val; 
fprintf('  - 制造阻塞: 线路1限额 %.2f kW\n', limit_val);

risk_p.beta = 10; 
[H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_robust_network(...
    P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, Physical_AC_Up, Physical_EV_Up, ...
    cost_params, risk_p, net_params_C);

[x_opt_C, ~, exitflag_C] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);

if exitflag_C > 0
    P_AC_C = x_opt_C(info.idx_P_AC);
    P_EV_C = x_opt_C(info.idx_P_EV);
    
    Sens_AC = net_params_C.PTDF * net_params_C.AcDist;
    Sens_EV = net_params_C.PTDF * net_params_C.EvDist;
    Flow_L1 = Sens_AC(1)*P_AC_C + Sens_EV(1)*P_EV_C;
    
    figure('Name', '场景C_网络阻塞', 'Color', 'w', 'Position', [100, 550, 800, 400]);
    plot(1:T_steps, Flow_L1, 'b-', 'LineWidth', 1.5); hold on;
    yline(limit_val, 'r--', 'LineWidth', 2, 'Label', '线路限额');
    yline(-limit_val, 'r--', 'LineWidth', 2);
    ylabel('线路1 潮流 (kW)'); xlabel('时间步');
    % title('网络约束下的潮流控制效果'); % [已删除标题]
    grid on;
    
    % 保存为高DPI PNG
    print(gcf, '网络阻塞管理测试.png', '-dpng', '-r300');
    
    fprintf('  - 求解成功。最大潮流: %.2f (限额 %.2f)\n', max(abs(Flow_L1)), limit_val);
else
    fprintf('  - 求解失败。\n');
end

%% ================= 场景 D: 鲁棒性测试 =================
fprintf('\n>>> 场景 D: 极端场景鲁棒性测试 <<<\n');

if isempty(strategies{1}) || isempty(strategies{3})
    fprintf('  - 警告: 策略数据缺失，跳过。\n');
else
    % 1. 挑选“黑天鹅”场景（总容量最小的场景）
    Total_Cap_Scen = Scenarios_AC_Up + Scenarios_EV_Up;
    [~, worst_idx] = min(sum(Total_Cap_Scen)); % 全时段总能量最小
    fprintf('  - 最恶劣场景: #%d\n', worst_idx);

    Real_Cap_AC = Scenarios_AC_Up(:, worst_idx);
    Real_Cap_EV = Scenarios_EV_Up(:, worst_idx);

    % 2. 策略对比
    % 计算物理违约量：调度指令超过该特定场景实际能力的部分
    calc_viol = @(P_ac, P_ev) sum(max(0, P_ac - Real_Cap_AC) + max(0, P_ev - Real_Cap_EV));
    
    viol_neutral = calc_viol(strategies{1}.P_AC, strategies{1}.P_EV);
    viol_robust  = calc_viol(strategies{3}.P_AC, strategies{3}.P_EV);
    
    fprintf('  - 中性策略(Beta=0) 违约量: %.2f kW\n', viol_neutral);
    fprintf('  - 规避策略(Beta=100) 违约量: %.2f kW\n', viol_robust);
    
    figure('Name', '场景D_鲁棒性', 'Color', 'w', 'Position', [600, 300, 500, 400]);
    b = bar([viol_neutral, viol_robust], 0.5);
    b.FaceColor = 'flat';
    b.CData(1,:) = [0.8 0.2 0.2]; % Red
    b.CData(2,:) = [0.2 0.6 0.2]; % Green
    
    set(gca, 'XTickLabel', {'中性 (\beta=0)', '规避 (\beta=100)'});
    ylabel('极端场景实际违约量 (kW)');
    % title(['黑天鹅场景(#', num2str(worst_idx), ')下的实际执行违约对比']); % [已删除标题]
    grid on;
    
    % 保存为高DPI PNG
    print(gcf, '极端场景鲁棒性测试.png', '-dpng', '-r300');
end

fprintf('\n测试结束。\n');