%% test_suite_comprehensive_dispatch_ieee5.m (IEEE 5节点系统 + 逻辑修正版)
% 修正物理边界逻辑，确保风险权衡有效；采用 IEEE PJM 5-Bus 拓扑。
clear; close all; clc;

%% ================= 1. 全局初始化 =================
fprintf('正在加载场景数据...\n');
data_file = 'reliable_regulation_domain_1000_mix_01.mat';
if ~exist(data_file, 'file')
    error('数据文件缺失！请先运行 main_scenario_generation_diff.m');
end
load(data_file);

% 数据维度对齐
[T_steps, N_scenarios] = size(Scenarios_AC_Up);

% --- [修改] 统一时间轴设置 ---
% 确保时间轴是从 6:00 到 30:00
if exist('time_points', 'var')
    t_axis = time_points;
else
    % 如果数据中没有 time_points，手动生成 (假设dt=5min)
    dt_sim = 5/60; 
    t_axis = 6 : dt_sim : (6 + (T_steps-1)*dt_sim);
end

% 定义通用的时间轴刻度和标签
x_ticks_set = [6, 12, 18, 24, 30];
x_labels_set = {'06:00', '12:00', '18:00', '24:00', '06:00(+1)'};

% --- 修正1: 物理边界设置 ---
% 使用场景集最大值作为物理硬约束，而非 95% 分位数，以便探索风险区域
Physical_AC_Up = max(Scenarios_AC_Up, [], 2);
Physical_EV_Up = max(Scenarios_EV_Up, [], 2);
Reliable_AC_Up = Reliable_AC_Up(:);
Reliable_EV_Up = Reliable_EV_Up(:);

% --- 基础参数 ---
cost_params.c1_ac = 0.5;  cost_params.c2_ac = 0.01;
cost_params.c1_ev = 0.4;  cost_params.c2_ev = 0.01;
cost_params.c_slack = 200; % 切负荷惩罚 (单位: 元/kW)

% --- 修正2: 构造风险权衡需求 ---
% 需求设为 [可靠容量, 物理最大容量] 之间，强制触发 CVaR 权衡
rng(105); 
Total_Reliable = Reliable_AC_Up + Reliable_EV_Up;
Total_Physical = Physical_AC_Up + Physical_EV_Up;
P_grid_demand = Total_Reliable + 0.40 * (Total_Physical - Total_Reliable) .* rand(T_steps, 1);

% --- 网络参数 (IEEE/PJM 5-Bus Standard) ---
N_bus = 5; 
N_line = 6;

% 1. 线路参数: From, To, X (p.u.)
LineData = [
    1, 2, 0.0281;
    1, 4, 0.0304;
    1, 5, 0.0064;
    2, 3, 0.0108;
    3, 4, 0.0297;
    4, 5, 0.0297
];

% 2. 计算 PTDF (直流潮流)
B_bus = zeros(N_bus, N_bus);
B_line = zeros(N_line, N_bus);

for l = 1:N_line
    f = LineData(l, 1);
    t = LineData(l, 2);
    x_val = LineData(l, 3);
    b_val = 1/x_val;
    
    % 更新导纳矩阵 B_bus
    B_bus(f, f) = B_bus(f, f) + b_val;
    B_bus(t, t) = B_bus(t, t) + b_val;
    B_bus(f, t) = B_bus(f, t) - b_val;
    B_bus(t, f) = B_bus(t, f) - b_val;
    
    % 更新关联矩阵 B_line
    B_line(l, f) = b_val;
    B_line(l, t) = -b_val;
end

% 选节点1为平衡节点，计算 PTDF
B_bus_reduced = B_bus(2:end, 2:end);
PTDF_reduced = B_line(:, 2:end) / B_bus_reduced;
net_params.PTDF = [zeros(N_line, 1), PTDF_reduced];

% 3. 资源分布
net_params.AcDist = [0.1; 0.15; 0.3; 0.3; 0.15]; % AC 分布
net_params.EvDist = [0.0; 0.2; 0.4; 0.4; 0.0];  % EV 分布 (负荷中心)

% 4. 线路容量与基础潮流
net_params.BaseFlow = zeros(N_line, T_steps);
Capacity_Base = 2.5e5; 
net_params.LineLimit = [
    Capacity_Base;       % 1-2
    Capacity_Base;       % 1-4
    Capacity_Base * 2;   % 1-5
    Capacity_Base * 0.8; % 2-3
    Capacity_Base;       % 3-4
    Capacity_Base * 0.6  % 4-5 (薄弱断面)
];

fprintf('初始化完成。最大物理容量: %.2f kW, 最大需求: %.2f kW\n', max(Total_Physical), max(P_grid_demand));

%% ================= 场景 B: 风险偏好灵敏度分析 =================
fprintf('\n>>> 场景 B: 风险偏好灵敏度分析 <<<\n');

beta_values = [0, 10, 100];
b_run_cost = nan(1, length(beta_values));
b_risk_val = nan(1, length(beta_values));
b_slack_sum = nan(1, length(beta_values));
strategies = cell(1, length(beta_values));

for i = 1:length(beta_values)
    beta = beta_values(i);
    fprintf('  工况 %d (Beta=%d): ', i, beta);
    
    risk_p.beta = beta;
    risk_p.confidence = 0.95;
    risk_p.rho_pen = 50; 
    
    % 关键: 传入 Physical_AC_Up 作为 UB
    [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_robust_network_improve(...
        P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
        Physical_AC_Up, Physical_EV_Up, ... 
        cost_params, risk_p, net_params);
    
    options = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'interior-point-convex');
    [x_opt, fval, exitflag] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
    
    if exitflag > 0
        P_AC = x_opt(info.idx_P_AC);
        P_EV = x_opt(info.idx_P_EV);
        P_Slack = x_opt(info.idx_Slack);
        eta_val = x_opt(info.idx_eta);
        z_val = x_opt(info.idx_z);
        
        % 1. 发电成本
        cost_gen = sum(cost_params.c1_ac*P_AC + cost_params.c2_ac*P_AC.^2 + ...
                       cost_params.c1_ev*P_EV + cost_params.c2_ev*P_EV.^2);
                    
        % 2. 切负荷成本
        cost_slack = sum(cost_params.c_slack * abs(P_Slack));

        % 3. 总运行成本 (用于验证理论：风险规避策略虽然发电成本低，但总成本高)
        total_real_cost = cost_gen + cost_slack;
        
        % 4. CVaR 风险值
        cvar_val = eta_val + (1 / (N_scenarios * (1 - risk_p.confidence))) * sum(z_val);
        
        % 记录结果
        b_run_cost(i) = cost_gen;       
        b_slack_sum(i) = sum(abs(P_Slack));
        b_risk_val(i) = cvar_val;
        
        strategies{i}.P_AC = P_AC;
        strategies{i}.P_EV = P_EV;
        strategies{i}.P_Slack = P_Slack; 
        
        % [修改] 增加总成本显示
        fprintf('发电成本: %.2f, 切负荷量: %.2f kW, 总运行成本: %.2f, CVaR风险: %.2f\n', ...
            cost_gen, b_slack_sum(i), total_real_cost, cvar_val);
    else
        fprintf('失败 (Exitflag %d)\n', exitflag);
    end
end

% --- 绘图 B1: 风险灵敏度统计 (柱状图) ---
if any(~isnan(b_run_cost))
    figure('Name', '场景B_风险灵敏度', 'Color', 'w', 'Position', [100, 100, 900, 400]);
    
    yyaxis left; 
    bar(1:3, b_slack_sum, 0.5, 'FaceColor', [0.8 0.3 0.3]); 
    ylabel('总切负荷量 (kW) [越保守越高]');
    set(gca, 'XTick', 1:3, 'XTickLabel', beta_values);
    
    yyaxis right; 
    plot(1:3, b_risk_val, 'b-o', 'LineWidth', 2, 'MarkerSize', 8); 
    ylabel('CVaR 潜在违约风险 (kW) [越激进越高]');
    
    xlabel('风险厌恶系数 \beta');
    legend('切负荷量 (安全性)', '潜在违约风险 (经济性代价)', 'Location', 'best');
    grid on;
    
    print(gcf, '风险偏好灵敏度分析.png', '-dpng', '-r300');
end

%% ================= [新增] 绘制具体调度方案 (以 Beta=10 为例) =================
fprintf('\n>>> 正在绘制详细调度方案 (选择 Beta=10) ...\n');
idx_plot = 2; % 对应 beta_values = [0, 10, 100] 中的 10 (适中策略)

if idx_plot <= length(strategies) && ~isempty(strategies{idx_plot})
    % 提取数据
    P_AC_opt = strategies{idx_plot}.P_AC;
    P_EV_opt = strategies{idx_plot}.P_EV;
    
    % 计算实际功率缺额 (Slack)
    % 优先使用优化器返回的 Slack，或者重新计算 P_req - P_supply
    if isfield(strategies{idx_plot}, 'P_Slack')
        P_Slack_calc = strategies{idx_plot}.P_Slack;
    else
        P_Supply = P_AC_opt + P_EV_opt;
        P_Slack_calc = P_grid_demand - P_Supply;
    end
    P_Slack_calc(P_Slack_calc < 1e-5) = 0; % 修正微小数值误差
    
    % --- 绘图: 协同调度堆叠图 & 边界校验 ---
    figure('Name', '协同调度详细方案 (Beta=10)', 'Color', 'w', 'Position', [150, 150, 1200, 800]);
    
    % === 子图 1: 功率平衡堆叠图 (Stacked Area Plot) ===
    subplot(2, 1, 1);
    hold on;
    
    % 1. 准备堆叠矩阵: [AC, EV, Slack]
    Y_Stack = [P_AC_opt, P_EV_opt, P_Slack_calc];
    
    % 2. 使用 area 绘制堆叠图
    h_area = area(t_axis, Y_Stack);
    
    % 3. 调整颜色和样式
    % AC (底层): 蓝色
    h_area(1).FaceColor = [0.00, 0.45, 0.74]; 
    h_area(1).EdgeColor = 'none'; 
    
    % EV (中层): 绿色
    h_area(2).FaceColor = [0.47, 0.67, 0.19]; 
    h_area(2).EdgeColor = 'none';
    
    % Slack (顶层): 红色
    h_area(3).FaceColor = [0.85, 0.33, 0.10]; 
    h_area(3).EdgeColor = 'none';
    h_area(3).FaceAlpha = 0.8;
    
    % 4. 绘制总需求曲线 (黑色虚线)
    plot(t_axis, P_grid_demand, 'k--', 'LineWidth', 2.0, 'DisplayName', '电网总需求 (Target)');
    
    ylabel('功率 (kW)');
    title(sprintf('源荷协同调度功率堆叠图 (风险厌恶系数 \\beta=%d)', beta_values(idx_plot)), 'FontSize', 14);
    
    legend([h_area(1), h_area(2), h_area(3)], ...
           {'空调出力 (AC)', '电动汽车出力 (EV)', '功率缺额 (Slack)'}, ...
           'Location', 'northwest', 'FontSize', 10);
           
    grid on;
    xlim([6, 30]); % 锁定时间轴范围
    % 设置X轴刻度
    set(gca, 'XTick', x_ticks_set, 'XTickLabel', x_labels_set);
    ylim([0, max(P_grid_demand) * 1.15]);
    
    % === 子图 2: 调度指令与可靠边界对比 (验证鲁棒性) ===
    subplot(2, 1, 2);
    hold on;
    
    % AC 对比 (蓝色系)
    plot(t_axis, Reliable_AC_Up, '--', 'Color', [0.6, 0.6, 1.0], 'LineWidth', 1.5, 'DisplayName', 'AC 可靠上界 (95%)');
    plot(t_axis, P_AC_opt, '-', 'Color', [0.00, 0.45, 0.74], 'LineWidth', 2.0, 'DisplayName', 'AC 实际调度');
    
    % EV 对比 (绿色系)
    plot(t_axis, Reliable_EV_Up, '--', 'Color', [0.6, 0.9, 0.6], 'LineWidth', 1.5, 'DisplayName', 'EV 可靠上界 (95%)');
    plot(t_axis, P_EV_opt, '-', 'Color', [0.47, 0.67, 0.19], 'LineWidth', 2.0, 'DisplayName', 'EV 实际调度');
    
    ylabel('调节功率 (kW)'); 
    xlabel('时间');
    title('调度指令 vs 可靠调节域边界 (鲁棒性校验)', 'FontSize', 14);
    
    legend('Location', 'best', 'NumColumns', 2, 'FontSize', 10);
    grid on;
    xlim([6, 30]); % 锁定时间轴范围
    set(gca, 'XTick', x_ticks_set, 'XTickLabel', x_labels_set);
    
    % 保存图像
    print(gcf, '协同调度详细方案_Beta10_Optimized.png', '-dpng', '-r300');
    fprintf('  >>> 优化后的调度方案图已保存: 协同调度详细方案_Beta10_Optimized.png\n');
else
    warning('未找到 Beta=10 的策略数据，跳过详细绘图。');
end

%% ================= 场景 C: 网络阻塞管理 =================
fprintf('\n>>> 场景 C: 网络阻塞管理测试 <<<\n');

net_params_C = net_params;
% 制造阻塞: 限制为平均需求的 20%
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
    
    % [修改] 使用 t_axis (6-30) 绘制潮流曲线
    plot(t_axis, Flow_L1, 'b-', 'LineWidth', 1.5); hold on;
    yline(limit_val, 'r--', 'LineWidth', 2, 'Label', '线路限额');
    yline(-limit_val, 'r--', 'LineWidth', 2);
    
    ylabel('线路1 潮流 (kW)'); 
    xlabel('时间');
    grid on;
    xlim([6, 30]);
    set(gca, 'XTick', x_ticks_set, 'XTickLabel', x_labels_set);
    
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
    % 1. 挑选最恶劣场景 (总能量最小)
    Total_Cap_Scen = Scenarios_AC_Up + Scenarios_EV_Up;
    [~, worst_idx] = min(sum(Total_Cap_Scen));
    fprintf('  - 最恶劣场景: #%d\n', worst_idx);

    Real_Cap_AC = Scenarios_AC_Up(:, worst_idx);
    Real_Cap_EV = Scenarios_EV_Up(:, worst_idx);

    % 2. 策略对比 (计算物理违约量)
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
    grid on;
    
    print(gcf, '极端场景鲁棒性测试.png', '-dpng', '-r300');
end

fprintf('\n测试结束。\n');