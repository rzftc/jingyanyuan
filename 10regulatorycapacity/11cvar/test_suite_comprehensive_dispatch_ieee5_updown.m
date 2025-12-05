%% test_suite_comprehensive_dispatch_ieee5_optimized.m (随机需求 + 恢复Area绘图版)
% 修正说明：
% 1. 物理边界：同时提取上调(Up)和下调(Down)边界。
% 2. 需求生成：使用随机信号(randn)模拟电网24小时内随机的上下调节需求。
% 3. 绘图修复：恢复使用 area 堆叠图展示调度结果（根据方向自动上下翻转），保持原有风格。

clear; close all; clc;

%% ================= 1. 全局初始化 =================
fprintf('正在加载场景数据...\n');
data_file = 'reliable_regulation_domain_1000_mix_01_new.mat';
if ~exist(data_file, 'file')
    error('数据文件缺失！请先运行 main_scenario_generation_diff.m');
end
load(data_file);

% 数据维度对齐
[T_steps, N_scenarios] = size(Scenarios_AC_Up);

% --- 时间轴设置 (统一 06:00 - 30:00) ---
if exist('time_points', 'var')
    t_axis = time_points;
else
    dt_sim = 5/60; 
    t_axis = 6 : dt_sim : (6 + (T_steps-1)*dt_sim);
end
x_ticks_set = [6, 12, 18, 24, 30];
x_labels_set = {'06:00', '12:00', '18:00', '24:00', '06:00(+1)'};

% --- 修正1: 物理边界设置 (提取双向边界) ---
% [上调 Up]
Physical_AC_Up = max(Scenarios_AC_Up, [], 2);
Physical_EV_Up = max(Scenarios_EV_Up, [], 2);
Reliable_AC_Up = Reliable_AC_Up(:);
Reliable_EV_Up = Reliable_EV_Up(:);

% [下调 Down] (取绝对值计算幅度)
Physical_AC_Down = max(abs(Scenarios_AC_Down), [], 2);
Physical_EV_Down = max(abs(Scenarios_EV_Down), [], 2);
Reliable_AC_Down = abs(Reliable_AC_Down(:));
Reliable_EV_Down = abs(Reliable_EV_Down(:));

% --- 基础成本参数 ---
cost_params.c1_ac = 0.5;  cost_params.c2_ac = 0.01;
cost_params.c1_ev = 0.4;  cost_params.c2_ev = 0.01;
cost_params.c_slack = 200; % 切负荷惩罚

% --- 互补性与相关性优化权重 ---
lambda_SDCI = 10;   % 互补性惩罚权重
lambda_Rho  = 10;   % 相关性惩罚权重
Max_Iter    = 3;    % 迭代次数

% --- 修正2: 构造混合方向的风险权衡需求 (随机波动) ---
rng(105); 

% 初始化求解器所需的有效输入向量
P_grid_demand = zeros(T_steps, 1);
Effective_Scen_AC = zeros(T_steps, N_scenarios);
Effective_Scen_EV = zeros(T_steps, N_scenarios);
Effective_Phys_AC = zeros(T_steps, 1);
Effective_Phys_EV = zeros(T_steps, 1);
Effective_Reliable_AC = zeros(T_steps, 1); 
Effective_Reliable_EV = zeros(T_steps, 1); 

% [修改] 定义方向信号：使用随机噪声模拟电网需求的随机波动
% randn 生成服从标准正态分布的随机数
% >0 为上调需求 (VPP出力), <0 为下调需求 (VPP吸纳)
direction_signal = randn(T_steps, 1); 

fprintf('生成混合需求: 基于随机信号生成 (上调/下调随机切换)。\n');

for t = 1:T_steps
    is_up_regulation = direction_signal(t) >= 0;
    
    if is_up_regulation
        % === 上调时段 (Up Regulation) ===
        cap_rel = Reliable_AC_Up(t) + Reliable_EV_Up(t);
        cap_phy = Physical_AC_Up(t) + Physical_EV_Up(t);
        
        % 映射上调数据
        Effective_Scen_AC(t,:) = Scenarios_AC_Up(t,:);
        Effective_Scen_EV(t,:) = Scenarios_EV_Up(t,:);
        Effective_Phys_AC(t) = Physical_AC_Up(t);
        Effective_Phys_EV(t) = Physical_EV_Up(t);
        Effective_Reliable_AC(t) = Reliable_AC_Up(t);
        Effective_Reliable_EV(t) = Reliable_EV_Up(t);
    else
        % === 下调时段 (Down Regulation) ===
        cap_rel = Reliable_AC_Down(t) + Reliable_EV_Down(t);
        cap_phy = Physical_AC_Down(t) + Physical_EV_Down(t);
        
        % 映射下调数据 (使用绝对值幅值)
        Effective_Scen_AC(t,:) = abs(Scenarios_AC_Down(t,:));
        Effective_Scen_EV(t,:) = abs(Scenarios_EV_Down(t,:));
        Effective_Phys_AC(t) = Physical_AC_Down(t);
        Effective_Phys_EV(t) = Physical_EV_Down(t);
        Effective_Reliable_AC(t) = Reliable_AC_Down(t);
        Effective_Reliable_EV(t) = Reliable_EV_Down(t);
    end
    
    % 生成需求幅值：设定在 [可靠容量, 物理极限] 之间，触发风险
    % 需求 = 可靠容量 + 40% * (风险容量区间)
    req_magnitude = cap_rel + 0.40 * (cap_phy - cap_rel) * rand();
    P_grid_demand(t) = req_magnitude;
end

% 将“拼接”好的有效数据赋值给原变量名，以便传入 construct 函数
Scenarios_AC_Up = Effective_Scen_AC;
Scenarios_EV_Up = Effective_Scen_EV;
Physical_AC_Up  = Effective_Phys_AC;
Physical_EV_Up  = Effective_Phys_EV;
Reliable_AC_Up  = Effective_Reliable_AC;
Reliable_EV_Up  = Effective_Reliable_EV;


% --- 网络参数 (IEEE 5节点) ---
N_bus = 5; N_line = 6;
LineData = [
    1, 2, 0.0281; 1, 4, 0.0304; 1, 5, 0.0064;
    2, 3, 0.0108; 3, 4, 0.0297; 4, 5, 0.0297
];

% 计算 PTDF
B_bus = zeros(N_bus, N_bus); B_line = zeros(N_line, N_bus);
for l = 1:N_line
    f = LineData(l, 1); t = LineData(l, 2); x_val = LineData(l, 3); b_val = 1/x_val;
    B_bus(f, f) = B_bus(f, f) + b_val; B_bus(t, t) = B_bus(t, t) + b_val;
    B_bus(f, t) = B_bus(f, t) - b_val; B_bus(t, f) = B_bus(t, f) - b_val;
    B_line(l, f) = b_val; B_line(l, t) = -b_val;
end
B_bus_reduced = B_bus(2:end, 2:end);
PTDF_reduced = B_line(:, 2:end) / B_bus_reduced;
net_params.PTDF = [zeros(N_line, 1), PTDF_reduced];
net_params.AcDist = [0.1; 0.15; 0.3; 0.3; 0.15];
net_params.EvDist = [0.0; 0.2; 0.4; 0.4; 0.0];
net_params.BaseFlow = zeros(N_line, T_steps);
net_params.LineLimit = ones(N_line, 1) * 2.5e5; 

fprintf('初始化完成。\n');

%% ================= 场景 B: 风险偏好灵敏度分析 (含互补性优化) =================
fprintf('\n>>> 场景 B: 风险偏好灵敏度分析 (含 SDCI/Rho 优化) <<<\n');

beta_values = [0, 10, 100];
b_run_cost = nan(1, length(beta_values));
b_risk_val = nan(1, length(beta_values));
b_slack_sum = nan(1, length(beta_values));
strategies = cell(1, length(beta_values));

options = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'interior-point-convex');

for i = 1:length(beta_values)
    beta = beta_values(i);
    fprintf('  工况 %d (Beta=%d): \n', i, beta);
    
    risk_p.beta = beta;
    risk_p.confidence = 0.95;
    risk_p.rho_pen = 50; 
    
    % --- 迭代优化循环 ---
    P_AC_prev = zeros(T_steps, 1);
    P_EV_prev = zeros(T_steps, 1);
    
    for iter = 1:Max_Iter
        % 1. 构建基础 QP 模型
        [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_robust_network_improve(...
            P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
            Physical_AC_Up, Physical_EV_Up, ... 
            cost_params, risk_p, net_params);
        
        % 2. 动态添加互补性与相关性惩罚项
        if iter > 1
            penalty_vec_AC = lambda_SDCI * P_EV_prev;
            penalty_vec_EV = lambda_SDCI * P_AC_prev;
            
            P_EV_centered = P_EV_prev - mean(P_EV_prev);
            P_AC_centered = P_AC_prev - mean(P_AC_prev);
            penalty_rho_AC = lambda_Rho * P_EV_centered;
            penalty_rho_EV = lambda_Rho * P_AC_centered;
            
            f(info.idx_P_AC) = f(info.idx_P_AC) + penalty_vec_AC + penalty_rho_AC;
            f(info.idx_P_EV) = f(info.idx_P_EV) + penalty_vec_EV + penalty_rho_EV;
        end
        
        % 3. 求解 QP
        [x_opt, fval, exitflag] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
        
        if exitflag > 0
            P_AC_curr = x_opt(info.idx_P_AC);
            P_EV_curr = x_opt(info.idx_P_EV);
            P_Slack = x_opt(info.idx_Slack);
            eta_val = x_opt(info.idx_eta);
            z_val = x_opt(info.idx_z);
            
            P_AC_prev = P_AC_curr;
            P_EV_prev = P_EV_curr;

            cost_gen = sum(cost_params.c1_ac*P_AC_curr + cost_params.c2_ac*P_AC_curr.^2 + ...
                           cost_params.c1_ev*P_EV_curr + cost_params.c2_ev*P_EV_curr.^2);
            cost_slack = sum(cost_params.c_slack * abs(P_Slack));
            total_real_cost = cost_gen + cost_slack;
            
            cvar_val = eta_val + (1 / (N_scenarios * (1 - risk_p.confidence))) * sum(z_val);
            
            n_dummy = ones(T_steps, 1);
            val_SDCI = calculate_SDCI_local(n_dummy, n_dummy, P_AC_curr, P_EV_curr);
            val_Rho  = calculate_Rho_local(n_dummy, P_AC_curr, n_dummy, P_EV_curr);

            b_run_cost(i) = cost_gen;      
            b_slack_sum(i) = sum(abs(P_Slack));
            b_risk_val(i) = cvar_val;
            
            strategies{i}.P_AC = P_AC_curr;
            strategies{i}.P_EV = P_EV_curr;
            strategies{i}.P_Slack = P_Slack; 
            
            if iter == Max_Iter
                % 【恢复】保持原始的输出格式
                fprintf('发电成本: %.2f, 切负荷量: %.2f, 总运行成本: %.2f, CVaR风险: %.2f, rho: %.4f, sdci: %.4f\n', ...
                    cost_gen, b_slack_sum(i), total_real_cost, cvar_val, val_Rho, val_SDCI);
            end
        else
            fprintf('    Iter %d: 求解失败 (Exitflag %d)\n', iter, exitflag);
            break; 
        end
    end
end

% --- 绘图 B: 风险偏好灵敏度 ---
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

%% ================= [绘制详细调度方案 (Beta=10)] =================
fprintf('\n>>> 正在绘制详细调度方案 (选择 Beta=10, 混合上/下调) ...\n');
idx_plot = 2; 

if idx_plot <= length(strategies) && ~isempty(strategies{idx_plot})
    P_AC_opt = strategies{idx_plot}.P_AC;
    P_EV_opt = strategies{idx_plot}.P_EV;
    if isfield(strategies{idx_plot}, 'P_Slack')
        P_Slack_calc = strategies{idx_plot}.P_Slack;
    else
        P_Slack_calc = P_grid_demand - (P_AC_opt + P_EV_opt);
    end
    P_Slack_calc(P_Slack_calc < 1e-5) = 0; 
    
    % --- 绘图处理：将下调时段的数据转为负值显示 ---
    direction_sign = sign(direction_signal);
    direction_sign(direction_sign == 0) = 1; 
    
    % 【关键】恢复使用 area 绘图
    % 准备带符号的堆叠数据 [AC, EV, Slack]
    % 因为 area 在同一时刻要求所有分量同号，所以必须乘以 direction_sign
    Y_Stack_Plot = [P_AC_opt, P_EV_opt, P_Slack_calc] .* direction_sign;
    
    % 准备带符号的边界线
    Reliable_AC_plot = Reliable_AC_Up .* direction_sign;
    Reliable_EV_plot = Reliable_EV_Up .* direction_sign;
    P_Demand_plot = P_grid_demand .* direction_sign;
    
    % --- 图 1: 功率堆叠 (独立窗口) ---
    fig_stack = figure('Name', '协同调度功率堆叠', 'Color', 'w', 'Position', [150, 150, 1000, 600]);
    hold on;
    
    % 使用 area 绘制 (自动处理正负)
    h_area = area(t_axis, Y_Stack_Plot);
    
    % 恢复原来的颜色设置
    h_area(1).FaceColor = [0.00, 0.45, 0.74]; h_area(1).EdgeColor = 'none'; % AC (蓝色)
    h_area(2).FaceColor = [0.47, 0.67, 0.19]; h_area(2).EdgeColor = 'none'; % EV (绿色)
    h_area(3).FaceColor = [0.85, 0.33, 0.10]; h_area(3).EdgeColor = 'none'; h_area(3).FaceAlpha = 0.8; % Slack (橙红)
    
    % 绘制需求曲线
    plot(t_axis, P_Demand_plot, 'k--', 'LineWidth', 2.0, 'DisplayName', '电网总需求');
    
    yline(0, 'k-');
    ylabel('功率 (kW) [正=上调, 负=下调]', 'FontSize', 14, 'FontName', 'Microsoft YaHei');
    legend([h_area(1), h_area(2), h_area(3)], {'空调 (AC)', '电动汽车 (EV)', '缺额 (Slack)'}, ...
           'Location', 'northwest', 'FontSize', 12, 'FontName', 'Microsoft YaHei');
    grid on; 
    xlim([6, 30]); 
    set(gca, 'XTick', x_ticks_set, 'XTickLabel', x_labels_set, ...
        'FontSize', 12, 'FontName', 'Microsoft YaHei', 'LineWidth', 1.2);
    
    print(fig_stack, '协同调度功率堆叠_Optimized.png', '-dpng', '-r600');
    
    % --- 图 2: 互补性展示 (边界校验) ---
    fig_comp = figure('Name', '调度指令 vs 可靠边界', 'Color', 'w', 'Position', [200, 200, 1000, 600]);
    hold on;
    % 将 P_AC_opt 和 P_EV_opt 也转为带符号的用于 plot 绘制
    P_AC_line = P_AC_opt .* direction_sign;
    P_EV_line = P_EV_opt .* direction_sign;

    plot(t_axis, Reliable_AC_plot, 'b:', 'LineWidth', 1.5, 'DisplayName', 'AC 可靠边界');
    plot(t_axis, P_AC_line, 'b-', 'LineWidth', 2.0, 'DisplayName', 'AC 调度');
    
    plot(t_axis, Reliable_EV_plot, 'g:', 'LineWidth', 1.5, 'DisplayName', 'EV 可靠边界');
    plot(t_axis, P_EV_line, 'g-', 'LineWidth', 2.0, 'DisplayName', 'EV 调度');
    
    yline(0, 'k-', 'HandleVisibility', 'off');
    ylabel('调节功率 (kW)', 'FontSize', 14, 'FontName', 'Microsoft YaHei'); 
    xlabel('时间', 'FontSize', 14, 'FontName', 'Microsoft YaHei');
    legend('Location', 'best', 'FontSize', 12, 'FontName', 'Microsoft YaHei'); 
    grid on; 
    xlim([6, 30]); 
    set(gca, 'XTick', x_ticks_set, 'XTickLabel', x_labels_set, ...
        'FontSize', 12, 'FontName', 'Microsoft YaHei', 'LineWidth', 1.2);
    
    print(fig_comp, 'AC与EV时序调度量对比_Optimized.png', '-dpng', '-r600');
end

%% ================= 场景 C: 网络阻塞管理 =================
fprintf('\n>>> 场景 C: 网络阻塞管理测试 <<<\n');
net_params_C = net_params;
limit_val = mean(P_grid_demand) * 0.2; 
net_params_C.LineLimit(1) = limit_val; 
fprintf('  - 制造阻塞: 线路1限额 %.2f kW\n', limit_val);

risk_p.beta = 10; 
[H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_robust_network_improve(...
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
    plot(t_axis, Flow_L1, 'b-', 'LineWidth', 1.5); hold on;
    yline(limit_val, 'r--', 'LineWidth', 2, 'Label', '线路限额');
    ylabel('线路1 潮流影响 (kW)'); xlabel('时间');
    grid on; xlim([6, 30]); set(gca, 'XTick', x_ticks_set, 'XTickLabel', x_labels_set);
    
    print(gcf, '网络阻塞管理测试.png', '-dpng', '-r300');
else
    fprintf('  - 求解失败。\n');
end

%% ================= 场景 D: 鲁棒性测试 =================
fprintf('\n>>> 场景 D: 极端场景鲁棒性测试 <<<\n');
if ~isempty(strategies{1}) && ~isempty(strategies{3})
    Total_Cap_Scen = Scenarios_AC_Up + Scenarios_EV_Up;
    [~, worst_idx] = min(sum(Total_Cap_Scen));
    fprintf('  - 最恶劣场景 (当前混合方向): #%d\n', worst_idx);

    Real_Cap_AC = Scenarios_AC_Up(:, worst_idx);
    Real_Cap_EV = Scenarios_EV_Up(:, worst_idx);

    calc_viol = @(P_ac, P_ev) sum(max(0, P_ac - Real_Cap_AC) + max(0, P_ev - Real_Cap_EV));
    viol_neutral = calc_viol(strategies{1}.P_AC, strategies{1}.P_EV);
    viol_robust  = calc_viol(strategies{3}.P_AC, strategies{3}.P_EV);
    
    fprintf('  - 中性策略(Beta=0) 违约量: %.2f kW\n', viol_neutral);
    fprintf('  - 规避策略(Beta=100) 违约量: %.2f kW\n', viol_robust);
    
    figure('Name', '场景D_鲁棒性', 'Color', 'w', 'Position', [600, 300, 500, 400]);
    b = bar([viol_neutral, viol_robust], 0.5);
    b.FaceColor = 'flat'; b.CData(1,:) = [0.8 0.2 0.2]; b.CData(2,:) = [0.2 0.6 0.2];
    set(gca, 'XTickLabel', {'中性 (\beta=0)', '规避 (\beta=100)'});
    ylabel('极端场景实际违约量 (kW)');
    grid on;
    print(gcf, '极端场景鲁棒性测试.png', '-dpng', '-r300');
end

fprintf('\n所有测试结束。\n');

%% ================= 本地辅助函数 =================
function SDCI = calculate_SDCI_local(n_AC, n_EV, deltaP_AC, deltaP_EV)
    AC = n_AC .* deltaP_AC;
    EV = n_EV .* deltaP_EV;
    min_vals = min(AC, EV);
    max_vals = max(AC, EV);
    total_max = sum(max_vals);
    if total_max == 0
        SDCI = 0;
    else
        SDCI = sum(min_vals) / total_max;
    end
end

function rho = calculate_Rho_local(n_AC, deltaP_AC, n_EV, deltaP_EV)
    AC_series = n_AC .* deltaP_AC;
    EV_series = n_EV .* deltaP_EV;
    if length(AC_series) > 2
        rho = corr(AC_series, EV_series, 'Type', 'Spearman');
    else
        rho = 0;
    end
end