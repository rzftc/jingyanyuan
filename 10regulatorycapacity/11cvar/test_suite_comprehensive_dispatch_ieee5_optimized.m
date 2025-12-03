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

% --- 时间轴设置 (统一 06:00 - 30:00) ---
if exist('time_points', 'var')
    t_axis = time_points;
else
    dt_sim = 5/60; 
    t_axis = 6 : dt_sim : (6 + (T_steps-1)*dt_sim);
end
x_ticks_set = [6, 12, 18, 24, 30];
x_labels_set = {'06:00', '12:00', '18:00', '24:00', '06:00(+1)'};

% --- 物理边界设置 ---
Physical_AC_Up = max(Scenarios_AC_Up, [], 2);
Physical_EV_Up = max(Scenarios_EV_Up, [], 2);
Reliable_AC_Up = Reliable_AC_Up(:);
Reliable_EV_Up = Reliable_EV_Up(:);

% --- 基础成本参数 ---
cost_params.c1_ac = 0.5;  cost_params.c2_ac = 0.01;
cost_params.c1_ev = 0.4;  cost_params.c2_ev = 0.01;
cost_params.c_slack = 200; % 切负荷惩罚

% --- [新增] 互补性与相关性优化权重 ---
% 这些参数控制对 SDCI 和 Rho 的优化力度
lambda_SDCI = 10;   % 互补性惩罚权重 (针对同向重叠)
lambda_Rho  = 10;   % 相关性惩罚权重 (针对同向趋势)
Max_Iter    = 3;    % 迭代次数 (通常3-5次即可收敛)

% --- 构造风险权衡需求 ---
rng(105); 
Total_Reliable = Reliable_AC_Up + Reliable_EV_Up;
Total_Physical = Physical_AC_Up + Physical_EV_Up;
P_grid_demand = Total_Reliable + 0.40 * (Total_Physical - Total_Reliable) .* rand(T_steps, 1);

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

fprintf('初始化完成。最大物理容量: %.2f kW\n', max(Total_Physical));

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
    
    % --- [核心改造开始]：迭代优化循环 ---
    % 初始化上一轮的功率结果 (用于计算惩罚项)
    P_AC_prev = zeros(T_steps, 1);
    P_EV_prev = zeros(T_steps, 1);
    
    for iter = 1:Max_Iter
        % 1. 构建基础 QP 模型 (成本 + 风险 + 物理约束)
        [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_robust_network_improve(...
            P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
            Physical_AC_Up, Physical_EV_Up, ... 
            cost_params, risk_p, net_params);
        
        % 2. [新增] 动态添加互补性与相关性惩罚项
        if iter > 1
            % 策略：如果上一轮 EV 出力高，则本轮增加 AC 在该时刻的成本 (迫使错峰)
            % SDCI 惩罚向量: lambda * P_other
            penalty_vec_AC = lambda_SDCI * P_EV_prev;
            penalty_vec_EV = lambda_SDCI * P_AC_prev;
            
            % Rho 惩罚向量: lambda * (P_other - mean) (抑制同步波动)
            P_EV_centered = P_EV_prev - mean(P_EV_prev);
            P_AC_centered = P_AC_prev - mean(P_AC_prev);
            penalty_rho_AC = lambda_Rho * P_EV_centered;
            penalty_rho_EV = lambda_Rho * P_AC_centered;
            
            % 将惩罚叠加到线性成本向量 f 中
            f(info.idx_P_AC) = f(info.idx_P_AC) + penalty_vec_AC + penalty_rho_AC;
            f(info.idx_P_EV) = f(info.idx_P_EV) + penalty_vec_EV + penalty_rho_EV;
            
            % fprintf('    Iter %d: 更新互补性惩罚...\n', iter);
        end
        
        % 3. 求解 QP
        [x_opt, fval, exitflag] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
        
        if exitflag > 0
            % 提取当前结果
            P_AC_curr = x_opt(info.idx_P_AC);
            P_EV_curr = x_opt(info.idx_P_EV);
            P_Slack = x_opt(info.idx_Slack);
            eta_val = x_opt(info.idx_eta);
            z_val = x_opt(info.idx_z);
            
            % 更新上一轮结果
            P_AC_prev = P_AC_curr;
            P_EV_prev = P_EV_curr;

            % 计算真实成本 (不含惩罚项)
            cost_gen = sum(cost_params.c1_ac*P_AC_curr + cost_params.c2_ac*P_AC_curr.^2 + ...
                           cost_params.c1_ev*P_EV_curr + cost_params.c2_ev*P_EV_curr.^2);
            cost_slack = sum(cost_params.c_slack * abs(P_Slack));
            total_real_cost = cost_gen + cost_slack;
            
            % CVaR 风险
            cvar_val = eta_val + (1 / (N_scenarios * (1 - risk_p.confidence))) * sum(z_val);
            
            % [新增] 计算最终互补性指标
            n_dummy = ones(T_steps, 1);
            val_SDCI = calculate_SDCI_local(n_dummy, n_dummy, P_AC_curr, P_EV_curr);
            val_Rho  = calculate_Rho_local(n_dummy, P_AC_curr, n_dummy, P_EV_curr);

            b_run_cost(i) = cost_gen;      
            b_slack_sum(i) = sum(abs(P_Slack));
            b_risk_val(i) = cvar_val;
            
            strategies{i}.P_AC = P_AC_curr;
            strategies{i}.P_EV = P_EV_curr;
            strategies{i}.P_Slack = P_Slack; 
            
            % ===================================================================
            % 【修改点】：仅在最后一次迭代输出指定格式的结果
            % ===================================================================
            if iter == Max_Iter
                fprintf('发电成本: %.2f, 切负荷量: %.2f, 总运行成本: %.2f, CVaR风险: %.2f, rho: %.4f, sdci: %.4f\n', ...
                    cost_gen, b_slack_sum(i), total_real_cost, cvar_val, val_Rho, val_SDCI);
            end
            % ===================================================================

        else
            fprintf('    Iter %d: 求解失败 (Exitflag %d)\n', iter, exitflag);
            break; 
        end
    end
    % --- [核心改造结束] ---
end

% --- 绘图 B: 风险偏好与互补性效果 ---
if any(~isnan(b_run_cost))
    figure('Name', '场景B_风险灵敏度', 'Color', 'w', 'Position', [100, 100, 900, 400]);
    
    yyaxis left; 
    bar(1:3, b_slack_sum, 0.5, 'FaceColor', [0.8 0.3 0.3]); 
    ylabel('总切负荷量 (kW)');
    set(gca, 'XTick', 1:3, 'XTickLabel', beta_values);
    
    yyaxis right; 
    plot(1:3, b_risk_val, 'b-o', 'LineWidth', 2, 'MarkerSize', 8); 
    ylabel('CVaR 风险 (kW)');
    
    xlabel('风险厌恶系数 \beta');
    legend('切负荷量 (安全性)', '潜在违约风险 (经济性)', 'Location', 'best');
    grid on;
    
    print(gcf, '风险偏好灵敏度分析.png', '-dpng', '-r300');
end

%% ================= [新增] 绘制详细调度方案 (Beta=10) =================
fprintf('\n>>> 正在绘制详细调度方案 (选择 Beta=10) ...\n');
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
    
    % --- 图 1: 功率堆叠 (独立窗口) ---
    fig_stack = figure('Name', '协同调度功率堆叠', 'Color', 'w', 'Position', [150, 150, 1000, 600]);
    hold on;
    Y_Stack = [P_AC_opt, P_EV_opt, P_Slack_calc];
    h_area = area(t_axis, Y_Stack);
    
    h_area(1).FaceColor = [0.00, 0.45, 0.74]; h_area(1).EdgeColor = 'none'; % AC (蓝色)
    h_area(2).FaceColor = [0.47, 0.67, 0.19]; h_area(2).EdgeColor = 'none'; % EV (绿色)
    h_area(3).FaceColor = [0.85, 0.33, 0.10]; h_area(3).EdgeColor = 'none'; h_area(3).FaceAlpha = 0.8; % Slack (橙红)
    
    plot(t_axis, P_grid_demand, 'k--', 'LineWidth', 2.0, 'DisplayName', '电网总需求');
    
    ylabel('功率 (kW)', 'FontSize', 14, 'FontName', 'Microsoft YaHei');
    % title(...); % [已移除标题]
    
    legend([h_area(1), h_area(2), h_area(3)], {'空调 (AC)', '电动汽车 (EV)', '缺额 (Slack)'}, ...
           'Location', 'northwest', 'FontSize', 12, 'FontName', 'Microsoft YaHei');
    grid on; 
    xlim([6, 30]); 
    set(gca, 'XTick', x_ticks_set, 'XTickLabel', x_labels_set, ...
        'FontSize', 12, 'FontName', 'Microsoft YaHei', 'LineWidth', 1.2);
    
    % 保存图 1
    print(fig_stack, '协同调度功率堆叠_Optimized.png', '-dpng', '-r600');
    fprintf('  >>> 图1已保存: 协同调度功率堆叠_Optimized.png\n');
    
    
    % --- 图 2: 互补性展示 (AC vs EV 曲线) (独立窗口) ---
    fig_comp = figure('Name', 'AC与EV时序出力对比', 'Color', 'w', 'Position', [200, 200, 1000, 600]);
    hold on;
    plot(t_axis, P_AC_opt, '-', 'Color', [0.00, 0.45, 0.74], 'LineWidth', 2.5, 'DisplayName', 'AC 出力');
    plot(t_axis, P_EV_opt, '-', 'Color', [0.47, 0.67, 0.19], 'LineWidth', 2.5, 'DisplayName', 'EV 出力');
    
    ylabel('调节功率 (kW)', 'FontSize', 14, 'FontName', 'Microsoft YaHei'); 
    xlabel('时间', 'FontSize', 14, 'FontName', 'Microsoft YaHei');
    % title(...); % [已移除标题]
    
    legend('Location', 'best', 'FontSize', 12, 'FontName', 'Microsoft YaHei'); 
    grid on; 
    xlim([6, 30]); 
    set(gca, 'XTick', x_ticks_set, 'XTickLabel', x_labels_set, ...
        'FontSize', 12, 'FontName', 'Microsoft YaHei', 'LineWidth', 1.2);
    
    % 保存图 2
    print(fig_comp, 'AC与EV时序出力对比_Optimized.png', '-dpng', '-r600');
    fprintf('  >>> 图2已保存: AC与EV时序出力对比_Optimized.png\n');
end
%% ================= 场景 C: 网络阻塞管理 =================
fprintf('\n>>> 场景 C: 网络阻塞管理测试 <<<\n');

net_params_C = net_params;
limit_val = mean(P_grid_demand) * 0.2; 
net_params_C.LineLimit(1) = limit_val; 
fprintf('  - 制造阻塞: 线路1限额 %.2f kW\n', limit_val);

risk_p.beta = 10; 
% 这里简单运行一次 QP，不做迭代，主要看网络约束是否生效
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
    yline(-limit_val, 'r--', 'LineWidth', 2);
    ylabel('线路1 潮流 (kW)'); xlabel('时间');
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
    fprintf('  - 最恶劣场景: #%d\n', worst_idx);

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