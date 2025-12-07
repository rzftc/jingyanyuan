%% 修正说明：
% 1. [核心修正] 成本计算引入时间积分 (dt)，将单位从 元/MW 修正为 元/MWh。
% 2. [核心修正] 调整二次成本系数 c2，使其符合物理实际，避免天文数字成本。
% 3. [核心修正] 将功率单位统一为 MW。
% 4. 物理边界、需求生成、绘图逻辑保持不变。

clear; close all; clc;

%% ================= 1. 全局初始化 =================
fprintf('正在加载场景数据...\n');
data_file = 'reliable_regulation_domain_1000_mix_01_new2.mat';
if ~exist(data_file, 'file')
    error('数据文件缺失！请先运行 main_scenario_generation_diff.m');
end
load(data_file);

% 数据维度对齐
[T_steps, N_scenarios] = size(Scenarios_AC_Up);

% --- [修改] 时间参数 ---
if exist('time_points', 'var')
    t_axis = time_points;
else
    dt_sim = 5/60; 
    t_axis = 6 : dt_sim : (6 + (T_steps-1)*dt_sim);
end
dt = 5/60; % 时间步长 (小时)
fprintf('时间步长 dt = %.4f 小时\n', dt);

% --- [修改] 单位转换系数 (kW -> MW) ---
unit_scale = 1/1000; 

% --- 修正1: 物理边界设置 (提取双向边界) & [修改] 单位转换 ---
Scenarios_AC_Up = Scenarios_AC_Up * unit_scale; 
Scenarios_EV_Up = Scenarios_EV_Up * unit_scale; 
Physical_AC_Up = max(Scenarios_AC_Up, [], 2);
Physical_EV_Up = max(Scenarios_EV_Up, [], 2);
Reliable_AC_Up = Reliable_AC_Up(:) * unit_scale;
Reliable_EV_Up = Reliable_EV_Up(:) * unit_scale;

Scenarios_AC_Down = Scenarios_AC_Down * unit_scale; 
Scenarios_EV_Down = Scenarios_EV_Down * unit_scale; 
Physical_AC_Down = max(abs(Scenarios_AC_Down), [], 2);
Physical_EV_Down = max(abs(Scenarios_EV_Down), [], 2);
Reliable_AC_Down = abs(Reliable_AC_Down(:)) * unit_scale;
Reliable_EV_Down = abs(Reliable_EV_Down(:)) * unit_scale;

% --- [修改] 成本参数 (符合中国电力市场现状) ---
% 说明：成本计算公式将改为 sum( (c1*P + c2*P^2) * dt )
% 单位：元/MWh
cost_params.c1_ac = 500;      % 500 元/MWh (0.5元/kWh) - 空调调节成本
cost_params.c2_ac = 50;       % 50 元/(MW)^2h - 二次项系数大幅降低，避免成本失真
cost_params.c1_ev = 400;      % 400 元/MWh (0.4元/kWh) - EV调节成本(略低于AC)
cost_params.c2_ev = 50;       % 50 元/(MW)^2h
cost_params.c_slack = 30000;  % 30,000 元/MWh (30元/kWh) - 切负荷惩罚 (VoLL)

% --- 互补性与相关性优化权重 ---
lambda_SDCI = 10;   
lambda_Rho  = 10;   
Max_Iter    = 3;    

% --- 修正2: 构造混合需求 ---
rng(105); 
P_grid_demand = zeros(T_steps, 1);
Effective_Scen_AC = zeros(T_steps, N_scenarios);
Effective_Scen_EV = zeros(T_steps, N_scenarios);
Effective_Phys_AC = zeros(T_steps, 1);
Effective_Phys_EV = zeros(T_steps, 1);
Effective_Reliable_AC = zeros(T_steps, 1); 
Effective_Reliable_EV = zeros(T_steps, 1); 

direction_signal = randn(T_steps, 1); 

for t = 1:T_steps
    is_up_regulation = direction_signal(t) >= 0;
    
    if is_up_regulation
        % 上调
        cap_rel = Reliable_AC_Up(t) + Reliable_EV_Up(t);
        cap_phy = Physical_AC_Up(t) + Physical_EV_Up(t);
        Effective_Scen_AC(t,:) = Scenarios_AC_Up(t,:);
        Effective_Scen_EV(t,:) = Scenarios_EV_Up(t,:);
        Effective_Phys_AC(t) = Physical_AC_Up(t);
        Effective_Phys_EV(t) = Physical_EV_Up(t);
        Effective_Reliable_AC(t) = Reliable_AC_Up(t);
        Effective_Reliable_EV(t) = Reliable_EV_Up(t);
    else
        % 下调
        cap_rel = Reliable_AC_Down(t) + Reliable_EV_Down(t);
        cap_phy = Physical_AC_Down(t) + Physical_EV_Down(t);
        Effective_Scen_AC(t,:) = abs(Scenarios_AC_Down(t,:));
        Effective_Scen_EV(t,:) = abs(Scenarios_EV_Down(t,:));
        Effective_Phys_AC(t) = Physical_AC_Down(t);
        Effective_Phys_EV(t) = Physical_EV_Down(t);
        Effective_Reliable_AC(t) = Reliable_AC_Down(t);
        Effective_Reliable_EV(t) = Reliable_EV_Down(t);
    end
    
    % 需求生成：在可靠容量基础上增加波动，模拟高风险场景
    req_magnitude = cap_rel + 0.40 * (cap_phy - cap_rel) * rand();
    P_grid_demand(t) = req_magnitude;
end

Scenarios_AC_Up = Effective_Scen_AC;
Scenarios_EV_Up = Effective_Scen_EV;
Physical_AC_Up  = Effective_Phys_AC;
Physical_EV_Up  = Effective_Phys_EV;
Reliable_AC_Up  = Effective_Reliable_AC;
Reliable_EV_Up  = Effective_Reliable_EV;

% --- 网络参数 ---
N_bus = 5; N_line = 6;
LineData = [
    1, 2, 0.0281; 1, 4, 0.0304; 1, 5, 0.0064;
    2, 3, 0.0108; 3, 4, 0.0297; 4, 5, 0.0297
];
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

% 1. 定义标准工况 (MW)
Base_Gen  = [250;    0;    0;    0;  750]; % G1, G5 发电
Base_Load = [  0;  300;  300;  400;    0]; % L2, L3, L4 负荷
P_inj_peak = Base_Gen - Base_Load;

% 2. 定义日负荷曲线形状 (峰值在第19小时)
% T_steps 已经在前面定义
t_vec = linspace(0, 24, T_steps);
load_curve = 0.55 + 0.45 * exp(-((t_vec - 19).^2) / 12) + 0.1 * exp(-((t_vec - 10).^2) / 20);
load_curve = load_curve / max(load_curve); % 归一化

% 3. 计算时变基础潮流
Flow_peak = net_params.PTDF * P_inj_peak; % (N_line x 1)
net_params.BaseFlow = Flow_peak * load_curve;  % (N_line x T)
% --- [关键修正：动态设置线路限额] ---
% 计算所有时刻、所有线路的最大背景潮流
max_base_flow = max(max(abs(net_params.BaseFlow)));
fprintf('系统最大背景潮流: %.2f MW\n', max_base_flow);

% 将线路限额设定为：max(250, 最大背景潮流 * 1.2)
% 这样既保留了原定的250MW下限，又防止了背景潮流导致的初始越限
% safe_limit = max(250, max_base_flow * 1.2); 
safe_limit = max_base_flow + 2.0;
net_params.LineLimit = ones(N_line, 1) * safe_limit;

fprintf('初始化完成。\n');

%% ================= 场景 B: 风险偏好灵敏度分析 =================
fprintf('\n>>> 场景 B: 风险偏好灵敏度分析 <<<\n');
x_ticks_set = [6, 12, 18, 24, 30];
x_labels_set = {'06:00', '12:00', '18:00', '24:00', '06:00(+1)'};
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
    risk_p.rho_pen = 50 * 1000 * unit_scale; % 保持惩罚的相对量级 (这里直接设为50也行，取决于对风险项的定义)
    % 修正：为了让风险项与新的成本项(元)匹配，rho_pen 需要重新标定。
    % 假设 rho_pen 是“每MW越限的惩罚价格”，设为 5000 元/MW 比较合适(略高于成本，低于切负荷)
    risk_p.rho_pen = 5000; 
    
    P_AC_prev = zeros(T_steps, 1);
    P_EV_prev = zeros(T_steps, 1);
    
    for iter = 1:Max_Iter
        % 1. 构建模型 (注意：construct函数内部不包含dt乘法，需要在外部处理f和H)
        [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast(...
            P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
            Physical_AC_Up, Physical_EV_Up, ... 
            cost_params, risk_p, net_params);
        
        % --- [关键修改] 乘以 dt，将功率成本转化为电能量成本 ---
        % 对 H 中的功率变量对应的对角元素 * dt
        % 对 f 中的功率变量对应的元素 * dt
        
        % 功率变量索引
        idx_pow = [info.idx_P_AC, info.idx_P_EV, info.idx_Slack];
        
        % 修正 H (二次项: c2 * P^2 * dt)
        for idx = idx_pow
             H(idx, idx) = H(idx, idx) * dt;
        end
        
        % 修正 f (一次项: c1 * P * dt)
        f(idx_pow) = f(idx_pow) * dt;
        
        % 2. 添加互补性惩罚 (惩罚项也建议乘以 dt 以保持量级一致)
        if iter > 1
            penalty_vec_AC = lambda_SDCI * P_EV_prev * dt; % * dt
            penalty_vec_EV = lambda_SDCI * P_AC_prev * dt; % * dt
            
            P_EV_centered = P_EV_prev - mean(P_EV_prev);
            P_AC_centered = P_AC_prev - mean(P_AC_prev);
            penalty_rho_AC = lambda_Rho * P_EV_centered * dt; % * dt
            penalty_rho_EV = lambda_Rho * P_AC_centered * dt; % * dt
            
            f(info.idx_P_AC) = f(info.idx_P_AC) + penalty_vec_AC + penalty_rho_AC;
            f(info.idx_P_EV) = f(info.idx_P_EV) + penalty_vec_EV + penalty_rho_EV;
        end
        
        % 3. 求解
        [x_opt, fval, exitflag] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
        
        if exitflag > 0
            P_AC_curr = x_opt(info.idx_P_AC);
            P_EV_curr = x_opt(info.idx_P_EV);
            P_Slack = x_opt(info.idx_Slack);
            eta_val = x_opt(info.idx_eta);
            z_val = x_opt(info.idx_z);
            
            P_AC_prev = P_AC_curr;
            P_EV_prev = P_EV_curr;

            % 计算真实成本 (含 dt)
            cost_gen = sum((cost_params.c1_ac*P_AC_curr + cost_params.c2_ac*P_AC_curr.^2) * dt + ...
                           (cost_params.c1_ev*P_EV_curr + cost_params.c2_ev*P_EV_curr.^2) * dt);
            cost_slack = sum(cost_params.c_slack * abs(P_Slack) * dt);
            total_real_cost = cost_gen + cost_slack;
            
            cvar_val = eta_val + (1 / (N_scenarios * (1 - risk_p.confidence))) * sum(z_val);
            
            % 记录指标
            n_dummy = ones(T_steps, 1);
            val_SDCI = calculate_SDCI_local(n_dummy, n_dummy, P_AC_curr, P_EV_curr);
            val_Rho  = calculate_Rho_local(n_dummy, P_AC_curr, n_dummy, P_EV_curr);

            b_run_cost(i) = cost_gen;      
            b_slack_sum(i) = sum(abs(P_Slack)) * dt; % 转换为 MWh
            b_risk_val(i) = cvar_val; % CVaR 本质是功率值(MW)的统计量，通常不乘dt，或者乘dt变为能量风险
            
            strategies{i}.P_AC = P_AC_curr;
            strategies{i}.P_EV = P_EV_curr;
            strategies{i}.P_Slack = P_Slack; 
            
            if iter == 1
                strategies{i}.SDCI_History = zeros(Max_Iter, 1);
                strategies{i}.Rho_History = zeros(Max_Iter, 1);
            end
            strategies{i}.SDCI_History(iter) = val_SDCI;
            strategies{i}.Rho_History(iter) = val_Rho;

            if iter == Max_Iter
                fprintf('发电成本: %.2f (元), 切负荷量: %.2f (MWh), 总成本: %.2f (元), CVaR风险: %.2f (MW), rho: %.4f, sdci: %.4f\n', ...
                    cost_gen, b_slack_sum(i), total_real_cost, cvar_val, val_Rho, val_SDCI);
            end
        else
            fprintf('    Iter %d: 求解失败 (Exitflag %d)\n', iter, exitflag);
            break; 
        end
    end
end

%% 
% --- 绘图 B: 风险偏好灵敏度分析 (含数据标注) ---
if any(~isnan(b_run_cost))
    figure('Name', '场景B_风险灵敏度', 'Color', 'w', 'Position', [100, 100, 900, 400]);
    
    % --- 左轴：切负荷 (柱状图) ---
    yyaxis left; 
    b = bar(1:3, b_slack_sum, 0.5, 'FaceColor', [0.8 0.3 0.3]); 
    ylabel('总切负荷电量 (MWh)'); 
    set(gca, 'XTick', 1:3, 'XTickLabel', beta_values);
    % 添加左轴数据标注
    for i = 1:length(b_slack_sum)
        text(i, b_slack_sum(i), sprintf('%.2f', b_slack_sum(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 12, 'Color', [0.6 0.1 0.1], 'FontWeight', 'bold');
    end
    
    % --- 右轴：风险 (折线图) ---
    yyaxis right; 
    plot(1:3, b_risk_val, 'b-o', 'LineWidth', 2, 'MarkerSize', 8); 
    ylabel('CVaR 潜在违约风险 (MW)'); 
    % 添加右轴数据标注
    for i = 1:length(b_risk_val)
        text(i, b_risk_val(i), sprintf('%.2f', b_risk_val(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 12, 'Color', 'b', 'FontWeight', 'bold');
    end
    
    xlabel('风险厌恶系数 \beta');
    legend('切负荷量 (安全性)', '潜在违约风险 (经济性)', 'Location', 'best');
    grid on;
    print(gcf, '风险偏好灵敏度分析.png', '-dpng', '-r300');
end

%% ================= [绘制详细调度方案 (Beta=10)] =================
fprintf('\n>>> 正在绘制详细调度方案 (Beta=10)... \n');
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
    
    direction_sign = sign(direction_signal);
    direction_sign(direction_sign == 0) = 1; 
    
    Y_Stack_Plot = [P_AC_opt, P_EV_opt, P_Slack_calc] .* direction_sign;
    Reliable_AC_plot = Reliable_AC_Up .* direction_sign;
    Reliable_EV_plot = Reliable_EV_Up .* direction_sign;
    P_Demand_plot = P_grid_demand .* direction_sign;
    
    % --- 图 1: 功率堆叠 ---
    fig_stack = figure('Name', '协同调度功率堆叠', 'Color', 'w', 'Position', [150, 150, 1000, 600]);
    hold on;
    h_area = area(t_axis, Y_Stack_Plot);
    h_area(1).FaceColor = [0.00, 0.45, 0.74]; h_area(1).EdgeColor = 'none';
    h_area(2).FaceColor = [0.47, 0.67, 0.19]; h_area(2).EdgeColor = 'none';
    h_area(3).FaceColor = [0.85, 0.33, 0.10]; h_area(3).EdgeColor = 'none'; h_area(3).FaceAlpha = 0.8;
    plot(t_axis, P_Demand_plot, 'k--', 'LineWidth', 2.0, 'DisplayName', '电网总需求');
    ylabel('功率 (MW) [正=上调, 负=下调]', 'FontSize', 14, 'FontName', 'Microsoft YaHei');
    legend([h_area(1), h_area(2), h_area(3)], {'空调 (AC)', '电动汽车 (EV)', '缺额 (Slack)'}, ...
           'Location', 'northwest', 'FontSize', 12, 'FontName', 'Microsoft YaHei');
    grid on; xlim([6, 30]); 
    set(gca, 'XTick', x_ticks_set, 'XTickLabel', x_labels_set, 'FontSize', 12, 'FontName', 'Microsoft YaHei');
    print(fig_stack, '协同调度功率堆叠_Optimized.png', '-dpng', '-r600');
    
    % --- 图 2: 互补性展示 ---
    fig_comp = figure('Name', '调度指令 vs 可靠边界', 'Color', 'w', 'Position', [200, 200, 1000, 600]);
    hold on;
    P_AC_line = P_AC_opt .* direction_sign;
    P_EV_line = P_EV_opt .* direction_sign;
    plot(t_axis, Reliable_AC_plot, 'b:', 'LineWidth', 1.5, 'DisplayName', 'AC 可靠边界');
    plot(t_axis, P_AC_line, 'b-', 'LineWidth', 2.0, 'DisplayName', 'AC 调度');
    plot(t_axis, Reliable_EV_plot, 'g:', 'LineWidth', 1.5, 'DisplayName', 'EV 可靠边界');
    plot(t_axis, P_EV_line, 'g-', 'LineWidth', 2.0, 'DisplayName', 'EV 调度');
    yline(0, 'k-', 'HandleVisibility', 'off');
    ylabel('调节功率 (MW)', 'FontSize', 14, 'FontName', 'Microsoft YaHei'); 
    xlabel('时间', 'FontSize', 14, 'FontName', 'Microsoft YaHei');
    legend('Location', 'best', 'FontSize', 12, 'FontName', 'Microsoft YaHei'); 
    grid on; xlim([6, 30]); 
    set(gca, 'XTick', x_ticks_set, 'XTickLabel', x_labels_set, 'FontSize', 12, 'FontName', 'Microsoft YaHei');
    print(fig_comp, 'AC与EV时序调度量对比_Optimized.png', '-dpng', '-r600');
    
    if isfield(strategies{idx_plot}, 'SDCI_History') && isfield(strategies{idx_plot}, 'Rho_History')
        % --- SDCI 对比图 ---
        fig_sdci = figure('Name', 'SDCI 迭代对比', 'Color', 'w', 'Position', [300, 300, 600, 400]);
        sdci_vals = [strategies{idx_plot}.SDCI_History(1), strategies{idx_plot}.SDCI_History(end)];
        bar(sdci_vals, 0.4, 'FaceColor', [0.2 0.6 0.8]);
        set(gca, 'XTickLabel', {'初始迭代', '最终迭代'}, 'FontSize', 12, 'FontName', 'Microsoft YaHei');
        ylabel('SDCI 指标 (互补性)', 'FontSize', 12, 'FontName', 'Microsoft YaHei');
        grid on;
        text(1:2, sdci_vals, num2str(sdci_vals', '%.4f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
        print(fig_sdci, 'SDCI对比.png', '-dpng', '-r600');

        % --- Rho 对比图 ---
        fig_rho = figure('Name', 'Rho 迭代对比', 'Color', 'w', 'Position', [400, 400, 600, 400]);
        rho_vals = [strategies{idx_plot}.Rho_History(1), strategies{idx_plot}.Rho_History(end)];
        bar(rho_vals, 0.4, 'FaceColor', [0.8 0.4 0.2]);
        set(gca, 'XTickLabel', {'初始迭代', '最终迭代'}, 'FontSize', 12, 'FontName', 'Microsoft YaHei');
        ylabel('Spearman Rho 指标 (相关性)', 'FontSize', 12, 'FontName', 'Microsoft YaHei');
        grid on;
        text(1:2, rho_vals, num2str(rho_vals', '%.4f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
        print(fig_rho, 'Rho对比.png', '-dpng', '-r600');
        
        fprintf('  >>> 迭代对比图已保存: SDCI对比.png, Rho对比.png\n');
    end
end

%% ================= 场景 C: 网络阻塞管理 =================
fprintf('\n>>> 场景 C: 网络阻塞管理测试 <<<\n');
net_params_C = net_params;

% [修正逻辑]：针对有背景潮流的情况制造阻塞
% 选取线路 1 在峰值时刻的背景潮流
peak_flow_L1 = max(abs(net_params.BaseFlow(1, :)));

% 将限额设定为比峰值潮流略小 (例如 90% 或 95%)，从而制造拥堵
% 注意：要确保 VPP 的调节能力(约10-20MW)足够消除这个阻塞，不要设得太低
limit_val = peak_flow_L1 * 0.95; 

% 如果背景潮流本来就很小(例如<5MW)，则使用旧逻辑
if limit_val < 5
     limit_val = mean(P_grid_demand) * 0.5;
end

net_params_C.LineLimit(1) = limit_val; 
fprintf('  - 制造阻塞: 线路1 原峰值流量 %.2f MW -> 新限额 %.2f MW\n', peak_flow_L1, limit_val);

risk_p.beta = 10; 
[H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast(...
    P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, Physical_AC_Up, Physical_EV_Up, ...
    cost_params, risk_p, net_params_C);

% 同样需要应用 dt 修正
idx_pow = [info.idx_P_AC, info.idx_P_EV, info.idx_Slack];
for idx = idx_pow
     H(idx, idx) = H(idx, idx) * dt;
end
f(idx_pow) = f(idx_pow) * dt;

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
    ylabel('线路1 潮流 (MW)'); xlabel('时间');
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

    % [修改] 违约量也改为电量 (MWh) 可能更直观，或者保持 MW 峰值?
    % 这里保持为累积违约功率 (MW) 以便与之前一致，但数值会小很多
    calc_viol = @(P_ac, P_ev) sum(max(0, P_ac - Real_Cap_AC) + max(0, P_ev - Real_Cap_EV));
    
    viol_neutral = calc_viol(strategies{1}.P_AC, strategies{1}.P_EV);
    viol_robust  = calc_viol(strategies{3}.P_AC, strategies{3}.P_EV);
    
    fprintf('  - 中性策略(Beta=0) 累积违约量: %.2f MW\n', viol_neutral);
    fprintf('  - 规避策略(Beta=100) 累积违约量: %.2f MW\n', viol_robust);
    
    figure('Name', '场景D_鲁棒性', 'Color', 'w', 'Position', [600, 300, 500, 400]);
    b = bar([viol_neutral, viol_robust], 0.5);
    b.FaceColor = 'flat'; b.CData(1,:) = [0.8 0.2 0.2]; b.CData(2,:) = [0.2 0.6 0.2];
    set(gca, 'XTickLabel', {'中性 (\beta=0)', '规避 (\beta=100)'});
    ylabel('极端场景累积违约量 (MW)');
    grid on;
    print(gcf, '极端场景鲁棒性测试.png', '-dpng', '-r300');
end

fprintf('\n所有测试结束。\n');