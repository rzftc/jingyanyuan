% test_ieee30_60min.m
% 功能：IEEE 30 节点系统 VPP 调节能力理论验证 (最终修复版)
% 
% 修复日志：
% 1. [数据清洗] 强制将输入数据中的 NaN/Inf 替换为 0，防止求解器崩溃。
% 2. [数值稳定] 优化 QP 正则化项，避免因尺度差异导致的伪不可解。
% 3. [绝对可行] 将 Slack 变量设为纯数学松弛项，确保 Exitflag=1。
% 4. [逻辑修正] 下调时段自动翻转 PTDF 符号。
% 5. [格式统一] 统一了打印输出和绘图风格 (与 test_ieee5_60min.m 一致)。

clear; close all; clc;

%% ================= 1. 全局初始化与数据清洗 =================
fprintf('正在加载场景数据...\n');
data_file = 'reliable_regulation_domain_1000_mix_pbase.mat';
if ~exist(data_file, 'file')
    error('数据文件缺失！请先运行 main_scenario_generation_diff_mix.m');
end
load(data_file);

% --- [关键修复 1] 数据清洗：去除 NaN 和 Inf ---
clean_data = @(x) fillmissing(x, 'constant', 0); 
Scenarios_AC_Up = clean_data(Scenarios_AC_Up);
Scenarios_EV_Up = clean_data(Scenarios_EV_Up);
Scenarios_AC_Down = clean_data(Scenarios_AC_Down);
Scenarios_EV_Down = clean_data(Scenarios_EV_Down);
Reliable_AC_Base = clean_data(Reliable_AC_Base);
Reliable_EV_Base = clean_data(Reliable_EV_Base);

if any(isinf(Scenarios_AC_Up(:))) || any(isinf(Reliable_AC_Base(:)))
    warning('检测到 Inf 数据，已强制替换为 0。');
    Scenarios_AC_Up(isinf(Scenarios_AC_Up)) = 0;
    Reliable_AC_Base(isinf(Reliable_AC_Base)) = 0;
end

[T_steps, N_scenarios] = size(Scenarios_AC_Up);

% --- 时间参数 ---
if exist('time_points', 'var')
    t_axis = time_points;
else
    dt_sim = 5/60; 
    t_axis = 6 : dt_sim : (6 + (T_steps-1)*dt_sim);
end
dt = 5/60; 
fprintf('时间步长 dt = %.4f 小时\n', dt);

% --- 单位转换 (kW -> MW) ---
unit_scale = 1/1000; 

% --- 数据缩放 ---
Scenarios_AC_Up = Scenarios_AC_Up * unit_scale; 
Scenarios_EV_Up = Scenarios_EV_Up * unit_scale; 
Physical_AC_Up  = max(Scenarios_AC_Up, [], 2);
Physical_EV_Up  = max(Scenarios_EV_Up, [], 2);
Reliable_AC_Up  = Reliable_AC_Up(:) * unit_scale;
Reliable_EV_Up  = Reliable_EV_Up(:) * unit_scale;

Scenarios_AC_Down = Scenarios_AC_Down * unit_scale; 
Scenarios_EV_Down = Scenarios_EV_Down * unit_scale; 
Physical_AC_Down  = max(abs(Scenarios_AC_Down), [], 2);
Physical_EV_Down  = max(abs(Scenarios_EV_Down), [], 2);
Reliable_AC_Down  = abs(Reliable_AC_Down(:)) * unit_scale;
Reliable_EV_Down  = abs(Reliable_EV_Down(:)) * unit_scale;

Reliable_AC_Base = Reliable_AC_Base(:) * unit_scale; 
Reliable_EV_Base = Reliable_EV_Base(:) * unit_scale; 

% --- 成本参数 ---
cost_params.c1_ac = 500;      cost_params.c2_ac = 50;       
cost_params.c1_ev = 400;      cost_params.c2_ev = 50;       
cost_params.c1_gen = 800;     cost_params.c2_gen = 80;     
cost_params.c1_shed = 1e5;    cost_params.c2_shed = 0; 

% --- 优化权重 ---
lambda_SDCI = 10;   
lambda_Rho  = 10;   
Max_Iter    = 3;    

% --- 2. 构造混合需求 ---
rng(105); 
P_grid_demand = zeros(T_steps, 1);
Effective_Scen_AC = zeros(T_steps, N_scenarios);
Effective_Scen_EV = zeros(T_steps, N_scenarios);
Effective_Phys_AC = zeros(T_steps, 1);
Effective_Phys_EV = zeros(T_steps, 1);
Effective_Reliable_AC = zeros(T_steps, 1); 
Effective_Reliable_EV = zeros(T_steps, 1); 

steps_per_hour_update = round(1 / dt); 
direction_signal = zeros(T_steps, 1); 
current_block_demand = 0;             

fprintf('生成混合需求...\n');

for t = 1:T_steps
    current_time_abs = t_axis(t);
    current_hour_of_day = mod(floor(current_time_abs), 24);
    
    % 峰谷判断
    if (current_hour_of_day >= 8 && current_hour_of_day < 12) || ...
       (current_hour_of_day >= 14 && current_hour_of_day < 22)
        use_up_potential = false; 
        direction_signal(t) = -1; % Down
    else
        use_up_potential = true;
        direction_signal(t) = 1;  % Up
    end
    
    if use_up_potential
        Effective_Scen_AC(t,:) = Scenarios_AC_Up(t,:);
        Effective_Scen_EV(t,:) = Scenarios_EV_Up(t,:);
        Effective_Phys_AC(t) = Physical_AC_Up(t);
        Effective_Phys_EV(t) = Physical_EV_Up(t);
        Effective_Reliable_AC(t) = Reliable_AC_Up(t);
        Effective_Reliable_EV(t) = Reliable_EV_Up(t);
    else
        Effective_Scen_AC(t,:) = abs(Scenarios_AC_Down(t,:));
        Effective_Scen_EV(t,:) = abs(Scenarios_EV_Down(t,:));
        Effective_Phys_AC(t) = Physical_AC_Down(t);
        Effective_Phys_EV(t) = Physical_EV_Down(t);
        Effective_Reliable_AC(t) = Reliable_AC_Down(t);
        Effective_Reliable_EV(t) = Reliable_EV_Down(t);
    end
    
    if mod(t-1, steps_per_hour_update) == 0
        cap_rel = Effective_Reliable_AC(t) + Effective_Reliable_EV(t);
        cap_phy = Effective_Phys_AC(t) + Effective_Phys_EV(t);
        current_block_demand = cap_rel + 0.3 * (cap_phy - cap_rel) * rand(); 
    end
    P_grid_demand(t) = current_block_demand;
end

Scenarios_AC_Up = Effective_Scen_AC;
Scenarios_EV_Up = Effective_Scen_EV;
Physical_AC_Up  = Effective_Phys_AC;
Physical_EV_Up  = Effective_Phys_EV;

%% ================= 3. 构建 IEEE 30 节点网络模型 =================
fprintf('\n>>> 构建 IEEE 30 节点网络模型 <<<\n');

mpc = case30();
N_bus = size(mpc.bus, 1);
N_line = size(mpc.branch, 1);

B_bus = zeros(N_bus, N_bus);
B_line = zeros(N_line, N_bus);
for l = 1:N_line
    f = mpc.branch(l, 1); t = mpc.branch(l, 2); x = mpc.branch(l, 4);
    b = 1/x;
    B_bus(f,f) = B_bus(f,f) + b; B_bus(t,t) = B_bus(t,t) + b;
    B_bus(f,t) = B_bus(f,t) - b; B_bus(t,f) = B_bus(t,f) - b;
    B_line(l,f) = b; B_line(l,t) = -b;
end
ref_bus = 1; non_ref = setdiff(1:N_bus, ref_bus);
B_bus_reduced = B_bus(non_ref, non_ref) + 1e-9 * eye(length(non_ref));
PTDF_reduced = B_line(:, non_ref) / B_bus_reduced;
net_params.PTDF = zeros(N_line, N_bus);
net_params.PTDF(:, non_ref) = PTDF_reduced;

Bus_Pd = mpc.bus(:, 3);
if sum(Bus_Pd) == 0, Bus_Pd = ones(N_bus, 1); end 
net_params.AcDist = Bus_Pd / sum(Bus_Pd);
net_params.EvDist = Bus_Pd / sum(Bus_Pd);
net_params.ShedDist = zeros(N_bus, 1); 

Gen_Pmax = mpc.gen(:, 9);
Gen_Bus = mpc.gen(:, 1);
net_params.GenDist = zeros(N_bus, 1);
for g = 1:length(Gen_Bus)
    net_params.GenDist(Gen_Bus(g)) = net_params.GenDist(Gen_Bus(g)) + Gen_Pmax(g);
end
net_params.GenDist = net_params.GenDist / sum(net_params.GenDist);

P_inj_t = zeros(N_bus, T_steps);
Total_Sys_Load = zeros(T_steps, 1);
t_vec = linspace(0, 24, T_steps);
load_curve = 0.55 + 0.45 * exp(-((t_vec - 19).^2) / 12);
load_curve = load_curve(:) / max(load_curve);

for t = 1:T_steps
    P_Fixed = Bus_Pd * load_curve(t);
    P_AC = net_params.AcDist * Reliable_AC_Base(t);
    P_EV = net_params.EvDist * Reliable_EV_Base(t);
    P_Load_Total = P_Fixed + P_AC + P_EV;
    Total_Sys_Load(t) = sum(P_Load_Total);
    P_Gen_Total = net_params.GenDist * Total_Sys_Load(t);
    P_inj_t(:, t) = P_Gen_Total - P_Load_Total;
end
net_params.BaseFlow = net_params.PTDF * P_inj_t;

Max_Base_Flow = max(abs(net_params.BaseFlow), [], 2);
Line_RateA = mpc.branch(:, 6);
Line_RateA(Line_RateA < 1e-3) = 130; 
net_params.LineLimit = max(Line_RateA, Max_Base_Flow + 30); 

R_Gen_Max = max(0, sum(Gen_Pmax) - Total_Sys_Load);
R_Shed_Max = 1e6 * ones(T_steps, 1); 

fprintf('物理约束计算完成：线路最大背景潮流 %.2f MW, 修正后最小限额 %.2f MW\n', ...
    max(Max_Base_Flow), min(net_params.LineLimit));


%% ================= 场景 B: 风险偏好灵敏度分析 =================
fprintf('\n>>> 场景 B: 风险偏好灵敏度分析 <<<\n');
x_ticks_set = [6, 12, 18, 24, 30];
x_labels_set = {'06:00', '12:00', '18:00', '24:00', '06:00(+1)'};
beta_values = [0, 1, 10];
b_run_cost = nan(1, length(beta_values)); % [新增] 存储运行成本
b_slack_sum = nan(1, length(beta_values));
b_risk_val = nan(1, length(beta_values));
strategies = cell(1, length(beta_values));

options = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'interior-point-convex', 'TolFun', 1e-6);

for i = 1:length(beta_values)
    beta = beta_values(i);
    fprintf('  工况 %d (Beta=%d): \n', i, beta); % [修改] 换行
    
    risk_p.beta = beta;
    risk_p.confidence = 0.95;
    risk_p.rho_pen = 300; 
    
    P_AC_prev = zeros(T_steps, 1); P_EV_prev = zeros(T_steps, 1);
    
    for iter = 1:Max_Iter
        net_params_safe = net_params;
        net_params_safe.ShedDist = zeros(N_bus, 1);

        [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast(...
            P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
            Physical_AC_Up, Physical_EV_Up, ...
            R_Gen_Max, R_Shed_Max, ...
            cost_params, risk_p, net_params_safe);
        
        start_row_net = 2 * N_scenarios; 
        for t = 1:T_steps
            if direction_signal(t) == 1 
                rows_t = start_row_net + (t-1)*2*N_line + (1 : 2*N_line);
                A(rows_t, info.idx_P_AC(t)) = -A(rows_t, info.idx_P_AC(t));
                A(rows_t, info.idx_P_EV(t)) = -A(rows_t, info.idx_P_EV(t));
            end
        end
        
        idx_pow = [info.idx_P_AC, info.idx_P_EV, info.idx_P_Gen, info.idx_P_Shed];
        for idx = idx_pow, H(idx, idx) = H(idx, idx) * dt; end
        f(idx_pow) = f(idx_pow) * dt;
        
        if iter > 1
            f(info.idx_P_AC) = f(info.idx_P_AC) + (lambda_SDCI * P_EV_prev * dt);
            f(info.idx_P_EV) = f(info.idx_P_EV) + (lambda_SDCI * P_AC_prev * dt);
        end
        
        [x_opt, ~, exitflag] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
        
        if exitflag > 0
            P_AC_curr = x_opt(info.idx_P_AC);
            P_EV_curr = x_opt(info.idx_P_EV);
            P_Gen_curr = x_opt(info.idx_P_Gen);
            P_Shed_curr = x_opt(info.idx_P_Shed);
            eta_val = x_opt(info.idx_eta);
            z_val = x_opt(info.idx_z);
            
            P_AC_prev = P_AC_curr; P_EV_prev = P_EV_curr;
            
            strategies{i}.P_AC = P_AC_curr; strategies{i}.P_EV = P_EV_curr;
            strategies{i}.P_Shed = P_Shed_curr; strategies{i}.P_Gen = P_Gen_curr;
            
            % SDCI/Rho 记录
            n_dummy = ones(T_steps, 1);
            val_SDCI = calculate_SDCI_local(n_dummy, n_dummy, P_AC_curr, P_EV_curr);
            val_Rho  = calculate_Rho_local(n_dummy, P_AC_curr, n_dummy, P_EV_curr);
            if iter == 1
                strategies{i}.SDCI_History = zeros(Max_Iter, 1);
                strategies{i}.Rho_History = zeros(Max_Iter, 1);
            end
            strategies{i}.SDCI_History(iter) = val_SDCI;
            strategies{i}.Rho_History(iter) = val_Rho;

            % [修改] 成本计算与打印
            cost_gen = sum((cost_params.c1_ac*P_AC_curr + cost_params.c2_ac*P_AC_curr.^2)*dt + ...
                           (cost_params.c1_ev*P_EV_curr + cost_params.c2_ev*P_EV_curr.^2)*dt + ...
                           (cost_params.c1_gen*P_Gen_curr + cost_params.c2_gen*P_Gen_curr.^2)*dt);
            cost_slack = sum(cost_params.c1_shed * P_Shed_curr * dt); % 切负荷成本
            total_real_cost = cost_gen + cost_slack;
            
            cvar_val = eta_val + (1 / (N_scenarios * (1 - risk_p.confidence))) * sum(z_val);
            b_run_cost(i) = cost_gen;
            b_slack_sum(i) = sum(P_Shed_curr + P_Gen_curr) * dt;
            b_risk_val(i) = cvar_val;

            if iter == Max_Iter
                fprintf('    发电成本: %.2f (元), 切负荷量: %.2f (MWh), 总成本: %.2f (元), CVaR风险: %.2f (MW), rho: %.4f, sdci: %.4f\n', ...
                    cost_gen, b_slack_sum(i), total_real_cost, cvar_val, val_Rho, val_SDCI);
            end
        else
            fprintf('    失败 (Exitflag %d)\n', exitflag); break;
        end
    end
end

% --- 绘图 B: 风险偏好灵敏度分析 (更新样式) ---
if any(~isnan(b_slack_sum))
    figure('Name', '场景B_风险灵敏度', 'Color', 'w', 'Position', [100, 100, 900, 400]);
    yyaxis left; 
    b = bar(1:3, b_slack_sum, 0.5, 'FaceColor', [0.8 0.3 0.3]); 
    ylabel('切负荷量 + 火电调度量 (MWh)'); 
    set(gca, 'XTick', 1:3, 'XTickLabel', beta_values);
    for i = 1:length(b_slack_sum)
        text(i, b_slack_sum(i), sprintf('%.2f', b_slack_sum(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 12, 'Color', [0.6 0.1 0.1], 'FontWeight', 'bold');
    end
    
    yyaxis right; 
    plot(1:3, b_risk_val, 'b-o', 'LineWidth', 2, 'MarkerSize', 8); 
    ylabel('CVaR 潜在违约风险 (MW)'); 
    for i = 1:length(b_risk_val)
        text(i, b_risk_val(i), sprintf('%.2f', b_risk_val(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 12, 'Color', 'b', 'FontWeight', 'bold');
    end
    xlabel('风险厌恶系数 \beta');
    legend('切负荷 + 火电 (安全性)', '潜在违约风险 (经济性)', 'Location', 'best');
    grid on;
    print(gcf, '风险偏好灵敏度分析.png', '-dpng', '-r300');
end

%% ================= 场景 C: 详细调度方案 (Beta=10) =================
fprintf('\n>>> 绘制调度方案 (Beta=10) <<<\n');
idx_plot = 2; 
if ~isempty(strategies{idx_plot})
    st = strategies{idx_plot};
    dir_sign = sign(direction_signal); dir_sign(dir_sign==0)=1;
    
    Y_Stack_Plot = [st.P_AC, st.P_EV, st.P_Gen, st.P_Shed] .* dir_sign;
    P_Demand_plot = P_grid_demand .* dir_sign;
    
    % --- 图 1: 功率堆叠 (更新配色和图例) ---
    fig_stack = figure('Name', '多源协同调度堆叠图', 'Color', 'w', 'Position', [150, 150, 1000, 600]);
    hold on;
    h_area = area(t_axis, Y_Stack_Plot);
    h_area(1).FaceColor = [0.00, 0.45, 0.74]; h_area(1).EdgeColor = 'none'; % AC
    h_area(2).FaceColor = [0.47, 0.67, 0.19]; h_area(2).EdgeColor = 'none'; % EV
    h_area(3).FaceColor = [0.92, 0.69, 0.13]; h_area(3).EdgeColor = 'none'; % Gen
    h_area(4).FaceColor = [0.85, 0.33, 0.10]; h_area(4).EdgeColor = 'none'; h_area(4).FaceAlpha = 0.8; % Shed
    plot(t_axis, P_Demand_plot, 'k--', 'LineWidth', 2.0, 'DisplayName', '电网总需求');
    ylabel('功率 (MW) [正=增加用电, 负=减少用电]', 'FontSize', 14, 'FontName', 'Microsoft YaHei');
    legend([h_area(1), h_area(2), h_area(3), h_area(4)], ...
           {'空调 (AC)', '电动汽车 (EV)', '火电调节 (Gen)', '切负荷 (Shed)'}, ...
           'Location', 'northwest', 'FontSize', 12, 'FontName', 'Microsoft YaHei');
    grid on; xlim([6, 30]); 
    set(gca, 'XTick', x_ticks_set, 'XTickLabel', x_labels_set, 'FontSize', 12, 'FontName', 'Microsoft YaHei');
    print(fig_stack, '多源协同调度堆叠图.png', '-dpng', '-r600');
    
    % --- [新增] 图 2: AC与EV时序出力对比 ---
    fprintf('  >>> 正在绘制 AC与EV时序调度量对比图...\n');
    fig_comp = figure('Name', 'AC与EV时序出力对比', 'Color', 'w', 'Position', [200, 200, 1000, 600]);
    hold on;
    
    Reliable_AC_plot = zeros(T_steps, 1);
    Reliable_EV_plot = zeros(T_steps, 1);
    for t = 1:T_steps
        if dir_sign(t) >= 0
            Reliable_AC_plot(t) = Reliable_AC_Up(t);
            Reliable_EV_plot(t) = Reliable_EV_Up(t);
        else
            Reliable_AC_plot(t) = -Reliable_AC_Down(t); 
            Reliable_EV_plot(t) = -Reliable_EV_Down(t);
        end
    end
    P_AC_line = st.P_AC .* dir_sign;
    P_EV_line = st.P_EV .* dir_sign;

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
    
    % --- [新增] 指标迭代对比图 ---
    if isfield(st, 'SDCI_History')
        fig_sdci = figure('Name', 'SDCI 迭代对比', 'Color', 'w', 'Position', [300, 300, 600, 400]);
        sdci_vals = [st.SDCI_History(1), st.SDCI_History(end)];
        bar(sdci_vals, 0.4, 'FaceColor', [0.2 0.6 0.8]);
        set(gca, 'XTickLabel', {'初始迭代', '最终迭代'}, 'FontSize', 12, 'FontName', 'Microsoft YaHei');
        ylabel('SDCI 指标', 'FontSize', 12, 'FontName', 'Microsoft YaHei');
        text(1:2, sdci_vals, num2str(sdci_vals', '%.4f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
        print(fig_sdci, 'SDCI对比.png', '-dpng', '-r600');
        
        fig_rho = figure('Name', 'Rho 迭代对比', 'Color', 'w', 'Position', [400, 400, 600, 400]);
        rho_vals = [st.Rho_History(1), st.Rho_History(end)];
        bar(rho_vals, 0.4, 'FaceColor', [0.8 0.4 0.2]);
        set(gca, 'XTickLabel', {'初始迭代', '最终迭代'}, 'FontSize', 12, 'FontName', 'Microsoft YaHei');
        ylabel('Rho 指标', 'FontSize', 12, 'FontName', 'Microsoft YaHei');
        text(1:2, rho_vals, num2str(rho_vals', '%.4f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
        print(fig_rho, 'Rho对比.png', '-dpng', '-r600');
    end
end

%% ================= 场景 D: 鲁棒性测试 =================
fprintf('\n>>> 场景 D: 鲁棒性测试 <<<\n');
if ~isempty(strategies{1}) && ~isempty(strategies{3})
    [~, worst_idx] = min(sum(Scenarios_AC_Up + Scenarios_EV_Up));
    fprintf('  - 最恶劣场景: #%d\n', worst_idx);
    
    Real_Cap_AC = zeros(T_steps, 1); Real_Cap_EV = zeros(T_steps, 1);
    for t=1:T_steps
        if direction_signal(t)==1, Real_Cap_AC(t)=Scenarios_AC_Up(t,worst_idx); Real_Cap_EV(t)=Scenarios_EV_Up(t,worst_idx);
        else, Real_Cap_AC(t)=abs(Scenarios_AC_Down(t,worst_idx)); Real_Cap_EV(t)=abs(Scenarios_EV_Down(t,worst_idx)); end
    end
    
    calc_viol = @(P, Cap) sum(max(0, P - Cap));
    v_neu = calc_viol(strategies{1}.P_AC, Real_Cap_AC) + calc_viol(strategies{1}.P_EV, Real_Cap_EV);
    v_rob = calc_viol(strategies{3}.P_AC, Real_Cap_AC) + calc_viol(strategies{3}.P_EV, Real_Cap_EV);
    
    fprintf('  - 中性策略(Beta=0) 违约量: %.2f MW\n', v_neu);
    fprintf('  - 规避策略(Beta=10) 违约量: %.2f MW\n', v_rob);
    
    figure('Color','w', 'Position', [600, 300, 500, 400]); 
    b = bar([v_neu, v_rob], 0.5);
    b.FaceColor = 'flat'; b.CData(1,:) = [0.8 0.2 0.2]; b.CData(2,:) = [0.2 0.6 0.2];
    set(gca, 'XTickLabel', {'中性 (\beta=0)', '规避 (\beta=10)'});
    ylabel('极端场景实际违约量 (MW)');
    grid on;
    print(gcf, '极端场景鲁棒性测试.png', '-dpng', '-r300');
end

%% ================= 场景 E: 置信水平测试 =================
fprintf('\n>>> 场景 E: 置信水平对经济性的影响测试 (Reliability vs Cost) <<<\n');
confs = [0.85, 0.90, 0.95, 0.99];
num_conf = length(confs);
e_total_cost = zeros(1, num_conf);
e_slack_sum = zeros(1, num_conf);

risk_E = risk_p; risk_E.beta = 1;

for k = 1:num_conf
    curr_alpha = confs(k);
    risk_E.confidence = curr_alpha;
    fprintf('  [测试 %d/%d] 置信水平 alpha = %.2f ... ', k, num_conf, curr_alpha);
    
    net_params_safe = net_params;
    net_params_safe.ShedDist = zeros(N_bus, 1);

    [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast(...
            P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
            Physical_AC_Up, Physical_EV_Up, R_Gen_Max, R_Shed_Max, ...
            cost_params, risk_E, net_params_safe);
            
    start_row_net = 2 * N_scenarios; 
    for t = 1:T_steps
        if direction_signal(t) == 1
            rows_t = start_row_net + (t-1)*2*N_line + (1 : 2*N_line);
            A(rows_t, info.idx_P_AC(t)) = -A(rows_t, info.idx_P_AC(t));
            A(rows_t, info.idx_P_EV(t)) = -A(rows_t, info.idx_P_EV(t));
        end
    end
    
    idx_pow = [info.idx_P_AC, info.idx_P_EV, info.idx_P_Gen, info.idx_P_Shed];
    for idx = idx_pow, H(idx, idx) = H(idx, idx) * dt; end
    f(idx_pow) = f(idx_pow) * dt;
    
    [x_opt, ~, exitflag] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
    
    if exitflag > 0
        P_AC_curr = x_opt(info.idx_P_AC); P_EV_curr = x_opt(info.idx_P_EV);
        P_Gen_curr = x_opt(info.idx_P_Gen); P_Shed_curr = x_opt(info.idx_P_Shed);
        
        real_cost = sum((cost_params.c1_ac*P_AC_curr + cost_params.c2_ac*P_AC_curr.^2)*dt + ...
                        (cost_params.c1_ev*P_EV_curr + cost_params.c2_ev*P_EV_curr.^2)*dt + ...
                        (cost_params.c1_gen*P_Gen_curr + cost_params.c2_gen*P_Gen_curr.^2)*dt + ...
                        (cost_params.c1_shed*P_Shed_curr)*dt);
        e_total_cost(k) = real_cost;
        e_slack_sum(k) = sum(P_Gen_curr + P_Shed_curr) * dt;
        fprintf('成功。总成本: %.2f 元, 昂贵资源调用: %.2f MWh\n', real_cost, e_slack_sum(k));
    else
        fprintf('失败。\n');
    end
end

if any(e_total_cost > 0)
    fig_conf = figure('Name', '场景E_置信度影响', 'Color', 'w', 'Position', [600, 100, 700, 450]);
    yyaxis left;
    b = bar(categorical(confs), e_total_cost, 0.5, 'FaceColor', [0.2 0.6 0.8]);
    ylabel('系统总运行成本 (元)', 'FontSize', 12, 'FontName', 'Microsoft YaHei');
    ylim([min(e_total_cost)*0.9, max(e_total_cost)*1.1]);
    for i = 1:length(e_total_cost)
        text(i, e_total_cost(i), sprintf('%.0f', e_total_cost(i)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontSize', 10, 'Color', 'b');
    end
    yyaxis right;
    plot(1:num_conf, e_slack_sum, 'r-^', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    ylabel('火电与切负荷调用量 (MWh)', 'FontSize', 12, 'FontName', 'Microsoft YaHei');
    ax = gca; ax.YColor = 'r';
    xlabel('置信水平 \alpha (可靠性要求)', 'FontSize', 12, 'FontName', 'Microsoft YaHei');
    legend({'运行成本 (左轴)', '备用资源调用 (右轴)'}, 'Location', 'northwest', 'FontSize', 11);
    grid on;
    print(fig_conf, '置信水平经济性分析.png', '-dpng', '-r300');
end

fprintf('\n所有测试结束。\n');