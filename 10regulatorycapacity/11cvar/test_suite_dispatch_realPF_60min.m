clear; close all; clc;

%% ================= 1. 全局初始化 =================
fprintf('正在加载场景数据...\n');
data_file = 'reliable_regulation_domain_1000_mix_pbase.mat';
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

% --- [新增] 加载并处理基线功率 ---
if ~exist('Reliable_AC_Base', 'var') || ~exist('Reliable_EV_Base', 'var')
    % 如果 mat 文件是旧版，没有基线数据，给予警告或报错
    error('数据文件中缺失 AC/EV 基线功率数据 (Reliable_AC_Base)。请重新生成数据。');
end
Reliable_AC_Base = Reliable_AC_Base(:) * unit_scale; % T x 1 (MW)
Reliable_EV_Base = Reliable_EV_Base(:) * unit_scale; % T x 1 (MW)

% --- [修改] 成本参数 (符合中国电力市场现状) ---
% 说明：成本计算公式将改为 sum( (c1*P + c2*P^2) * dt )
% 单位：元/MWh
cost_params.c1_ac = 500;      % 500 元/MWh (0.5元/kWh) - 空调调节成本
cost_params.c2_ac = 50;       % 50 元/(MW)^2h - 二次项系数大幅降低，避免成本失真
cost_params.c1_ev = 400;      % 400 元/MWh (0.4元/kWh) - EV调节成本(略低于AC)
cost_params.c2_ev = 50;       % 50 元/(MW)^2h
cost_params.c_slack = 30000;  % 30,000 元/MWh (30元/kWh) - 切负荷惩罚 (VoLL)
% 火电与切负荷成本 (用于 construct_risk_constrained_qp_fast 内部)
cost_params.c1_gen = 800;      cost_params.c2_gen = 80;     % 火电调节成本
cost_params.c1_shed = 30000;   cost_params.c2_shed = 0;     % 切负荷惩罚

% --- 互补性与相关性优化权重 ---
lambda_SDCI = 10;   
lambda_Rho  = 10;   
Max_Iter    = 3;    

% --- 修正2: 构造混合需求 (修改为每1小时生成一次，符合国网规律) ---
rng(105); 
P_grid_demand = zeros(T_steps, 1);
Effective_Scen_AC = zeros(T_steps, N_scenarios);
Effective_Scen_EV = zeros(T_steps, N_scenarios);
Effective_Phys_AC = zeros(T_steps, 1);
Effective_Phys_EV = zeros(T_steps, 1);
Effective_Reliable_AC = zeros(T_steps, 1); 
Effective_Reliable_EV = zeros(T_steps, 1); 

% [修改开始]：设置每 1 小时变化一次
steps_per_hour_update = round(1 / dt); % 1小时对应的仿真步数 (例如 12步)

direction_signal = zeros(T_steps, 1); % 用于存储展开后的方向信号 (供后续绘图)
current_block_demand = 0;             % 当前块的需求值缓存

fprintf('生成混合需求: 指令每 1 小时更新一次 (符合国网峰谷规律)。\n');
fprintf('  - 上调 (Up, >0): 增加用电 (填谷)\n');
fprintf('  - 下调 (Down, <0): 减少用电 (削峰)\n');

for t = 1:T_steps
    % 1. 获取当前时刻 (绝对时间，例如 6.5, 14.0)
    current_time_abs = t_axis(t);
    current_hour_of_day = mod(floor(current_time_abs), 24); % 获取当前小时数 (0-23)
    
    % 2. 判断调节方向 (State Grid Scheduling Pattern)
    % 用户定义：上调=增加用电，下调=减少用电
    
    % 高峰时段 (Peak): 08:00-12:00, 14:00-22:00
    % 电网需求：减少负荷 (削峰) -> 调用用户的 [下调潜力] (Down)
    if (current_hour_of_day >= 8 && current_hour_of_day < 12) || ...
       (current_hour_of_day >= 14 && current_hour_of_day < 22)
        use_up_potential = false; 
        direction_signal(t) = -1; % 负号表示减少用电 (下调)
        
    % 低谷时段 (Valley): 22:00-08:00, 12:00-14:00
    % 电网需求：增加负荷 (填谷) -> 调用用户的 [上调潜力] (Up)
    else
        use_up_potential = true;
        direction_signal(t) = 1;  % 正号表示增加用电 (上调)
    end
    
    % 3. 根据方向映射有效数据
    if use_up_potential
        % 使用上调数据 (增加用电)
        Effective_Scen_AC(t,:) = Scenarios_AC_Up(t,:);
        Effective_Scen_EV(t,:) = Scenarios_EV_Up(t,:);
        Effective_Phys_AC(t) = Physical_AC_Up(t);
        Effective_Phys_EV(t) = Physical_EV_Up(t);
        Effective_Reliable_AC(t) = Reliable_AC_Up(t);
        Effective_Reliable_EV(t) = Reliable_EV_Up(t);
    else
        % 使用下调数据 (减少用电)，取绝对值计算
        Effective_Scen_AC(t,:) = abs(Scenarios_AC_Down(t,:));
        Effective_Scen_EV(t,:) = abs(Scenarios_EV_Down(t,:));
        Effective_Phys_AC(t) = Physical_AC_Down(t);
        Effective_Phys_EV(t) = Physical_EV_Down(t);
        Effective_Reliable_AC(t) = Reliable_AC_Down(t);
        Effective_Reliable_EV(t) = Reliable_EV_Down(t);
    end
    
    % 4. 计算并锁定每小时的需求幅度
    % 逻辑：如果是小时的开始 (t=1, 13, 25...)，计算一个新的需求值并锁定
    if mod(t-1, steps_per_hour_update) == 0
        cap_rel = Effective_Reliable_AC(t) + Effective_Reliable_EV(t);
        cap_phy = Effective_Phys_AC(t) + Effective_Phys_EV(t);
        
        % 需求 = 可靠容量 + 50% * (风险区间 * 随机因子)
        current_block_demand = cap_rel + 0.50 * (cap_phy - cap_rel) * rand();
    end
    
    % 在1小时内保持需求恒定
    P_grid_demand(t) = current_block_demand;
end
% [修改结束]

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

% --- [关键部分]：物理背景潮流与实时容量计算 ---

% 0. 定义节点分布向量
net_params.AcDist = [0.1; 0.15; 0.3; 0.3; 0.15]; 
net_params.EvDist = [0.0; 0.2; 0.4; 0.4; 0.0];
net_params.GenDist = [0.25; 0; 0; 0; 0.75];     % 火电节点 1, 5
net_params.ShedDist = [0; 0.3; 0.3; 0.4; 0];    % 负荷节点 2, 3, 4

% 1. 基础刚性负荷 (Fixed Load)
Base_Load_Fixed = [0; 300; 300; 400; 0]; % (MW)

% 2. 生成日负荷曲线 (用于刚性负荷)
t_vec = linspace(0, 24, T_steps);
load_curve = 0.55 + 0.45 * exp(-((t_vec - 19).^2) / 12) + 0.1 * exp(-((t_vec - 10).^2) / 20);
load_curve = load_curve(:) / max(load_curve);

% 3. [新增] 循环计算每一时刻的背景潮流和系统总负荷
P_inj_t = zeros(N_bus, T_steps);       % 节点净注入矩阵
Total_System_Load_t = zeros(T_steps, 1); % 系统总实时负荷

for t = 1:T_steps
    % (1) 计算各节点上的实时负荷
    % 刚性负荷 (按曲线波动)
    P_Fixed_Node = Base_Load_Fixed * load_curve(t); 
    
    % AC/EV 基线负荷 (按分布向量映射)
    P_AC_Node = net_params.AcDist * Reliable_AC_Base(t);
    P_EV_Node = net_params.EvDist * Reliable_EV_Base(t);
    
    % 该时刻各节点总负荷
    P_Load_Node_Total = P_Fixed_Node + P_AC_Node + P_EV_Node;
    
    % (2) 记录系统总负荷 (用于计算切负荷上限和火电基线)
    Total_System_Load_t(t) = sum(P_Load_Node_Total);
    
    % (3) 计算各节点发电 (假设火电完全平衡负荷)
    % 总发电量 = 总负荷
    P_Gen_Node_Total = net_params.GenDist * Total_System_Load_t(t);
    
    % (4) 计算节点净注入 (Generation - Load)
    P_inj_t(:, t) = P_Gen_Node_Total - P_Load_Node_Total;
end

% 4. 计算包含 AC/EV 基线的背景潮流 (N_line x T)
net_params.BaseFlow = net_params.PTDF * P_inj_t;

% 5. 设置线路限额 (基于最大背景潮流)
max_base_flow = max(max(abs(net_params.BaseFlow)));
net_params.LineLimit = ones(N_line, 1) * (max_base_flow + 5.0);

% 6. 计算实时的调节容量物理上限

% 火电装机容量 (MW) - G1(350) + G5(900)
Gen_Capacity_Installed = 1000; 

% 火电调节能力(上调) = 装机容量 - 当前基线出力
% 当前基线出力 = Total_System_Load_t (因为 Gen = Load)
R_Gen_Max = Gen_Capacity_Installed - Total_System_Load_t;
R_Gen_Max = max(0, R_Gen_Max); % 物理保护

% 切负荷能力 = 当前总负荷
R_Shed_Max = Total_System_Load_t;

fprintf('物理约束计算完成：\n');
fprintf('  - 系统峰值总负荷 (含AC/EV): %.2f MW\n', max(Total_System_Load_t));
fprintf('  - 线路最大背景潮流: %.2f MW\n', max_base_flow);
fprintf('  - 峰值时刻火电备用: %.2f MW\n', min(R_Gen_Max));


%% ================= 场景 B: 风险偏好灵敏度分析 =================
fprintf('\n>>> 场景 B: 风险偏好灵敏度分析 (含火电与切负荷) <<<\n');
x_ticks_set = [6, 12, 18, 24, 30];
x_labels_set = {'06:00', '12:00', '18:00', '24:00', '06:00(+1)'};
beta_values = [0, 1, 10];
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
    risk_p.rho_pen = 500; 
    
    P_AC_prev = zeros(T_steps, 1);
    P_EV_prev = zeros(T_steps, 1);
    
    for iter = 1:Max_Iter
        % 1. 构建模型 (传入实时的 R_Gen_Max 和 R_Shed_Max)
        [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast(...
            P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
            Physical_AC_Up, Physical_EV_Up, ...
            R_Gen_Max, R_Shed_Max, ...  % <--- 这里传入的是向量
            cost_params, risk_p, net_params);
        
        % 时间积分修正
        idx_pow = [info.idx_P_AC, info.idx_P_EV, info.idx_P_Gen, info.idx_P_Shed];
        for idx = idx_pow, H(idx, idx) = H(idx, idx) * dt; end
        f(idx_pow) = f(idx_pow) * dt;
        
        % 2. 添加互补性惩罚
        if iter > 1
            penalty_vec_AC = lambda_SDCI * P_EV_prev * dt;
            penalty_vec_EV = lambda_SDCI * P_AC_prev * dt;
            P_EV_centered = P_EV_prev - mean(P_EV_prev);
            P_AC_centered = P_AC_prev - mean(P_AC_prev);
            penalty_rho_AC = lambda_Rho * P_EV_centered * dt;
            penalty_rho_EV = lambda_Rho * P_AC_centered * dt;
            
            f(info.idx_P_AC) = f(info.idx_P_AC) + penalty_vec_AC + penalty_rho_AC;
            f(info.idx_P_EV) = f(info.idx_P_EV) + penalty_vec_EV + penalty_rho_EV;
        end
        
        % 3. 求解
        [x_opt, fval, exitflag] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
        
        if exitflag > 0
            P_AC_curr = x_opt(info.idx_P_AC);
            P_EV_curr = x_opt(info.idx_P_EV);
            P_Gen_curr = x_opt(info.idx_P_Gen);
            P_Shed_curr = x_opt(info.idx_P_Shed);
            eta_val = x_opt(info.idx_eta);
            z_val = x_opt(info.idx_z);
            
            P_AC_prev = P_AC_curr; P_EV_prev = P_EV_curr;

            strategies{i}.P_AC = P_AC_curr;
            strategies{i}.P_EV = P_EV_curr;
            strategies{i}.P_Gen = P_Gen_curr;
            strategies{i}.P_Shed = P_Shed_curr;
            
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

            cost_real = sum((cost_params.c1_ac*P_AC_curr + cost_params.c1_ev*P_EV_curr + ...
                             cost_params.c1_gen*P_Gen_curr + cost_params.c1_shed*P_Shed_curr)*dt);
            cvar_val = eta_val + (1 / (N_scenarios * (1 - risk_p.confidence))) * sum(z_val);
            b_run_cost(i) = cost_real;
            b_risk_val(i) = cvar_val;
            
            % ==================== [修改点] ====================
            % 将原来的仅切负荷 改为 切负荷 + 火电调度
            b_slack_sum(i) = sum(P_Shed_curr + P_Gen_curr) * dt; 
            % ==================================================
            
            if iter == Max_Iter
                fprintf('    VPP: %.2f MWh, 火电: %.2f MWh, 切负荷: %.2f MWh, 总成本: %.2f 元, CVaR: %.2f\n', ...
                    sum(P_AC_curr+P_EV_curr)*dt, sum(P_Gen_curr)*dt, sum(P_Shed_curr)*dt, cost_real, cvar_val);
            end
        else
            fprintf('    求解失败 (Exitflag %d)\n', exitflag);
            break; 
        end
    end
end

% --- 绘图 B: 风险偏好灵敏度分析 ---
if any(~isnan(b_run_cost))
    figure('Name', '场景B_风险灵敏度', 'Color', 'w', 'Position', [100, 100, 900, 400]);
    yyaxis left; 
    b = bar(1:3, b_slack_sum, 0.5, 'FaceColor', [0.8 0.3 0.3]); 
    
    % ==================== [修改点] ====================
    % 修改 Y 轴标签
    ylabel('切负荷量 + 火电调度量 (MWh)'); 
    % ==================================================
    
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
    
    % ==================== [修改点] ====================
    % 修改图例名称
    legend('切负荷 + 火电 (安全性)', '潜在违约风险 (经济性)', 'Location', 'best');
    % ==================================================
    
    grid on;
    print(gcf, '风险偏好灵敏度分析.png', '-dpng', '-r300');
end

%% ================= [绘制详细调度方案 (Beta=10)] =================
fprintf('\n>>> 正在绘制调度方案 (Beta=10)... \n');
idx_plot = 2; 

if idx_plot <= length(strategies) && ~isempty(strategies{idx_plot})
    P_AC_opt = strategies{idx_plot}.P_AC;
    P_EV_opt = strategies{idx_plot}.P_EV;
    P_Gen = strategies{idx_plot}.P_Gen;
    P_Shed = strategies{idx_plot}.P_Shed;
    if isfield(strategies{idx_plot}, 'P_Slack')
        P_Slack_calc = strategies{idx_plot}.P_Slack;
    else
        P_Slack_calc = P_grid_demand - (P_AC_opt + P_EV_opt);
    end
    P_Slack_calc(P_Slack_calc < 1e-5) = 0; 
    
    direction_sign = sign(direction_signal);
    direction_sign(direction_sign == 0) = 1; 
    
    Y_Stack_Plot = [P_AC_opt, P_EV_opt, P_Gen, P_Shed] .* direction_sign;
    P_Demand_plot = P_grid_demand .* direction_sign;
    
    % --- 图 1: 功率堆叠 (独立窗口) ---
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
    
    % ================= [本次修改] 新增 图 2: AC与EV时序出力对比 =================
    fprintf('  >>> 正在绘制 AC与EV时序调度量对比图...\n');
    fig_comp = figure('Name', 'AC与EV时序出力对比', 'Color', 'w', 'Position', [200, 200, 1000, 600]);
    hold on;
    
    % 构造可靠边界曲线 (根据方向选择 Up 或 Down 边界)
    Reliable_AC_plot = zeros(T_steps, 1);
    Reliable_EV_plot = zeros(T_steps, 1);
    for t = 1:T_steps
        if direction_sign(t) >= 0
            Reliable_AC_plot(t) = Reliable_AC_Up(t);
            Reliable_EV_plot(t) = Reliable_EV_Up(t);
        else
            Reliable_AC_plot(t) = -Reliable_AC_Down(t); % 负值表示下调容量边界
            Reliable_EV_plot(t) = -Reliable_EV_Down(t);
        end
    end

    % 构造调度出力曲线 (带符号)
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
    grid on; 
    xlim([6, 30]); 
    set(gca, 'XTick', x_ticks_set, 'XTickLabel', x_labels_set, ...
        'FontSize', 12, 'FontName', 'Microsoft YaHei', 'LineWidth', 1.2);
    
    % 保存图 2
    print(fig_comp, 'AC与EV时序调度量对比_Optimized.png', '-dpng', '-r600');
    fprintf('  >>> 图2已保存: AC与EV时序调度量对比_Optimized.png\n');
    % =========================================================================
    
    % --- 指标迭代对比图 ---
    if isfield(strategies{idx_plot}, 'SDCI_History')
        fig_sdci = figure('Name', 'SDCI 迭代对比', 'Color', 'w', 'Position', [300, 300, 600, 400]);
        sdci_vals = [strategies{idx_plot}.SDCI_History(1), strategies{idx_plot}.SDCI_History(end)];
        bar(sdci_vals, 0.4, 'FaceColor', [0.2 0.6 0.8]);
        set(gca, 'XTickLabel', {'初始迭代', '最终迭代'}, 'FontSize', 12, 'FontName', 'Microsoft YaHei');
        ylabel('SDCI 指标', 'FontSize', 12, 'FontName', 'Microsoft YaHei');
        text(1:2, sdci_vals, num2str(sdci_vals', '%.4f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
        print(fig_sdci, 'SDCI对比.png', '-dpng', '-r600');
        
        fig_rho = figure('Name', 'Rho 迭代对比', 'Color', 'w', 'Position', [400, 400, 600, 400]);
        rho_vals = [strategies{idx_plot}.Rho_History(1), strategies{idx_plot}.Rho_History(end)];
        bar(rho_vals, 0.4, 'FaceColor', [0.8 0.4 0.2]);
        set(gca, 'XTickLabel', {'初始迭代', '最终迭代'}, 'FontSize', 12, 'FontName', 'Microsoft YaHei');
        ylabel('Rho 指标', 'FontSize', 12, 'FontName', 'Microsoft YaHei');
        text(1:2, rho_vals, num2str(rho_vals', '%.4f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
        print(fig_rho, 'Rho对比.png', '-dpng', '-r600');
    end
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
    
    fprintf('  - 中性策略(Beta=0) 违约量: %.2f MW\n', viol_neutral);
    fprintf('  - 规避策略(Beta=100) 违约量: %.2f MW\n', viol_robust);
    
    figure('Name', '场景D_鲁棒性', 'Color', 'w', 'Position', [600, 300, 500, 400]);
    b = bar([viol_neutral, viol_robust], 0.5);
    b.FaceColor = 'flat'; b.CData(1,:) = [0.8 0.2 0.2]; b.CData(2,:) = [0.2 0.6 0.2];
    set(gca, 'XTickLabel', {'中性 (\beta=0)', '规避 (\beta=100)'});
    ylabel('极端场景实际违约量 (MW)');
    grid on;
    print(gcf, '极端场景鲁棒性测试.png', '-dpng', '-r300');
end

fprintf('\n所有测试结束。\n');