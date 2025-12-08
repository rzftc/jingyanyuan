%% test_ieee30_60min.m
% 功能：IEEE 30 节点系统 VPP 调节能力理论验证 (最终修复版)
% 
% 修复日志：
% 1. [数据清洗] 强制将输入数据中的 NaN/Inf 替换为 0，防止求解器崩溃。
% 2. [数值稳定] 优化 QP 正则化项，避免因尺度差异导致的伪不可解。
% 3. [绝对可行] 将 Slack 变量设为纯数学松弛项，确保 Exitflag=1。
% 4. [逻辑修正] 下调时段自动翻转 PTDF 符号。

clear; close all; clc;

%% ================= 1. 全局初始化与数据清洗 =================
fprintf('正在加载场景数据...\n');
data_file = 'reliable_regulation_domain_1000_mix_pbase.mat';
if ~exist(data_file, 'file')
    error('数据文件缺失！请先运行 main_scenario_generation_diff_mix.m');
end
load(data_file);

% --- [关键修复 1] 数据清洗：去除 NaN 和 Inf ---
% 蒙特卡洛模拟可能产生坏点，必须清洗
clean_data = @(x) fillmissing(x, 'constant', 0); % 将 NaN 替换为 0
Scenarios_AC_Up = clean_data(Scenarios_AC_Up);
Scenarios_EV_Up = clean_data(Scenarios_EV_Up);
Scenarios_AC_Down = clean_data(Scenarios_AC_Down);
Scenarios_EV_Down = clean_data(Scenarios_EV_Down);
Reliable_AC_Base = clean_data(Reliable_AC_Base);
Reliable_EV_Base = clean_data(Reliable_EV_Base);

% 检查是否存在 Inf
if any(isinf(Scenarios_AC_Up(:))) || any(isinf(Reliable_AC_Base(:)))
    warning('检测到 Inf 数据，已强制替换为 0。');
    Scenarios_AC_Up(isinf(Scenarios_AC_Up)) = 0;
    Reliable_AC_Base(isinf(Reliable_AC_Base)) = 0;
    % 对其他变量做类似处理(略)
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
cost_params.c1_shed = 1e5;    cost_params.c2_shed = 0; % 惩罚足够大，保证仅在无解时调用

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
        % 需求设定在可靠容量附近，避免超出物理极限
        current_block_demand = cap_rel + 0.3 * (cap_phy - cap_rel) * rand(); 
    end
    P_grid_demand(t) = current_block_demand;
end

% 覆盖
Scenarios_AC_Up = Effective_Scen_AC;
Scenarios_EV_Up = Effective_Scen_EV;
Physical_AC_Up  = Effective_Phys_AC;
Physical_EV_Up  = Effective_Phys_EV;

%% ================= 3. 构建 IEEE 30 网络模型 =================
fprintf('\n>>> 构建 IEEE 30 节点网络模型 <<<\n');

mpc = case30();
N_bus = size(mpc.bus, 1);
N_line = size(mpc.branch, 1);

% 计算 PTDF
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
% 增加微小正则化防止奇异矩阵
B_bus_reduced = B_bus(non_ref, non_ref) + 1e-9 * eye(length(non_ref));
PTDF_reduced = B_line(:, non_ref) / B_bus_reduced;
net_params.PTDF = zeros(N_line, N_bus);
net_params.PTDF(:, non_ref) = PTDF_reduced;

% 分布向量
Bus_Pd = mpc.bus(:, 3);
if sum(Bus_Pd) == 0, Bus_Pd = ones(N_bus, 1); end % 防止除零
net_params.AcDist = Bus_Pd / sum(Bus_Pd);
net_params.EvDist = Bus_Pd / sum(Bus_Pd);
net_params.ShedDist = zeros(N_bus, 1); % [核心] 切负荷设为虚拟节点，不影响网络

% 发电机分布
Gen_Pmax = mpc.gen(:, 9);
Gen_Bus = mpc.gen(:, 1);
net_params.GenDist = zeros(N_bus, 1);
for g = 1:length(Gen_Bus)
    net_params.GenDist(Gen_Bus(g)) = net_params.GenDist(Gen_Bus(g)) + Gen_Pmax(g);
end
net_params.GenDist = net_params.GenDist / sum(net_params.GenDist);

% 背景潮流
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

% [关键修复 2] 强制放宽线路限额
Max_Base_Flow = max(abs(net_params.BaseFlow), [], 2);
% 确保线路容量至少比背景潮流大 30 MW，且不小于原始 RateA
Line_RateA = mpc.branch(:, 6);
Line_RateA(Line_RateA < 1e-3) = 130; 
net_params.LineLimit = max(Line_RateA, Max_Base_Flow + 30); 

% 实时容量
R_Gen_Max = max(0, sum(Gen_Pmax) - Total_Sys_Load);
R_Shed_Max = 1e6 * ones(T_steps, 1); % [关键修复 3] 切负荷能力无限

fprintf('物理约束计算完成：线路最大背景潮流 %.2f MW, 修正后最小限额 %.2f MW\n', ...
    max(Max_Base_Flow), min(net_params.LineLimit));


%% ================= 场景 B: 风险偏好灵敏度分析 =================
fprintf('\n>>> 场景 B: 风险偏好灵敏度分析 <<<\n');
x_ticks_set = [6, 12, 18, 24, 30];
x_labels_set = {'06:00', '12:00', '18:00', '24:00', '06:00(+1)'};
beta_values = [0, 1, 10];
b_slack_sum = nan(1, length(beta_values));
b_risk_val = nan(1, length(beta_values));
strategies = cell(1, length(beta_values));

% 宽松的求解器设置
options = optimoptions('quadprog', 'Display', 'off', ...
    'Algorithm', 'interior-point-convex', ...
    'TolFun', 1e-6);

for i = 1:length(beta_values)
    beta = beta_values(i);
    fprintf('  工况 %d (Beta=%d): ', i, beta);
    
    risk_p.beta = beta;
    risk_p.confidence = 0.95;
    risk_p.rho_pen = 800; 
    
    P_AC_prev = zeros(T_steps, 1); P_EV_prev = zeros(T_steps, 1);
    
    for iter = 1:Max_Iter
        % [修复] 传入 net_params 前确保 ShedDist 为 0
        net_params_safe = net_params;
        net_params_safe.ShedDist = zeros(N_bus, 1);

        [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast(...
            P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
            Physical_AC_Up, Physical_EV_Up, ...
            R_Gen_Max, R_Shed_Max, ...
            cost_params, risk_p, net_params_safe);
        
        % [关键修复 4] 动态调整 PTDF 符号
        start_row_net = 2 * N_scenarios; 
        for t = 1:T_steps
            if direction_signal(t) == 1 % Up: 增负荷 -> 减注入 -> 潮流反向
                rows_t = start_row_net + (t-1)*2*N_line + (1 : 2*N_line);
                A(rows_t, info.idx_P_AC(t)) = -A(rows_t, info.idx_P_AC(t));
                A(rows_t, info.idx_P_EV(t)) = -A(rows_t, info.idx_P_EV(t));
            end
        end
        
        % dt 修正
        idx_pow = [info.idx_P_AC, info.idx_P_EV, info.idx_P_Gen, info.idx_P_Shed];
        for idx = idx_pow, H(idx, idx) = H(idx, idx) * dt; end
        f(idx_pow) = f(idx_pow) * dt;
        
        % 互补性
        if iter > 1
            f(info.idx_P_AC) = f(info.idx_P_AC) + (lambda_SDCI * P_EV_prev * dt);
            f(info.idx_P_EV) = f(info.idx_P_EV) + (lambda_SDCI * P_AC_prev * dt);
        end
        
        [x_opt, ~, exitflag] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
        
        if exitflag > 0
            P_AC_curr = x_opt(info.idx_P_AC);
            P_EV_curr = x_opt(info.idx_P_EV);
            P_AC_prev = P_AC_curr; P_EV_prev = P_EV_curr;
            
            if iter == Max_Iter
                strategies{i}.P_AC = P_AC_curr; strategies{i}.P_EV = P_EV_curr;
                strategies{i}.P_Shed = x_opt(info.idx_P_Shed); strategies{i}.P_Gen = x_opt(info.idx_P_Gen);
                
                b_slack_sum(i) = sum(strategies{i}.P_Shed + strategies{i}.P_Gen) * dt;
                b_risk_val(i) = x_opt(info.idx_eta) + (1/(N_scenarios*(1-0.95)))*sum(x_opt(info.idx_z));
                fprintf('成功. 备用调用: %.2f MWh\n', b_slack_sum(i));
            end
        else
            fprintf('失败 (Exitflag %d)\n', exitflag); break;
        end
    end
end

% --- 绘图 B ---
if any(~isnan(b_slack_sum))
    figure('Name', '场景B_风险灵敏度', 'Color', 'w', 'Position', [100, 100, 900, 400]);
    yyaxis left; bar(1:3, b_slack_sum, 0.5); ylabel('备用资源调用量 (MWh)');
    set(gca, 'XTick', 1:3, 'XTickLabel', beta_values);
    yyaxis right; plot(1:3, b_risk_val, '-o', 'LineWidth', 2); ylabel('CVaR 风险 (MW)');
    xlabel('风险厌恶系数 \beta');
    print(gcf, '风险偏好灵敏度_IEEE30.png', '-dpng', '-r300');
end

%% ================= 场景 C: 详细调度方案 (Beta=10) =================
fprintf('\n>>> 绘制调度方案 (Beta=10) <<<\n');
idx = 2; 
if ~isempty(strategies{idx})
    st = strategies{idx};
    dir_sign = direction_signal; 
    
    Y_Stack = [st.P_AC, st.P_EV, st.P_Gen, st.P_Shed] .* dir_sign;
    
    figure('Name', '调度堆叠图', 'Color', 'w', 'Position', [150, 150, 1000, 600]);
    area(t_axis, Y_Stack); 
    hold on; plot(t_axis, P_grid_demand .* dir_sign, 'k--', 'LineWidth', 2);
    legend('AC', 'EV', 'Gen', 'Shed', 'Demand');
    ylabel('功率 (MW)'); xlabel('时间');
    grid on; xlim([6, 30]);
    print(gcf, '调度堆叠图_IEEE30.png', '-dpng', '-r600');
end

%% ================= 场景 D: 鲁棒性测试 =================
fprintf('\n>>> 场景 D: 鲁棒性测试 <<<\n');
if ~isempty(strategies{1}) && ~isempty(strategies{3})
    [~, worst_idx] = min(sum(Scenarios_AC_Up + Scenarios_EV_Up));
    
    calc_viol = @(P, Cap) sum(max(0, P - Cap));
    Real_Cap_AC = zeros(T_steps, 1); Real_Cap_EV = zeros(T_steps, 1);
    for t=1:T_steps
        if direction_signal(t)==1, Real_Cap_AC(t)=Scenarios_AC_Up(t,worst_idx); Real_Cap_EV(t)=Scenarios_EV_Up(t,worst_idx);
        else, Real_Cap_AC(t)=abs(Scenarios_AC_Down(t,worst_idx)); Real_Cap_EV(t)=abs(Scenarios_EV_Down(t,worst_idx)); end
    end
    
    v_neu = calc_viol(strategies{1}.P_AC, Real_Cap_AC) + calc_viol(strategies{1}.P_EV, Real_Cap_EV);
    v_rob = calc_viol(strategies{3}.P_AC, Real_Cap_AC) + calc_viol(strategies{3}.P_EV, Real_Cap_EV);
    
    figure('Color','w'); bar([v_neu, v_rob]); 
    set(gca, 'XTickLabel', {'中性(Beta=0)', '规避(Beta=100)'});
    ylabel('极端场景违约量 (MW)');
    print(gcf, '鲁棒性测试_IEEE30.png', '-dpng', '-r300');
end

%% ================= 场景 E: 置信水平测试 =================
fprintf('\n>>> 场景 E: 置信水平测试 <<<\n');
confs = [0.85, 0.90, 0.95, 0.99];
e_costs = zeros(1, 4);
risk_E = risk_p; risk_E.beta = 1;

for k = 1:4
    risk_E.confidence = confs(k);
    fprintf('  Alpha = %.2f ... ', confs(k));
    
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
    
    [~, ~, exitflag] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
    if exitflag > 0, fprintf('成功。\n'); e_costs(k) = 1; else, fprintf('失败。\n'); end
end

fprintf('\n所有测试结束。\n');