%% 
% 功能：IEEE 30 节点系统 VPP 调节能力理论验证 (模块化重构版)
% 
% 重构内容：保留所有初始化与网络构建逻辑，将场景 B, C, D, E 封装为独立函数。

clear; close all; clc;

%% ================= 1. 全局初始化与数据清洗 =================
fprintf('正在加载场景数据...\n');
data_file = 'reliable_regulation_domain_1000_mix_pbase_8am.mat';
if ~exist(data_file, 'file')
    error('数据文件缺失！请先运行 main_scenario_generation_diff_mix.m');
end
load(data_file);

% --- 数据清洗：去除 NaN 和 Inf ---
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
    t_axis = 8 : dt_sim : (8 + (T_steps-1)*dt_sim);
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

% mpc = case30();
mpc = case_ieee30;
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
Bus_KV = mpc.bus(:, 10);
Voltage_Weight = ones(N_bus, 1);
Voltage_Weight(Bus_KV < 100) = 0.2; 
Raw_Dist = Bus_Pd .* Voltage_Weight;
if sum(Raw_Dist) == 0
    Raw_Dist = ones(N_bus, 1); % 防止全0的保护
end
net_params.AcDist = Raw_Dist / sum(Raw_Dist);
net_params.EvDist = Raw_Dist / sum(Raw_Dist);
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
% Line_RateA(Line_RateA < 1e-3) = 130; 
% net_params.LineLimit = max(Line_RateA, Max_Base_Flow + 30); 
estimated_limits = max(Max_Base_Flow * 1.5, Max_Base_Flow + 20);
net_params.LineLimit = max(Line_RateA, estimated_limits);

R_Gen_Max = max(0, sum(Gen_Pmax) - Total_Sys_Load);
R_Shed_Max = 1e6 * ones(T_steps, 1); 

fprintf('物理约束计算完成：线路最大背景潮流 %.2f MW, 修正后最小限额 %.2f MW\n', ...
    max(Max_Base_Flow), min(net_params.LineLimit));
beta_values = [0, 1, 10];
% 准备通用绘图参数
x_ticks_set = [8, 14, 20, 26, 32];
x_labels_set = {'08:00', '14:00', '20:00', '02:00 (次日)', '08:00 (次日)'};
options = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'interior-point-convex', 'TolFun', 1e-6);

%% ================= 4. 执行各个场景 (函数调用) =================

%% 1. 运行场景 B: 风险偏好灵敏度分析
% 返回计算出的策略 (strategies) 供后续场景使用
strategies = run_scenario_B(beta_values, Max_Iter, N_scenarios, N_bus, N_line, dt, ...
    P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, Physical_AC_Up, Physical_EV_Up, ...
    R_Gen_Max, R_Shed_Max, cost_params, net_params, direction_signal, lambda_SDCI, lambda_Rho, options);

%% 2. 运行场景 C: 详细调度方案 (Beta=10)
run_scenario_C(strategies, t_axis, P_grid_demand, direction_signal, ...
    Reliable_AC_Up, Reliable_EV_Up, Reliable_AC_Down, Reliable_EV_Down, ...
    x_ticks_set, x_labels_set);

%% 3. 运行场景 D: 鲁棒性测试
run_scenario_D_ramp(strategies, Scenarios_AC_Up, Scenarios_EV_Up, Scenarios_AC_Down, Scenarios_EV_Down, direction_signal);

%% 4. 运行场景 E: 置信水平测试
run_scenario_E(P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
    Physical_AC_Up, Physical_EV_Up, R_Gen_Max, R_Shed_Max, ...
    cost_params, net_params, direction_signal, dt, options, N_scenarios, N_line, N_bus);

fprintf('\n所有测试结束。\n');
