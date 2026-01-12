%% 
% 功能：IEEE 30 节点系统 VPP 调节能力理论验证 (日前调度 15min 版)
% 
% 修正说明：
% 1. 步长设定：确认输入数据已为 15min 分辨率，不再进行下采样。
% 2. 规模扩展：将 2000AC/1000EV 数据线性扩展至 10000AC/5000EV (Scale=5)。
% 3. 聚合参数处理：保留 AC 灰盒参数，**移除 EV 灰盒参数处理** (不再使用期望SOC约束)。
% 4. 自动对齐：增加数据长度自动截断逻辑，防止索引越界。
% 5. [新增] 引入累积能量约束 (Reliable_EV_E_Up/Down) 并传递给优化模型。
% 6. [新增] 引入分时变化调度成本 (Time-Varying Cost)，体现火电/EV/AC在不同时段的经济性差异。

clear; close all; clc;

%% ================= 1. 全局初始化与数据清洗 =================
fprintf('正在加载场景数据...\n');
data_file = 'reliable_regulation_domain_soc_2bound.mat';
if ~exist(data_file, 'file')
    error('数据文件缺失！请先运行 main_scenario_generation_diff_mix.m');
end
load(data_file);

% --- 数据清洗 ---
clean_data = @(x) fillmissing(x, 'constant', 0); 
Scenarios_AC_Up = clean_data(Scenarios_AC_Up);
Scenarios_EV_Up = clean_data(Scenarios_EV_Up);
Scenarios_AC_Down = clean_data(Scenarios_AC_Down);
Scenarios_EV_Down = clean_data(Scenarios_EV_Down);
Reliable_AC_Base = clean_data(Reliable_AC_Base_95);
Reliable_EV_Base = clean_data(Reliable_EV_Base_95);
Reliable_AC_Up  = clean_data(Reliable_AC_Up);
Reliable_EV_Up  = clean_data(Reliable_EV_Up);
Reliable_AC_Down  = abs(clean_data(Reliable_AC_Down));
Reliable_EV_Down  = abs(clean_data(Reliable_EV_Down));
Reliable_EV_E_Up = quantile(Scenarios_EV_E_Up, 0.2, 2);
Reliable_EV_E_Down = quantile(Scenarios_EV_E_Down, 0.2, 2);
% [新增] 能量边界清洗
if exist('Reliable_EV_E_Up', 'var')
    Reliable_EV_E_Up = clean_data(Reliable_EV_E_Up);
    Reliable_EV_E_Down = abs(clean_data(Reliable_EV_E_Down));
else
    warning('未找到 Reliable_EV_E_Up/Down 变量，将使用无穷大边界(不约束能量)');
    Reliable_EV_E_Up = inf(size(Reliable_EV_Up));
    Reliable_EV_E_Down = inf(size(Reliable_EV_Up));
end

if any(isinf(Scenarios_AC_Up(:)))
    Scenarios_AC_Up(isinf(Scenarios_AC_Up)) = 0;
end

% ================= [关键] 数据对齐与截断 =================
% 目标长度：取 Reliable_EV_Params.A 的长度 (通常为 96)
if isfield(Reliable_EV_Params, 'A')
    target_len = length(Reliable_EV_Params.A);
    current_len = size(Scenarios_AC_Up, 1);
    
    if current_len > target_len
        fprintf('>>> 自动对齐: 将场景数据长度从 %d 截断为 %d 以匹配参数长度 <<<\n', current_len, target_len);
        
        % 1. 截断场景功率数据
        Scenarios_AC_Up = Scenarios_AC_Up(1:target_len, :);
        Scenarios_EV_Up = Scenarios_EV_Up(1:target_len, :);
        Scenarios_AC_Down = Scenarios_AC_Down(1:target_len, :);
        Scenarios_EV_Down = Scenarios_EV_Down(1:target_len, :);
        
        % 2. 截断可靠性边界与基线
        Reliable_AC_Base = Reliable_AC_Base(1:target_len);
        Reliable_EV_Base = Reliable_EV_Base(1:target_len);
        Reliable_AC_Up   = Reliable_AC_Up(1:target_len);
        Reliable_EV_Up   = Reliable_EV_Up(1:target_len);
        Reliable_AC_Down = Reliable_AC_Down(1:target_len);
        Reliable_EV_Down = Reliable_EV_Down(1:target_len);
        
        % [新增] 截断能量边界
        Reliable_EV_E_Up   = Reliable_EV_E_Up(1:target_len);
        Reliable_EV_E_Down = Reliable_EV_E_Down(1:target_len);
    end
end
% =======================================================

% --- [关键] 规模缩放与参数处理 ---
fprintf('正在处理聚合模型参数 (Scale x5, kW->MW)...\n');

% 规模缩放因子
Scale_AC = 5; 
Scale_EV = 10;

% 1. 应用缩放到功率数据
Scenarios_AC_Up = Scale_AC * Scenarios_AC_Up;
Scenarios_EV_Up = Scale_EV * Scenarios_EV_Up;
Scenarios_AC_Down = Scale_AC * Scenarios_AC_Down;
Scenarios_EV_Down = Scale_EV * Scenarios_EV_Down;

Reliable_AC_Base = Scale_AC * Reliable_AC_Base;
Reliable_EV_Base = Scale_EV * Reliable_EV_Base;

Reliable_AC_Up = Scale_AC * Reliable_AC_Up;
Reliable_EV_Up = Scale_EV * Reliable_EV_Up;
Reliable_AC_Down = Scale_AC * Reliable_AC_Down;
Reliable_EV_Down = Scale_EV * Reliable_EV_Down;

% [新增] 应用缩放到能量边界 (能量与功率线性相关，Scale相同)
Reliable_EV_E_Up = Scale_EV * Reliable_EV_E_Up;
Reliable_EV_E_Down = Scale_EV * Reliable_EV_E_Down;

% 2. 聚合参数处理 (仅处理 AC)
% AC 参数 (标量)
AC_Params_Scaled.A = Reliable_AC_Params.A;
AC_Params_Scaled.B = (Reliable_AC_Params.B / Scale_AC) * 1000; % kW -> MW, Scale corrected
AC_Params_Scaled.C = Reliable_AC_Params.C;

% [修改] 不再处理 EV 参数，因为我们放弃了 EV 的期望 SOC 约束
% 仅保留空结构体或必要信息防止报错 (如果下游需要)
EV_Params_Scaled = []; 

% 获取时间步数
[T_steps, N_scenarios] = size(Scenarios_AC_Up);

% --- 时间参数 ---
dt_sim = 15/60; % 0.25 小时
t_axis = 8 : dt_sim : (8 + (T_steps-1)*dt_sim);
dt = 15/60; 
fprintf('时间步长 dt = %.4f 小时 (15分钟)\n', dt);
fprintf('调度点数 T_steps = %d (日前96点)\n', T_steps);

% --- 单位转换 (kW -> MW) ---
unit_scale = 1/1000; 
Scenarios_AC_Up = Scenarios_AC_Up * unit_scale; 
Scenarios_EV_Up = Scenarios_EV_Up * unit_scale; 
Scenarios_AC_Down = Scenarios_AC_Down * unit_scale; 
Scenarios_EV_Down = Scenarios_EV_Down * unit_scale; 
Reliable_AC_Base = Reliable_AC_Base(:) * unit_scale; 
Reliable_EV_Base = Reliable_EV_Base(:) * unit_scale; 
Reliable_AC_Up = Reliable_AC_Up(:) * unit_scale;
Reliable_EV_Up = Reliable_EV_Up(:) * unit_scale;
Reliable_AC_Down = Reliable_AC_Down(:) * unit_scale;
Reliable_EV_Down = Reliable_EV_Down(:) * unit_scale;

% [新增] 能量边界单位转换 (kWh -> MWh)
Reliable_EV_E_Up = Reliable_EV_E_Up(:) * unit_scale;
Reliable_EV_E_Down = Reliable_EV_E_Down(:) * unit_scale;

Physical_AC_Up  = max(Scenarios_AC_Up, [], 2);
Physical_EV_Up  = max(Scenarios_EV_Up, [], 2);
Physical_AC_Down  = max(abs(Scenarios_AC_Down), [], 2);
Physical_EV_Down  = max(abs(Scenarios_EV_Down), [], 2);

% --- [关键修改] 成本参数 (分时变化) ---
% 1. 生成小时向量 (处理跨天逻辑: 8..23..0..8)
Hour_Vector = mod(t_axis, 24)'; 

% 2. 火电 (Gen): 模拟电力市场，早晚高峰高，深夜谷值低
C1_Gen_Vec = 800 * ones(T_steps, 1); % 基准
idx_gen_peak   = (Hour_Vector >= 10 & Hour_Vector < 12) | ...
                 (Hour_Vector >= 18 & Hour_Vector < 21);
idx_gen_valley = (Hour_Vector >= 23) | (Hour_Vector < 5); % 覆盖跨天午夜
C1_Gen_Vec(idx_gen_peak)   = 1000;
C1_Gen_Vec(idx_gen_valley) = 600;

% 3. 电动汽车 (EV): 模拟用户紧迫度，下班返家高，深夜充电低
C1_EV_Vec = 500 * ones(T_steps, 1); % 基准
idx_ev_urgent = (Hour_Vector >= 17 & Hour_Vector < 21);
idx_ev_low    = (Hour_Vector >= 23) | (Hour_Vector < 7); % 覆盖跨天睡眠期
C1_EV_Vec(idx_ev_urgent) = 800;
C1_EV_Vec(idx_ev_low)    = 300;

% 4. 空调 (AC): 模拟舒适度敏感性，午后高温高，夜间低
C1_AC_Vec = 400 * ones(T_steps, 1); % 基准
idx_ac_hot  = (Hour_Vector >= 13 & Hour_Vector < 16);
idx_ac_cool = (Hour_Vector >= 22) | (Hour_Vector < 8); % 覆盖跨天
C1_AC_Vec(idx_ac_hot)  = 600;
C1_AC_Vec(idx_ac_cool) = 300;

% 5. 赋值参数
cost_params.c1_ac = C1_AC_Vec;      cost_params.c2_ac = 10;     
cost_params.c1_ev = C1_EV_Vec;      cost_params.c2_ev = 10;       
cost_params.c1_gen = C1_Gen_Vec;    cost_params.c2_gen = 80; 
cost_params.c1_shed = 2e5;          cost_params.c2_shed = 0; 

fprintf('>>> 成本参数已更新为分时变化模式 (Time-Varying Costs) <<<\n');

% --- 优化权重 ---
lambda_SDCI = 10;   
lambda_Rho  = 10;   
Max_Iter    = 3;    

%% ================= 2. 准备 IEEE 30 基础数据 =================
fprintf('读取 IEEE 30 节点基础数据...\n');
mpc = case_ieee30;
N_bus = size(mpc.bus, 1);
N_line = size(mpc.branch, 1);

% 获取 IEEE 30 的固定负荷基准 (MW)
Bus_Pd_Base = mpc.bus(:, 3); 
Total_IEEE_Base_Load = sum(Bus_Pd_Base);
fprintf('  - IEEE 30 固定负荷基准: %.2f MW\n', Total_IEEE_Base_Load);

%% ================= 3. 物理环境建模：负荷、光伏与调节指令 =================
fprintf('构建源-荷物理模型与真实调节指令...\n');
rng(2024); 

P_grid_demand = zeros(T_steps, 1);
Effective_Scen_AC = zeros(T_steps, N_scenarios);
Effective_Scen_EV = zeros(T_steps, N_scenarios);
Effective_Phys_AC = zeros(T_steps, 1);
Effective_Phys_EV = zeros(T_steps, 1);
Effective_Reliable_AC = zeros(T_steps, 1); 
Effective_Reliable_EV = zeros(T_steps, 1); 
direction_signal = zeros(T_steps, 1); 

% --- 1. 修正：构建符合物理规律的日负荷曲线 ---
load_curve = 0.55 ...                                     % 基础底荷
           + 0.20 * exp(-((t_axis - 10).^2) / 15) ...     % 当日早高峰 (10:00)
           + 0.45 * exp(-((t_axis - 20).^2) / 12) ...     % 当日晚高峰 (20:00)
           + 0.20 * exp(-((t_axis - 34).^2) / 15);        % 次日早高峰影响

load_curve = load_curve(:) / max(load_curve); % 归一化

P_IEEE_Fixed_Load = Total_IEEE_Base_Load * load_curve .* (1 + 0.01 * randn(T_steps, 1));
P_System_Total_Load = P_IEEE_Fixed_Load + Reliable_AC_Base + Reliable_EV_Base;

% --- 2. 构建光伏出力曲线 (PV Profile) ---
PV_Capacity = max(P_System_Total_Load) * 0.4; 
P_PV_Total = zeros(T_steps, 1);
for i = 1:T_steps
    t_val = mod(t_axis(i), 24);
    if t_val > 6 && t_val < 18
        x = (t_val - 6) / 12;
        intensity = beta_pdf_proxy(x, 2.5, 2.5); 
        cloud_factor = 1 - 0.25 * rand() * (rand > 0.6); 
        P_PV_Total(i) = PV_Capacity * intensity * cloud_factor;
    else
        P_PV_Total(i) = 0;
    end
end

% --- 3. 计算“实时净负荷” ---
P_Net_Load = P_System_Total_Load - P_PV_Total;

% --- 4. 制定“日前发电计划” ---
P_Gen_Schedule = movmean(P_Net_Load, floor(4.5/dt)); 

% --- 5. 生成调节指令 ---
P_System_Mismatch = P_Net_Load - P_Gen_Schedule;

% --- 6. 分配指令并确定方向 ---
for t = 1:T_steps
    raw_val = P_System_Mismatch(t);
    
    if raw_val < 0
        direction_signal(t) = 1; % Up Regulation (Charging/Cooling, Load Increase)
        P_grid_demand(t) = abs(raw_val); 
        Effective_Scen_AC(t,:) = Scenarios_AC_Up(t,:);
        Effective_Scen_EV(t,:) = Scenarios_EV_Up(t,:);
        Effective_Phys_AC(t) = Physical_AC_Up(t);
        Effective_Phys_EV(t) = Physical_EV_Up(t);
        Effective_Reliable_AC(t) = Reliable_AC_Up(t);
        Effective_Reliable_EV(t) = Reliable_EV_Up(t);
    else
        direction_signal(t) = -1; % Down Regulation (Discharging/Heating, Load Decrease)
        P_grid_demand(t) = abs(raw_val);
        Effective_Scen_AC(t,:) = abs(Scenarios_AC_Down(t,:));
        Effective_Scen_EV(t,:) = abs(Scenarios_EV_Down(t,:));
        Effective_Phys_AC(t) = Physical_AC_Down(t);
        Effective_Phys_EV(t) = Physical_EV_Down(t);
        Effective_Reliable_AC(t) = Reliable_AC_Down(t);
        Effective_Reliable_EV(t) = Reliable_EV_Down(t);
    end
    
    if P_grid_demand(t) < 1e-4, P_grid_demand(t) = 1e-4; end
end

Scenarios_AC_Up = Effective_Scen_AC;
Scenarios_EV_Up = Effective_Scen_EV;
Physical_AC_Up  = Effective_Phys_AC;
Physical_EV_Up  = Effective_Phys_EV;

%% ================= 绘图验证 (论文格式修改) =================

% --- 图1：系统功率平衡 ---
fig1 = figure('Name', 'System_Balance', 'Color', 'w', 'Position', [100, 100, 800, 400]);
hold on;
plot(t_axis, P_System_Total_Load, 'k-', 'LineWidth', 1.5, 'DisplayName', '系统总负荷 (含AC/EV)');
plot(t_axis, P_Net_Load, 'g--', 'LineWidth', 1.5, 'DisplayName', '净负荷 (扣除光伏)');
plot(t_axis, P_Gen_Schedule, 'b-.', 'LineWidth', 1.5, 'DisplayName', '日前发电计划');
legend('Location', 'best'); 
ylabel('功率 (MW)'); xlabel('时刻'); 
grid on;

% 修改 X 轴：[8, 32] 对应 08:00 到次日 08:00
xlim([8, 32]);
set(gca, 'XTick', [8, 12, 16, 20, 24, 28, 32]);
set(gca, 'XTickLabel', {'08:00', '12:00', '16:00', '20:00', '00:00', '04:00', '08:00'});

% 保存为中文文件名
print(fig1, '源荷平衡验证.png', '-dpng', '-r600');
fprintf('  > 已保存: 源荷平衡验证.png (600 DPI)\n');

% --- 图2：调节指令 ---
fig2 = figure('Name', 'Regulation_Instruction', 'Color', 'w', 'Position', [100, 550, 800, 400]);
hold on;
area(t_axis, P_grid_demand .* direction_signal, 'FaceColor', [0.7 0.7 0.9]);
ylabel('调节指令 (MW)'); xlabel('时刻'); 
grid on;

% 修改 X 轴：[8, 32] 对应 08:00 到次日 08:00
xlim([8, 32]);
set(gca, 'XTick', [8, 12, 16, 20, 24, 28, 32]);
set(gca, 'XTickLabel', {'08:00', '12:00', '16:00', '20:00', '00:00', '04:00', '08:00'});

% 保存为中文文件名
print(fig2, '调节指令验证.png', '-dpng', '-r600');
fprintf('  > 已保存: 调节指令验证.png (600 DPI)\n');
fprintf('  - 光伏峰值: %.2f MW\n', max(P_PV_Total));
fprintf('  - 系统总负荷峰值: %.2f MW\n', max(P_System_Total_Load));
fprintf('  - 调节指令峰值 (绝对值): %.2f MW\n', max(P_grid_demand));
close all;

%% ================= 4. 构建 IEEE 30 节点网络模型 =================
fprintf('\n>>> 构建 IEEE 30 节点网络模型 (含光伏与VPP基线) <<<\n');

% 构建 PTDF
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

% 节点分布权重
if sum(Bus_Pd_Base) == 0, Bus_Pd_Base = ones(N_bus, 1); end 
net_params.AcDist = Bus_Pd_Base / sum(Bus_Pd_Base); 
net_params.EvDist = Bus_Pd_Base / sum(Bus_Pd_Base); 
net_params.ShedDist = zeros(N_bus, 1); 

Gen_Pmax = mpc.gen(:, 9);
Gen_Pmin = mpc.gen(:, 10); 
Gen_Bus = mpc.gen(:, 1);
net_params.GenDist = zeros(N_bus, 1);
for g = 1:length(Gen_Bus)
    net_params.GenDist(Gen_Bus(g)) = net_params.GenDist(Gen_Bus(g)) + Gen_Pmax(g);
end
net_params.GenDist = net_params.GenDist / sum(net_params.GenDist);

% 光伏分布
PV_Nodes = 15:30;
net_params.PVDist = zeros(N_bus, 1);
net_params.PVDist(PV_Nodes) = 1; 
net_params.PVDist = net_params.PVDist / sum(net_params.PVDist);

% --- 关键：计算背景潮流 (Base Flow) ---
P_inj_t = zeros(N_bus, T_steps);

for t = 1:T_steps
    P_Node_Fixed = net_params.AcDist * P_IEEE_Fixed_Load(t);
    P_Node_AC_Base = net_params.AcDist * Reliable_AC_Base(t);
    P_Node_EV_Base = net_params.EvDist * Reliable_EV_Base(t);
    P_Node_Load_Total = P_Node_Fixed + P_Node_AC_Base + P_Node_EV_Base;
    P_Node_PV = net_params.PVDist * P_PV_Total(t);
    P_Node_Gen = net_params.GenDist * P_Gen_Schedule(t);
    P_inj_t(:, t) = P_Node_Gen + P_Node_PV - P_Node_Load_Total;
end

net_params.BaseFlow = net_params.PTDF * P_inj_t;
Max_Base_Flow = max(abs(net_params.BaseFlow), [], 2);
Line_RateA = mpc.branch(:, 6);
net_params.LineLimit = max(Line_RateA, Max_Base_Flow + 20); 

% --- [新增] 将灰盒聚合参数注入 net_params ---
net_params.AC_Params = AC_Params_Scaled;
net_params.EV_Params = []; % [修改] 置空 EV 参数
net_params.Reliable_EV_Base = Reliable_EV_Base; 
net_params.Direction_Signal = direction_signal; 

% [新增] 将能量边界和时间步长注入 net_params
net_params.Reliable_EV_E_Up = Reliable_EV_E_Up;
net_params.Reliable_EV_E_Down = Reliable_EV_E_Down;
net_params.dt = dt;

% --- 修正火电调节能力 ---
R_Gen_Up_Cap   = max(0, sum(Gen_Pmax) - P_Gen_Schedule);   
R_Gen_Down_Cap = max(0, P_Gen_Schedule - sum(Gen_Pmin));   

R_Gen_Effective = zeros(T_steps, 1);
for t = 1:T_steps
    if direction_signal(t) == 1 
        R_Gen_Effective(t) = R_Gen_Down_Cap(t);
    else
        R_Gen_Effective(t) = R_Gen_Up_Cap(t);
    end
end
R_Gen_Max = R_Gen_Effective; 
R_Shed_Max = 1e6 * ones(T_steps, 1); 

fprintf('物理约束计算完成：线路最大背景潮流 %.2f MW\n', max(Max_Base_Flow));

% --- Beta 扫描范围 ---
beta_values = [0, 1, 10]; 
x_ticks_set = [8, 14, 20, 26, 32];
x_labels_set = {'08:00', '14:00', '20:00', '02:00 (次日)', '08:00 (次日)'};
options = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'interior-point-convex', 'TolFun', 1e-6);

%% ================= 5. 执行各个场景 =================

%% 1. 运行场景 B: 风险偏好灵敏度分析
strategies = run_scenario_B_tly(beta_values, Max_Iter, N_scenarios, N_bus, N_line, dt, ...
    P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, Effective_Reliable_AC, Effective_Reliable_EV, ...
    R_Gen_Max, R_Shed_Max, cost_params, net_params, direction_signal, lambda_SDCI, lambda_Rho, options);

%% 2. 运行场景 C: 详细调度方案 (Beta=中等)
run_scenario_C(strategies, t_axis, P_grid_demand, direction_signal, ...
    Reliable_AC_Up, Reliable_EV_Up, Reliable_AC_Down, Reliable_EV_Down, ...
    x_ticks_set, x_labels_set);

%% 3. 运行场景 D: 鲁棒性测试
run_scenario_D_ramp(strategies, Scenarios_AC_Up, Scenarios_EV_Up, Scenarios_AC_Down, Scenarios_EV_Down, direction_signal);

%% 4. 运行场景 E: 置信水平测试
run_scenario_E_tly(P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
    Physical_AC_Up, Physical_EV_Up, R_Gen_Max, R_Shed_Max, ...
    cost_params, net_params, direction_signal, dt, options, N_scenarios, N_line, N_bus);

fprintf('\n所有测试结束。\n');
%% 5. 运行场景 F: 协同约束效益对比验证
beta_for_comparison = 46; 

run_scenario_F_comparison(beta_for_comparison, 2, N_scenarios, N_bus, N_line, dt, ...
    P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
    Effective_Reliable_AC, Effective_Reliable_EV, ... % <--- 注意：这里使用可靠边界！
    R_Gen_Max, R_Shed_Max, cost_params, net_params, direction_signal, 10, 10, options);
%% 5. 运行场景 G: 确定性 vs 随机优化 效益对比分析
beta_for_G = 100; 

run_scenario_G_comparison(beta_for_G, Max_Iter, N_scenarios, N_bus, N_line, dt, ...
    P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
    Reliable_AC_Up, Reliable_EV_Up, ...
    Reliable_AC_Down, Reliable_EV_Down, ...
    Reliable_AC_Base, Reliable_EV_Base, ...
    R_Gen_Max, R_Shed_Max, cost_params, net_params, direction_signal, options);
%% 6. 运行场景 H: 传统鲁棒优化 vs CVaR随机优化 对比验证
% 设定一个用于对比的风险厌恶系数，通常与场景G或F保持一致
beta_for_H = 10; 

run_scenario_H_robust_comparison(beta_for_H, N_scenarios, N_bus, N_line, dt, ...
    P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
    R_Gen_Max, R_Shed_Max, cost_params, net_params, direction_signal, options);
%% --- 辅助函数：模拟钟形曲线 ---
function y = beta_pdf_proxy(x, a, b)
    if x < 0 || x > 1
        y = 0;
    else
        y = x^(a-1) * (1-x)^(b-1);
        y = y / 0.28; 
    end
end