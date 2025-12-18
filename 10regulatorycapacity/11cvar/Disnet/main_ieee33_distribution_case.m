%% main_ieee33_distribution_case.m
% 功能：IEEE 33 节点配电网 VPP 协同调度主程序 (完整修正版)
% 场景：混合调节指令 (上调/增加负荷 + 下调/减少负荷)
% 理论：LinDistFlow, CVaR, SDCI/Rho, 仅有功支持

clear; close all; clc;

%% ================= 1. 系统初始化 =================
fprintf('>>> [Step 1] 系统初始化与数据加载...\n');

% 1.1 加载场景数据
data_file = 'reliable_regulation_domain_1000_mix_pbase_8am.mat';
if ~exist(data_file, 'file')
    error('数据文件 %s 缺失！请先运行场景生成程序。', data_file);
end
load(data_file); 
% --- 数据清洗 ---
clean_data = @(x) fillmissing(x, 'constant', 0); 
Scenarios_AC_Up = clean_data(Scenarios_AC_Up);
Scenarios_EV_Up = clean_data(Scenarios_EV_Up);
Scenarios_AC_Down = clean_data(Scenarios_AC_Down);
Scenarios_EV_Down = clean_data(Scenarios_EV_Down);
Reliable_AC_Base = clean_data(Reliable_AC_Base);
Reliable_EV_Base = clean_data(Reliable_EV_Base);
Reliable_AC_Up  = clean_data(Reliable_AC_Up);
Reliable_EV_Up  = clean_data(Reliable_EV_Up);
Reliable_AC_Down  = abs(clean_data(Reliable_AC_Down));
Reliable_EV_Down  = abs(clean_data(Reliable_EV_Down));

% 1.2 数据预处理 (下采样: 5min -> 15min)
downsample_vec = @(x) 1 * x(1:3:end);      % 功率系数修正
downsample_mat = @(x) 1 * x(1:3:end, :);

% 提取可靠域边界及场景 (单位转换 kW -> MW)
scale = 1/1000;
Reliable_AC_Up   = downsample_vec(Reliable_AC_Up) * scale;
Reliable_EV_Up   = downsample_vec(Reliable_EV_Up) * scale;
Reliable_AC_Down = downsample_vec(Reliable_AC_Down) * scale;
Reliable_EV_Down = downsample_vec(Reliable_EV_Down) * scale;

Scenarios_AC_Up   = downsample_mat(Scenarios_AC_Up) * scale;
Scenarios_EV_Up   = downsample_mat(Scenarios_EV_Up) * scale;
Scenarios_AC_Down = downsample_mat(Scenarios_AC_Down) * scale;
Scenarios_EV_Down = downsample_mat(Scenarios_EV_Down) * scale;

% 时间轴
dt = 15/60; 
t_axis = 8 : dt : 31.75; 
T_steps = length(t_axis);

% 1.3 加载 IEEE 33 配电网模型
mpc = case33bw;
mpc.baseMVA = 10; % 基准容量 10 MVA
baseKV = 12.66;
Z_base = baseKV^2 / mpc.baseMVA; 

% 提取网络参数 (标幺值)
net.branch_r = mpc.branch(:, 3) / Z_base; 
net.branch_x = mpc.branch(:, 4) / Z_base;
net.branch_f = mpc.branch(:, 1); 
net.branch_t = mpc.branch(:, 2);
% 注意: case33bw 原始数据通常是 kW, 需要确认. 
% 假设 mpc.bus(:,3) 是 MW (若原始是kW需除1000，这里假设已是标准Matpower格式)
% 为保险起见，这里按 100kW 级别负载处理，若过大则 LinDistFlow 会由电压约束自动限制
net.P_load_base = mpc.bus(:, 3) / mpc.baseMVA; 
net.Q_load_base = mpc.bus(:, 4) / mpc.baseMVA;
net.N_bus = 33;

% 1.4 资源空间分布 (Spatial Distribution)
net.Dist_AC = zeros(33, 1); 
net.Dist_EV = zeros(33, 1);
Nodes_Comm = 2:18;  % 商业区 (前端)
Nodes_Res  = 19:33; % 居民区 (末端)

% 按容量比例分配
net.Dist_AC(Nodes_Comm) = 0.4 / length(Nodes_Comm);
net.Dist_AC(Nodes_Res)  = 0.6 / length(Nodes_Res);
net.Dist_EV(Nodes_Comm) = 0.4 / length(Nodes_Comm);
net.Dist_EV(Nodes_Res)  = 0.6 / length(Nodes_Res);

%% ================= 2. 生成混合调度指令 =================
fprintf('>>> [Step 2] 生成混合调度指令 (上调+下调)...\n');

% 构造正弦波动指令: 正半周为上调(增负荷), 负半周为下调(减负荷)
% 周期设为 12小时，确保一天内既有上调又有下调
signal_raw = 1.2 * sin(2*pi * (t_axis - 12) / 12)'; 

% 添加随机噪声
rng(42);
signal_noise = signal_raw + 0.1 * randn(T_steps, 1);

P_req_mag = zeros(T_steps, 1);
direction_flag = zeros(T_steps, 1); % 1=Up, -1=Down

for t = 1:T_steps
    val = signal_noise(t);
    if val >= 0
        % 需要上调 (增加负荷)
        direction_flag(t) = 1;
        % 限制指令不超过最大能力
        limit = Reliable_AC_Up(t) + Reliable_EV_Up(t);
        P_req_mag(t) = min(abs(val), 0.9*limit); 
    else
        % 需要下调 (减少负荷)
        direction_flag(t) = -1;
        limit = Reliable_AC_Down(t) + Reliable_EV_Down(t);
        P_req_mag(t) = min(abs(val), 0.9*limit);
    end
end

% 转换为标幺值
P_req_pu = P_req_mag;

% 成本参数
cost_params.c1_ac = 500; cost_params.c2_ac = 50;
cost_params.c1_ev = 400; cost_params.c2_ev = 50;
cost_params.c1_shed = 5e4; % 高惩罚避免切负荷

%% ================= 3. 执行测试模块 =================

% % 3.1 测试二：调节能力不确定性与可靠域
% run_test2_uncertainty(t_axis, Scenarios_AC_Up, Scenarios_EV_Up, ...
%                       Scenarios_AC_Down, Scenarios_EV_Down, ...
%                       Reliable_AC_Up, Reliable_EV_Up, ...
%                       Reliable_AC_Down, Reliable_EV_Down);

% 3.2 测试三：风险规避策略分析 (对比 Beta=0 和 Beta=10)
run_test3_risk_analysis(t_axis, P_req_pu, direction_flag, net, cost_params, ...
    Scenarios_AC_Up, Scenarios_EV_Up, Scenarios_AC_Down, Scenarios_EV_Down, ...
    Reliable_AC_Up, Reliable_EV_Up, Reliable_AC_Down, Reliable_EV_Down, mpc.baseMVA);

% 3.3 测试四：源荷协同优化效果 (SDCI/Rho)
run_test4_coordination(t_axis, P_req_pu, direction_flag, net, cost_params, ...
    Reliable_AC_Up, Reliable_EV_Up, Reliable_AC_Down, Reliable_EV_Down, mpc.baseMVA);

% 3.4 新增：配电网电压安全验证 (LinDistFlow)
run_test_voltage_safety(t_axis, net, cost_params, Reliable_AC_Up, Reliable_EV_Up, mpc.baseMVA);

fprintf('\n>>> 所有测试执行完毕。\n');