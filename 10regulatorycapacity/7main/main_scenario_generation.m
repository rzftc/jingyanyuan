%% main_scenario_generation.m
% 实现大论文 4.2.1 节：基于蒙特卡洛模拟的用户行为不确定性场景生成
% 功能：循环调用修改后的仿真函数，生成 AC 和 EV 聚合体的调节能力场景集

clc; clear; close all;

%% 1. 参数设置
num_scenarios = 100;  % 设定场景生成数量 N (根据计算能力可调至 500 或 1000)
T_steps = 288;        % 24小时 + 1个点 (0到24小时，5分钟步长 -> 289个点)
                      % 注意：如果您的仿真逻辑产生 288 个点，请相应调整。
                      % 下方的逻辑会自动适配实际输出长度。

% 初始化存储矩阵 [时间步 x 场景数]
% 预分配时假设长度足够，之后会根据第一次运行结果截断或调整
Scenarios_AC_Up   = [];
Scenarios_AC_Down = [];
Scenarios_EV_Up   = [];
Scenarios_EV_Down = [];

fprintf('开始执行蒙特卡洛模拟，生成 %d 个场景...\n', num_scenarios);
tic; % 计时开始

%% 2. 蒙特卡洛模拟循环
% 使用 parfor 加速 (如果计算机支持并行工具箱)，否则改回 for
% 注意：parfor 内部无法直接动态改变 Scenarios 矩阵大小，需预先定义或使用 cell

temp_AC_Up = cell(1, num_scenarios);
temp_AC_Down = cell(1, num_scenarios);
temp_EV_Up = cell(1, num_scenarios);
temp_EV_Down = cell(1, num_scenarios);

parfor s = 1:num_scenarios
    
    % 打印进度
    if mod(s, 10) == 0
        fprintf('正在生成场景: %d / %d\n', s, num_scenarios);
    end
    
    % 生成随机种子，确保每个场景的用户行为参数不同
    % 种子控制：用户参与度、温度偏好偏差、EV入离网时间等不确定性
    current_seed = 2024 + s; 
    
    % --- 2.1 生成 AC 场景 ---
    % 调用修改后的 AC 仿真函数
    [ac_up, ac_down] = run_AC_simulation_MC(current_seed);
    
    % --- 2.2 生成 EV 场景 ---
    % 调用修改后的 EV 仿真函数
    [ev_up, ev_down] = run_EV_simulation_MC(current_seed);
    
    % 存储到临时 cell 中
    temp_AC_Up{s}   = ac_up;
    temp_AC_Down{s} = ac_down;
    temp_EV_Up{s}   = ev_up;
    temp_EV_Down{s} = ev_down;
end

toc; % 计时结束
fprintf('场景生成完成。正在整理数据...\n');

%% 3. 数据整理与保存
% 假设所有仿真的时间步长一致，取第一个场景确定维度
T_actual = length(temp_AC_Up{1});

Scenarios_AC_Up   = zeros(T_actual, num_scenarios);
Scenarios_AC_Down = zeros(T_actual, num_scenarios);
Scenarios_EV_Up   = zeros(T_actual, num_scenarios);
Scenarios_EV_Down = zeros(T_actual, num_scenarios);

for s = 1:num_scenarios
    % 简单的容错：如果某个场景长度不一致(极少见)，进行截断或补零
    len = length(temp_AC_Up{s});
    valid_len = min(len, T_actual);
    
    Scenarios_AC_Up(1:valid_len, s)   = temp_AC_Up{s}(1:valid_len);
    Scenarios_AC_Down(1:valid_len, s) = temp_AC_Down{s}(1:valid_len);
    
    len_ev = length(temp_EV_Up{s});
    valid_len_ev = min(len_ev, T_actual);
    Scenarios_EV_Up(1:valid_len_ev, s)   = temp_EV_Up{s}(1:valid_len_ev);
    Scenarios_EV_Down(1:valid_len_ev, s) = temp_EV_Down{s}(1:valid_len_ev);
end

% 保存生成的场景集
output_file = 'uncertainty_scenarios_dataset.mat';
save(output_file, 'Scenarios_AC_Up', 'Scenarios_AC_Down', ...
     'Scenarios_EV_Up', 'Scenarios_EV_Down', 'num_scenarios');
fprintf('场景集已保存至: %s\n', output_file);

%% 4. (可选) 简易可视化：绘制场景集云图
figure('Name', '不确定性场景集可视化', 'Position', [100, 100, 1200, 600], 'Color', 'w');
time_axis = linspace(6, 30, T_actual); % 6:00 到 次日 6:00

% 绘制 AC 上调能力
subplot(2, 2, 1); hold on;
plot(time_axis, Scenarios_AC_Up, 'Color', [0.6 0.8 1 0.1]); % 半透明淡蓝线
plot(time_axis, mean(Scenarios_AC_Up, 2), 'b-', 'LineWidth', 1.5); % 均值线
ylabel('AC 上调潜力 (kW)'); title('AC 场景集 (上调)'); xlim([6 30]);

% 绘制 AC 下调能力
subplot(2, 2, 2); hold on;
plot(time_axis, Scenarios_AC_Down, 'Color', [1 0.8 0.6 0.1]); % 半透明淡橙线
plot(time_axis, mean(Scenarios_AC_Down, 2), 'r-', 'LineWidth', 1.5);
ylabel('AC 下调潜力 (kW)'); title('AC 场景集 (下调)'); xlim([6 30]);

% 绘制 EV 上调能力
subplot(2, 2, 3); hold on;
plot(time_axis, Scenarios_EV_Up, 'Color', [0.6 0.8 1 0.1]);
plot(time_axis, mean(Scenarios_EV_Up, 2), 'b-', 'LineWidth', 1.5);
xlabel('时间 (h)'); ylabel('EV 上调潜力 (kW)'); title('EV 场景集 (上调)'); xlim([6 30]);

% 绘制 EV 下调能力
subplot(2, 2, 4); hold on;
plot(time_axis, Scenarios_EV_Down, 'Color', [1 0.8 0.6 0.1]);
plot(time_axis, mean(Scenarios_EV_Down, 2), 'r-', 'LineWidth', 1.5);
xlabel('时间 (h)'); ylabel('EV 下调潜力 (kW)'); title('EV 场景集 (下调)'); xlim([6 30]);