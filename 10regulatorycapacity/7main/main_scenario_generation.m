%% main_scenario_generation.m
% 对应论文章节：
% 4.2.1 基于蒙特卡洛模拟的用户行为不确定性场景生成
% 4.2.2 基于分位数回归的聚合体可靠调节域边界提取
% 功能：生成 AC 和 EV 的随机场景，并提取 95% 置信度的可靠调节边界。
% 【修复版】：解决了 parfor 变量分类错误，使用 cell 数组进行中间存储。

clc; clear; close all;

%% ================= 1. 参数设置 =================
num_scenarios = 100;    % 场景生成数量 N (论文建议 > 100)
alpha = 0.05;           % 置信水平 alpha (对应 95% 可靠度)

% 仿真时间设置 (需与 run_ 函数内部保持一致)
simulation_start_hour = 6;
simulation_end_hour   = 30; % 次日 6:00
dt = 5/60;                  % 时间步长 (小时)
time_points = simulation_start_hour:dt:simulation_end_hour;
T_steps = length(time_points);

fprintf('==========================================================\n');
fprintf('步骤 1: 执行蒙特卡洛模拟 (共 %d 个场景)...\n', num_scenarios);
fprintf('==========================================================\n');

tic; 

%% ================= 2. 蒙特卡洛模拟循环 (4.2.1 节) =================
% 使用临时 Cell 数组存储 parfor 的结果
% 避免 "无法对 parfor 循环体中的变量进行分类" 错误
temp_AC_Up   = cell(1, num_scenarios);
temp_AC_Down = cell(1, num_scenarios);
temp_EV_Up   = cell(1, num_scenarios);
temp_EV_Down = cell(1, num_scenarios);
temp_EV_Energy = zeros(1, num_scenarios); % 标量数组通常可以直接在 parfor 中切片

% 使用 parfor 加速
parfor s = 1:num_scenarios
    
    if mod(s, 10) == 0
        fprintf('正在生成场景: %d / %d\n', s, num_scenarios);
    end
    
    % 生成唯一随机种子，控制该场景下的所有用户行为不确定性
    % (影响：用户参与度、响应偏差、EV入离网时间等均受此种子影响)
    current_seed = 2024 + s; 
    
    % --- 2.1 生成 AC 场景 ---
    % 调用 run_AC_simulation_MC (基准指令为0，获取物理边界)
    [ac_up, ac_down] = run_AC_simulation_MC(current_seed);
    
    % --- 2.2 生成 EV 场景 ---
    % 调用 run_EV_simulation_MC (基准指令为0，获取物理边界)
    % 获取 EV_Power_Sum 以计算总能量需求
    [ev_up, ev_down, ev_power_profile] = run_EV_simulation_MC(current_seed);
    
    % 计算该场景下的 EV 总能量需求 (kWh)
    % E = sum(P(t) * dt)
    ev_energy = sum(ev_power_profile) * dt;
    
    % --- 2.3 存入临时 Cell ---
    % 注意：这里直接存储原始结果，数据对齐在循环外进行
    temp_AC_Up{s}   = ac_up;
    temp_AC_Down{s} = ac_down;
    temp_EV_Up{s}   = ev_up;
    temp_EV_Down{s} = ev_down;
    temp_EV_Energy(s) = ev_energy;
end

toc;
fprintf('场景生成完成。正在整理数据...\n');

%% ================= 3. 数据整理与对齐 =================
% 将 Cell 数组转换为矩阵，并确保维度一致
Scenarios_AC_Up   = zeros(T_steps, num_scenarios);
Scenarios_AC_Down = zeros(T_steps, num_scenarios);
Scenarios_EV_Up   = zeros(T_steps, num_scenarios);
Scenarios_EV_Down = zeros(T_steps, num_scenarios);
Scenarios_EV_Energy_Total = temp_EV_Energy;

for s = 1:num_scenarios
    % 获取当前场景的数据
    ac_up = temp_AC_Up{s};
    ac_down = temp_AC_Down{s};
    ev_up = temp_EV_Up{s};
    ev_down = temp_EV_Down{s};
    
    % 对齐长度 (截断或补零)
    % AC Up
    len = length(ac_up);
    if len >= T_steps
        Scenarios_AC_Up(:, s) = ac_up(1:T_steps);
    else
        Scenarios_AC_Up(1:len, s) = ac_up;
        % 不足部分保持为0
    end
    
    % AC Down
    len = length(ac_down);
    if len >= T_steps
        Scenarios_AC_Down(:, s) = ac_down(1:T_steps);
    else
        Scenarios_AC_Down(1:len, s) = ac_down;
    end
    
    % EV Up
    len = length(ev_up);
    if len >= T_steps
        Scenarios_EV_Up(:, s) = ev_up(1:T_steps);
    else
        Scenarios_EV_Up(1:len, s) = ev_up;
    end
    
    % EV Down
    len = length(ev_down);
    if len >= T_steps
        Scenarios_EV_Down(:, s) = ev_down(1:T_steps);
    else
        Scenarios_EV_Down(1:len, s) = ev_down;
    end
end

% 清理临时变量以释放内存
clear temp_AC_Up temp_AC_Down temp_EV_Up temp_EV_Down;

%% ================= 4. 可靠调节域边界提取 (4.2.2 节) =================
fprintf('\n==========================================================\n');
fprintf('步骤 2: 基于分位数回归提取可靠调节域 (置信度 %.0f%%)...\n', (1-alpha)*100);
fprintf('==========================================================\n');

% --- 4.1 提取 AC 可靠边界 ---
% 上调能力 (正值): 取 alpha (5%) 分位数。
% 物理意义：保证 95% 的场景下，实际可用上调能力 >= 可靠边界。
Reliable_AC_Up = quantile(Scenarios_AC_Up, alpha, 2);

% 下调能力 (负值): 取 1-alpha (95%) 分位数。
% 物理意义：由于数值为负，95%分位数是“最不负”（最接近0）的边界。
% 保证 95% 的场景下，实际可用下调能力（绝对值） >= 可靠边界（绝对值）。
Reliable_AC_Down = quantile(Scenarios_AC_Down, 1-alpha, 2);

% --- 4.2 提取 EV 可靠边界 ---
% 原理同上
Reliable_EV_Up = quantile(Scenarios_EV_Up, alpha, 2);
Reliable_EV_Down = quantile(Scenarios_EV_Down, 1-alpha, 2);

% --- 4.3 提取 EV 能量约束边界 ---
% 能量需求 (正值): 取 1-alpha (95%) 分位数。
% 物理意义：保证 95% 的场景下，预留的能量足以满足用户的充电需求（覆盖高耗能场景）。
Reliable_EV_Energy_Need = quantile(Scenarios_EV_Energy_Total, 1-alpha);

fprintf('提取完成。\n');
fprintf('  EV 95%% 置信度下的最小能量需求: %.2f kWh\n', Reliable_EV_Energy_Need);

%% ================= 5. 结果保存 =================
output_file = 'reliable_regulation_domain.mat';
save(output_file, ...
    'Scenarios_AC_Up', 'Scenarios_AC_Down', ...
    'Scenarios_EV_Up', 'Scenarios_EV_Down', ...
    'Reliable_AC_Up', 'Reliable_AC_Down', ...
    'Reliable_EV_Up', 'Reliable_EV_Down', ...
    'Reliable_EV_Energy_Need', ...
    'time_points', 'alpha', 'num_scenarios');
fprintf('所有数据已保存至: %s\n', output_file);

%% ================= 6. 可视化 (复现论文图表风格) =================
figure('Name', '可靠调节域提取结果', 'Position', [100, 100, 1200, 800], 'Color', 'w');

% --- 子图 1: AC 调节域 ---
subplot(2, 1, 1);
hold on;
% 绘制所有场景的随机分布 (浅色背景)
plot(time_points, Scenarios_AC_Up, 'Color', [0.6, 0.8, 1, 0.15], 'HandleVisibility', 'off');
plot(time_points, Scenarios_AC_Down, 'Color', [1, 0.6, 0.6, 0.15], 'HandleVisibility', 'off');
% 绘制可靠边界 (深色粗线)
p1 = plot(time_points, Reliable_AC_Up, 'b-', 'LineWidth', 2, 'DisplayName', '可靠上调边界 (95%置信)');
p2 = plot(time_points, Reliable_AC_Down, 'r-', 'LineWidth', 2, 'DisplayName', '可靠下调边界 (95%置信)');
% 绘制 0 线
yline(0, 'k--');

title('空调 (AC) 聚合体可靠调节域');
xlabel('时间 (h)'); ylabel('功率 (kW)');
legend([p1, p2], 'Location', 'best');
grid on;
xlim([simulation_start_hour, simulation_end_hour]);
set(gca, 'FontSize', 12);

% --- 子图 2: EV 调节域 ---
subplot(2, 1, 2);
hold on;
% 绘制所有场景的随机分布 (浅色背景)
plot(time_points, Scenarios_EV_Up, 'Color', [0.6, 0.8, 1, 0.15], 'HandleVisibility', 'off');
plot(time_points, Scenarios_EV_Down, 'Color', [1, 0.6, 0.6, 0.15], 'HandleVisibility', 'off');
% 绘制可靠边界 (深色粗线)
p3 = plot(time_points, Reliable_EV_Up, 'b-', 'LineWidth', 2, 'DisplayName', '可靠上调边界 (95%置信)');
p4 = plot(time_points, Reliable_EV_Down, 'r-', 'LineWidth', 2, 'DisplayName', '可靠下调边界 (95%置信)');
% 绘制 0 线
yline(0, 'k--');

title('电动汽车 (EV) 聚合体可靠调节域');
xlabel('时间 (h)'); ylabel('功率 (kW)');
legend([p3, p4], 'Location', 'best');
grid on;
xlim([simulation_start_hour, simulation_end_hour]);
set(gca, 'FontSize', 12);

fprintf('绘图完成。\n');