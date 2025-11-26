%% main_scenario_generation.m
% 对应论文章节：4.2.1 & 4.2.2
% 功能：基于蒙特卡洛模拟生成场景 (并行加速版)

clc; clear; close all;

%% ================= 1. 参数设置 =================
num_scenarios = 100;    % 场景生成数量
alpha = 0.05;           % 置信水平

% 仿真时间设置
simulation_start_hour = 6;
simulation_end_hour   = 30; 
dt = 5/60;                  
time_points = simulation_start_hour:dt:simulation_end_hour;
T_steps = length(time_points);

fprintf('==========================================================\n');
fprintf('步骤 1: 执行蒙特卡洛模拟 (并行加速: %d 个场景)...\n', num_scenarios);
fprintf('==========================================================\n');

tic; 

%% ================= 2. 蒙特卡洛模拟循环 (并行化) =================
temp_AC_Up     = cell(1, num_scenarios);
temp_AC_Down   = cell(1, num_scenarios);
temp_EV_Up     = cell(1, num_scenarios);
temp_EV_Down   = cell(1, num_scenarios);
temp_EV_Energy = zeros(1, num_scenarios); 

% 检查并行池
if isempty(gcp('nocreate')), parpool; end

% 使用 parfor 进行并行计算
parfor s = 1:num_scenarios
    
    % --- 1. 为当前 Worker 生成唯一的临时文件名 ---
    % 加上 Worker ID 或 循环索引 s 以防止文件冲突
    temp_ev_file = sprintf('ini_data_mont/EV_Data_Temp_%d.xlsx', s);
    temp_ac_file = sprintf('ini_data_mont/AC_Data_Temp_%d.xlsx', s);
    
    % 生成随机种子
    current_seed = 2024 + s; 
    
    % --- 2. 生成数据文件 (指定唯一文件名) ---
    % 注意：生成函数需要支持传入文件名参数
    
    % 生成 EV 数据 (generateEVParameters_real_6am 原生支持传入文件名)
    generateEVParameters_real_6am(temp_ev_file, 1000, 0, 'AreaType', '居民区');
    
    % 生成 AC 数据 (需要修改 generateExampleExcel_real_24 以支持传入文件名)
    % 这里假设您修改了该函数，或者使用下面的修改版 generateACParameters_Safe
    generateACParameters(temp_ac_file, 1000); 
    
    % --- 3. 执行仿真 (传入对应的唯一文件名) ---
    % 注意：run_ 函数也需要修改以接受文件名参数
    
    % AC 仿真
    [ac_up, ac_down] = run_AC_simulation_MC(current_seed, temp_ac_file);
    
    % EV 仿真
    [ev_up, ev_down, ev_power_profile] = run_EV_simulation_MC(current_seed, temp_ev_file);
    
    % 计算 EV 能量
    ev_energy = sum(ev_power_profile) * dt;
    
    % --- 4. 存储结果 ---
    temp_AC_Up{s}   = ac_up;
    temp_AC_Down{s} = ac_down;
    temp_EV_Up{s}   = ev_up;
    temp_EV_Down{s} = ev_down;
    temp_EV_Energy(s) = ev_energy;
    
    % --- 5. 清理临时文件 ---
    if exist(temp_ev_file, 'file'), delete(temp_ev_file); end
    if exist(temp_ac_file, 'file'), delete(temp_ac_file); end
    
    if mod(s, 10) == 0
        fprintf('已完成场景: %d / %d\n', s, num_scenarios);
    end
end

toc;
fprintf('场景生成完成。正在整理数据...\n');

%% ================= 3. 数据整理与对齐 (保持不变) =================
Scenarios_AC_Up   = zeros(T_steps, num_scenarios);
Scenarios_AC_Down = zeros(T_steps, num_scenarios);
Scenarios_EV_Up   = zeros(T_steps, num_scenarios);
Scenarios_EV_Down = zeros(T_steps, num_scenarios);
Scenarios_EV_Energy_Total = temp_EV_Energy;

for s = 1:num_scenarios
    % AC 对齐
    ac_up = temp_AC_Up{s}; ac_down = temp_AC_Down{s};
    len_ac = length(ac_up);
    if len_ac >= T_steps
        Scenarios_AC_Up(:, s) = ac_up(1:T_steps);
        Scenarios_AC_Down(:, s) = ac_down(1:T_steps);
    else
        Scenarios_AC_Up(1:len_ac, s) = ac_up;
        Scenarios_AC_Down(1:len_ac, s) = ac_down;
    end
    
    % EV 对齐
    ev_up = temp_EV_Up{s}; ev_down = temp_EV_Down{s};
    len_ev = length(ev_up);
    if len_ev >= T_steps
        Scenarios_EV_Up(:, s) = ev_up(1:T_steps);
        Scenarios_EV_Down(:, s) = ev_down(1:T_steps);
    else
        Scenarios_EV_Up(1:len_ev, s) = ev_up;
        Scenarios_EV_Down(1:len_ev, s) = ev_down;
    end
end

clear temp_AC_Up temp_AC_Down temp_EV_Up temp_EV_Down;

%% ================= 4. 可靠调节域边界提取 (保持不变) =================
fprintf('\n提取可靠调节域 (置信度 %.0f%%)...\n', (1-alpha)*100);
Reliable_AC_Up = quantile(Scenarios_AC_Up, alpha, 2);
Reliable_AC_Down = quantile(Scenarios_AC_Down, 1-alpha, 2);
Reliable_EV_Up = quantile(Scenarios_EV_Up, alpha, 2);
Reliable_EV_Down = quantile(Scenarios_EV_Down, 1-alpha, 2);
Reliable_EV_Energy_Need = quantile(Scenarios_EV_Energy_Total, 1-alpha);

fprintf('提取完成。EV 最小能量需求: %.2f kWh\n', Reliable_EV_Energy_Need);

%% ================= 5. 结果保存 (保持不变) =================
save('reliable_regulation_domain.mat', ...
    'Scenarios_AC_Up', 'Scenarios_AC_Down', ...
    'Scenarios_EV_Up', 'Scenarios_EV_Down', ...
    'Reliable_AC_Up', 'Reliable_AC_Down', ...
    'Reliable_EV_Up', 'Reliable_EV_Down', ...
    'Reliable_EV_Energy_Need', ...
    'time_points', 'alpha', 'num_scenarios');
fprintf('数据已保存。\n');

%% ================= 6. 可视化 (保持不变) =================
figure('Name', '可靠调节域提取结果', 'Position', [100, 100, 1200, 800], 'Color', 'w');
subplot(2, 1, 1); hold on;
plot(time_points, Scenarios_AC_Up, 'Color', [0.6, 0.8, 1, 0.15], 'HandleVisibility', 'off');
plot(time_points, Scenarios_AC_Down, 'Color', [1, 0.6, 0.6, 0.15], 'HandleVisibility', 'off');
p1 = plot(time_points, Reliable_AC_Up, 'b-', 'LineWidth', 2, 'DisplayName', '可靠上调边界');
p2 = plot(time_points, Reliable_AC_Down, 'r-', 'LineWidth', 2, 'DisplayName', '可靠下调边界');
yline(0, 'k--'); title('空调 (AC) 聚合体可靠调节域'); legend([p1, p2], 'Location', 'best'); grid on; xlim([simulation_start_hour, simulation_end_hour]);

subplot(2, 1, 2); hold on;
plot(time_points, Scenarios_EV_Up, 'Color', [0.6, 0.8, 1, 0.15], 'HandleVisibility', 'off');
plot(time_points, Scenarios_EV_Down, 'Color', [1, 0.6, 0.6, 0.15], 'HandleVisibility', 'off');
p3 = plot(time_points, Reliable_EV_Up, 'b-', 'LineWidth', 2, 'DisplayName', '可靠上调边界');
p4 = plot(time_points, Reliable_EV_Down, 'r-', 'LineWidth', 2, 'DisplayName', '可靠下调边界');
yline(0, 'k--'); title('电动汽车 (EV) 聚合体可靠调节域'); legend([p3, p4], 'Location', 'best'); grid on; xlim([simulation_start_hour, simulation_end_hour]);