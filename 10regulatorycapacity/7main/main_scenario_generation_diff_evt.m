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

% 文件夹准备
if ~exist('ini_data_mont', 'dir'), mkdir('ini_data_mont'); end

fprintf('==========================================================\n');
fprintf('步骤 0: 准备数据生成环境...\n');
fprintf('==========================================================\n');

% --- 0. 准备 EV 基准数据 (物理参数锁定) ---
% 我们生成一份基准文件，包含所有EV的固定物理参数(C, P_N, eta, E_ini等)
base_ev_file = 'ini_data_mont/EV_Base_Fixed.xlsx';

if ~exist(base_ev_file, 'file')
    % 生成 1000 辆 EV 的基准数据 (仅生成一次)
    generateEVParameters_real_6am(base_ev_file, 1000, 0, 'AreaType', '居民区');
    fprintf('已生成 EV 基准文件 (物理参数将固定): %s\n', base_ev_file);
else
    fprintf('使用现有 EV 基准文件: %s\n', base_ev_file);
end

% 将基准 EV 数据读入内存，避免在 parfor 中重复 IO，且作为随机化的基础
BaseEVTable = readtable(base_ev_file);


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
    
    % 生成随机种子 (控制本场景的随机性)
    current_seed = 2024 + s; 
    
    % 定义当前 worker 的临时文件名
    temp_ev_file = sprintf('ini_data_mont/EV_Data_Temp_%d.xlsx', s);
    temp_ac_file = sprintf('ini_data_mont/AC_Data_Temp_%d.xlsx', s);
    
    % --- 1. AC 数据生成 (保持原有规则：每次重新生成) ---
    % 每次调用都会生成全新的 AC 参数（ID、热阻等均可能变化，符合原逻辑）
    generateACParameters(temp_ac_file, 1000); 
    
    % --- 2. EV 数据生成 (修改规则：仅随机化时间) ---
    % 调用辅助函数，基于基准表仅修改 t_in 和 t_dep
    CurrentEVTable = randomize_ev_times_only(BaseEVTable, current_seed);
    % 保存为当前场景的输入文件
    writetable(CurrentEVTable, temp_ev_file);
    
    % --- 3. 执行仿真 ---
    
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
