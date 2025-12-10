clc; clear; close all;

%% ================= 1. 参数设置 =================
num_scenarios = 1000;    
alpha = 0.05;           

% 仿真时间设置
simulation_start_hour = 6;
simulation_end_hour   = 30; 
dt = 5/60;                  
time_points = simulation_start_hour:dt:simulation_end_hour;
T_steps = length(time_points);

if ~exist('ini_data_mont', 'dir'), mkdir('ini_data_mont'); end

fprintf('==========================================================\n');
fprintf('步骤 0: 准备数据生成环境...\n');
fprintf('==========================================================\n');

base_ev_file = 'ini_data_mont/EV_Base_Fixed_mix_01_1000.xlsx';
if ~exist(base_ev_file, 'file')
    generateEVParameters_real_6am(base_ev_file, 2000, 0, 'AreaType', '居民区');
end
BaseEVTable = readtable(base_ev_file);

fprintf('==========================================================\n');
fprintf('步骤 1: 执行蒙特卡洛模拟 (并行加速: %d 个场景)...\n', num_scenarios);
fprintf('==========================================================\n');

tic; 

%% ================= 2. 蒙特卡洛模拟循环 (并行化) =================
temp_AC_Up     = cell(1, num_scenarios);
temp_AC_Down   = cell(1, num_scenarios);
temp_AC_Base   = cell(1, num_scenarios); % 新增：AC基线

temp_EV_Up     = cell(1, num_scenarios);
temp_EV_Down   = cell(1, num_scenarios);
temp_EV_Base   = cell(1, num_scenarios); % 新增：EV基线 (即充电负荷)

if isempty(gcp('nocreate')), parpool; end

parfor s = 1:num_scenarios
    current_seed = 2024 + s; 
    temp_ev_file = sprintf('ini_data_mont/EV_Data_Temp_%d.xlsx', s);
    temp_ac_file = sprintf('ini_data_mont/AC_Data_Temp_%d.xlsx', s);
    
    generateACParameters_less_temp(temp_ac_file, 2000); 
    CurrentEVTable = randomize_ev_times_only_mix(BaseEVTable, current_seed);
    writetable(CurrentEVTable, temp_ev_file);
    
    % --- AC 仿真 (接收基线) ---
    [ac_up, ac_down, ac_base] = run_AC_simulation_MC(current_seed, temp_ac_file);
    
    % --- EV 仿真 (接收基线) ---
    % 注意：ev_power_profile 即为该场景下的 EV 基线充电负荷
    [ev_up, ev_down, ev_power_profile] = run_EV_simulation_MC(current_seed, temp_ev_file);
    
    % 存储结果
    temp_AC_Up{s}   = ac_up;
    temp_AC_Down{s} = ac_down;
    temp_AC_Base{s} = ac_base; % 保存AC基线
    
    temp_EV_Up{s}   = ev_up;
    temp_EV_Down{s} = ev_down;
    temp_EV_Base{s} = ev_power_profile; % 保存EV基线 (Power Profile)
    
    if exist(temp_ev_file, 'file'), delete(temp_ev_file); end
    if exist(temp_ac_file, 'file'), delete(temp_ac_file); end
    
    if mod(s, 10) == 0, fprintf('已完成场景: %d / %d\n', s, num_scenarios); end
end

toc;
fprintf('场景生成完成。正在整理数据...\n');

%% ================= 3. 数据整理与对齐 =================
Scenarios_AC_Up   = zeros(T_steps, num_scenarios);
Scenarios_AC_Down = zeros(T_steps, num_scenarios);
Scenarios_AC_Base = zeros(T_steps, num_scenarios); % 新增

Scenarios_EV_Up   = zeros(T_steps, num_scenarios);
Scenarios_EV_Down = zeros(T_steps, num_scenarios);
Scenarios_EV_Base = zeros(T_steps, num_scenarios); % 新增

for s = 1:num_scenarios
    % AC 对齐
    ac_up = temp_AC_Up{s}; ac_down = temp_AC_Down{s}; ac_base = temp_AC_Base{s};
    len_ac = length(ac_up);
    if len_ac >= T_steps
        Scenarios_AC_Up(:, s) = ac_up(1:T_steps);
        Scenarios_AC_Down(:, s) = ac_down(1:T_steps);
        Scenarios_AC_Base(:, s) = ac_base(1:T_steps);
    else
        Scenarios_AC_Up(1:len_ac, s) = ac_up;
        Scenarios_AC_Down(1:len_ac, s) = ac_down;
        Scenarios_AC_Base(1:len_ac, s) = ac_base;
    end
    
    % EV 对齐
    ev_up = temp_EV_Up{s}; ev_down = temp_EV_Down{s}; ev_base = temp_EV_Base{s};
    len_ev = length(ev_up);
    if len_ev >= T_steps
        Scenarios_EV_Up(:, s) = ev_up(1:T_steps);
        Scenarios_EV_Down(:, s) = ev_down(1:T_steps);
        Scenarios_EV_Base(:, s) = ev_base(1:T_steps);
    else
        Scenarios_EV_Up(1:len_ev, s) = ev_up;
        Scenarios_EV_Down(1:len_ev, s) = ev_down;
        Scenarios_EV_Base(1:len_ev, s) = ev_base;
    end
end

%% ================= 4. 可靠调节域 & 可靠基线提取 =================
fprintf('\n提取可靠调节域 (置信度 %.0f%%)...\n', (1-alpha)*100);

% 调节能力 (Up/Down) 取分位数 (保守估计)
Reliable_AC_Up = quantile(Scenarios_AC_Up, alpha, 2);
Reliable_AC_Down = quantile(Scenarios_AC_Down, 1-alpha, 2);
Reliable_EV_Up = quantile(Scenarios_EV_Up, alpha, 2);
Reliable_EV_Down = quantile(Scenarios_EV_Down, 1-alpha, 2);

% [新增] 基础功率 (Base) 取均值或中位数 (作为系统运行的基准参考)
% 使用均值能反映平均负荷水平，用于潮流计算更合理
Reliable_AC_Base = mean(Scenarios_AC_Base, 2); 
Reliable_EV_Base = mean(Scenarios_EV_Base, 2); 

fprintf('提取完成。\n');

%% ================= 5. 结果保存 =================
% 增加了 Reliable_AC_Base 和 Reliable_EV_Base
save('reliable_regulation_domain_1000_mix_pbase_less.mat', ...
    'Scenarios_AC_Up', 'Scenarios_AC_Down', 'Reliable_AC_Up', 'Reliable_AC_Down', ...
    'Scenarios_EV_Up', 'Scenarios_EV_Down', 'Reliable_EV_Up', 'Reliable_EV_Down', ...
    'Reliable_AC_Base', 'Reliable_EV_Base', ... % <--- 关键新增
    'time_points', 'alpha', 'num_scenarios');
fprintf('数据已保存。\n');