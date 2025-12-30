clc; clear; close all;

%% ================= 1. 参数设置 =================
num_scenarios = 1000;    
alpha = 0.05;           

% 仿真时间设置
simulation_start_hour = 8;
simulation_end_hour   = 32; 
dt = 15/60;                  
time_points = simulation_start_hour:dt:simulation_end_hour;
T_steps = length(time_points);

if ~exist('ini_data_mont', 'dir'), mkdir('ini_data_mont'); end

fprintf('==========================================================\n');
fprintf('步骤 0: 准备数据生成环境...\n');
fprintf('==========================================================\n');

base_ev_file = 'ini_data_mont/EV_Base_Fixed_mix_01_1000_8am.xlsx';
if ~exist(base_ev_file, 'file')
    generateEVParameters_real_8am(base_ev_file, 2000, 0, 'AreaType', '居民区');
end
BaseEVTable = readtable(base_ev_file);

fprintf('==========================================================\n');
fprintf('步骤 1: 执行蒙特卡洛模拟 (并行加速: %d 个场景)...\n', num_scenarios);
fprintf('==========================================================\n');

tic; 

%% ================= 2. 蒙特卡洛模拟循环 (并行化) =================
temp_AC_Up     = cell(1, num_scenarios);
temp_AC_Down   = cell(1, num_scenarios);
temp_AC_Base   = cell(1, num_scenarios); 
temp_AC_Params = cell(1, num_scenarios); % [新增]

temp_EV_Up     = cell(1, num_scenarios);
temp_EV_Down   = cell(1, num_scenarios);
temp_EV_Base   = cell(1, num_scenarios);
temp_EV_Params = cell(1, num_scenarios); % [新增]
% [新增] 能量边界临时存储
temp_EV_E_Up   = cell(1, num_scenarios);
temp_EV_E_Down = cell(1, num_scenarios);

if isempty(gcp('nocreate')), parpool; end

parfor s = 1:num_scenarios
    current_seed = 2024 + s; 
    temp_ev_file = sprintf('ini_data_mont/EV_Data_Temp_%d.xlsx', s);
    temp_ac_file = sprintf('ini_data_mont/AC_Data_Temp_%d.xlsx', s);
    
    generateACParameters_less_temp(temp_ac_file, 2000); 
    CurrentEVTable = randomize_ev_times_only_mix_8am(BaseEVTable, current_seed);
    writetable(CurrentEVTable, temp_ev_file);
    
    % --- AC 仿真 (接收基线和参数) ---
    [ac_up, ac_down, ac_base, ac_params] = run_AC_simulation_MC_soc(current_seed, temp_ac_file);
    
    % --- EV 仿真 (接收基线、参数以及新增的能量边界) ---
    [ev_up, ev_down, ev_power_profile, ev_params, ev_e_up, ev_e_down] ...
        = run_EV_simulation_MC_soc_bound(current_seed, temp_ev_file);
    
    % 存储结果
    temp_AC_Up{s}   = ac_up;
    temp_AC_Down{s} = ac_down;
    temp_AC_Base{s} = ac_base; 
    temp_AC_Params{s} = ac_params; 
    
    temp_EV_Up{s}   = ev_up;
    temp_EV_Down{s} = ev_down;
    temp_EV_Base{s} = ev_power_profile; 
    temp_EV_Params{s} = ev_params;
    
    % [新增] 存储能量边界
    temp_EV_E_Up{s}   = ev_e_up;
    temp_EV_E_Down{s} = ev_e_down;
    
    if exist(temp_ev_file, 'file'), delete(temp_ev_file); end
    if exist(temp_ac_file, 'file'), delete(temp_ac_file); end
    
    if mod(s, 10) == 0, fprintf('已完成场景: %d / %d\n', s, num_scenarios); end
end

toc;
fprintf('场景生成完成。正在整理数据...\n');

%% ================= 3. 数据整理与对齐 =================
Scenarios_AC_Up   = zeros(T_steps, num_scenarios);
Scenarios_AC_Down = zeros(T_steps, num_scenarios);
Scenarios_AC_Base = zeros(T_steps, num_scenarios); 

Scenarios_EV_Up   = zeros(T_steps, num_scenarios);
Scenarios_EV_Down = zeros(T_steps, num_scenarios);
Scenarios_EV_Base = zeros(T_steps, num_scenarios); 

% [新增] 能量边界数据整理
Scenarios_EV_E_Up   = zeros(T_steps, num_scenarios);
Scenarios_EV_E_Down = zeros(T_steps, num_scenarios);

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
    
    % EV 对齐 (包含功率与能量)
    ev_up = temp_EV_Up{s}; ev_down = temp_EV_Down{s}; ev_base = temp_EV_Base{s};
    ev_e_up = temp_EV_E_Up{s}; ev_e_down = temp_EV_E_Down{s}; % [新增]
    
    len_ev = length(ev_up);
    if len_ev >= T_steps
        Scenarios_EV_Up(:, s) = ev_up(1:T_steps);
        Scenarios_EV_Down(:, s) = ev_down(1:T_steps);
        Scenarios_EV_Base(:, s) = ev_base(1:T_steps);
        
        Scenarios_EV_E_Up(:, s) = ev_e_up(1:T_steps);     % [新增]
        Scenarios_EV_E_Down(:, s) = ev_e_down(1:T_steps); % [新增]
    else
        Scenarios_EV_Up(1:len_ev, s) = ev_up;
        Scenarios_EV_Down(1:len_ev, s) = ev_down;
        Scenarios_EV_Base(1:len_ev, s) = ev_base;
        
        Scenarios_EV_E_Up(1:len_ev, s) = ev_e_up;       % [新增]
        Scenarios_EV_E_Down(1:len_ev, s) = ev_e_down;   % [新增]
    end
end

%% ================= 4. 可靠调节域 & 可靠基线提取 =================
fprintf('\n提取可靠调节域 (置信度 %.0f%%)...\n', (1-alpha)*100);

Reliable_AC_Up = quantile(Scenarios_AC_Up, alpha, 2);
Reliable_AC_Down = quantile(Scenarios_AC_Down, 1-alpha, 2);
Reliable_EV_Up = quantile(Scenarios_EV_Up, alpha, 2);
Reliable_EV_Down = quantile(Scenarios_EV_Down, 1-alpha, 2);

Reliable_AC_Base = mean(Scenarios_AC_Base, 2); 
Reliable_EV_Base = mean(Scenarios_EV_Base, 2); 
Reliable_AC_Base_95 = quantile(Scenarios_AC_Base, 0.95, 2); 
Reliable_EV_Base_95 = quantile(Scenarios_EV_Base, 0.95, 2);

% [新增] 提取可靠能量包络 (Reliable Energy Envelope)
% 由于 E_Up 和 E_Down 都是正值代表"容量大小"，我们使用 alpha (例如0.05) 分位数
% 含义：在95%的场景下，我们至少拥有这么大的累积调节能量容量
Reliable_EV_E_Up = quantile(Scenarios_EV_E_Up, alpha, 2);
Reliable_EV_E_Down = quantile(Scenarios_EV_E_Down, alpha, 2); % 同样取保守下界(容量的最小值)

%% ================= 4.5 代表性场景聚类 (K-means) =================
fprintf('正在执行 K-means 聚类提取代表性基线场景...\n');
n_clusters = 5; 

X_features = [Scenarios_AC_Base; Scenarios_EV_Base]'; 

rng(2024); 
[idx_cluster, C_centers] = kmeans(X_features, n_clusters, 'Replicates', 5, 'MaxIter', 1000);

Reliable_AC_Base_Clusters = C_centers(:, 1:T_steps)'; 
Reliable_EV_Base_Clusters = C_centers(:, (T_steps+1):end)';

cluster_counts = groupcounts(idx_cluster);
fprintf('  > 聚类完成，各场景归类统计: %s\n', mat2str(cluster_counts'));

%% ================= 4.6 [新增] 计算聚合参数的期望值 (用于优化模型) =================
fprintf('正在计算聚合模型参数的期望值 (A, B, C)...\n');

% --- 处理 AC 参数 (通常是标量，但也可能是时变，取均值) ---
ac_A_list = cellfun(@(x) x.A, temp_AC_Params);
ac_B_list = cellfun(@(x) x.B, temp_AC_Params);
ac_C_list = cellfun(@(x) x.C, temp_AC_Params);

Reliable_AC_Params.A = mean(ac_A_list);
Reliable_AC_Params.B = mean(ac_B_list);
Reliable_AC_Params.C = mean(ac_C_list); 

% --- 处理 EV 参数 (时变向量) ---
all_ev_A = [temp_EV_Params{:}]; 
ev_A_mat = [all_ev_A.A]; 
ev_B_mat = [all_ev_A.B]; 
ev_C_mat = [all_ev_A.C];
ev_Cap_mat = [all_ev_A.C_total];

Reliable_EV_Params.A = mean(ev_A_mat, 2); 
Reliable_EV_Params.B = mean(ev_B_mat, 2); 
Reliable_EV_Params.C = mean(ev_C_mat, 2); 
Reliable_EV_Params.C_total = mean(ev_Cap_mat);

fprintf('  > AC 参数: A=%.4f, B=%.4f (Mean)\n', Reliable_AC_Params.A, Reliable_AC_Params.B);
fprintf('  > EV 参数: (时变向量已计算), C_total=%.2f kWh (Mean)\n', Reliable_EV_Params.C_total);

%% ================= 5. 结果保存 =================
save('reliable_regulation_domain_soc_bound.mat', ...
    'Scenarios_AC_Up', 'Scenarios_AC_Down', 'Reliable_AC_Up', 'Reliable_AC_Down', ...
    'Scenarios_EV_Up', 'Scenarios_EV_Down', 'Reliable_EV_Up', 'Reliable_EV_Down', ...
    'Scenarios_EV_E_Up', 'Scenarios_EV_E_Down', 'Reliable_EV_E_Up', 'Reliable_EV_E_Down', ... % [新增] 保存能量边界
    'Reliable_AC_Base', 'Reliable_EV_Base', 'Reliable_AC_Base_95','Reliable_EV_Base_95',... 
    'Reliable_AC_Base_Clusters', 'Reliable_EV_Base_Clusters', ...
    'Reliable_AC_Params', 'Reliable_EV_Params', ... 
    'time_points', 'alpha', 'num_scenarios');
fprintf('数据已保存 (包含聚合模型参数及能量边界)。\n');