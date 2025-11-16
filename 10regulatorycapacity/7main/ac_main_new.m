% AC_main_stateful_aggregation.m
% 这是一个重构版本，实现了大论文中 2.4.1 节的
% "面向SOC状态一致" 的聚合与指令分解仿真。
% (基于 AC_main.m 和 ac_simulation_block.m 重构)

clear; close all; clc;

tic; % 启动一个总计时器

%% 1. 系统初始化
rng(2023, 'Threefry'); % 固定随机种子
T_total = 24; % 总时长（小时）
dt = 5/60;    % 时间分辨率（小时）
time_points = 0:dt:T_total; % 仿真时间点
T_steps_total = length(time_points);
steps_per_hour = round(1/dt);
num_hours = floor(T_steps_total / steps_per_hour);

base_price = 30; % 基础电价（元/kWh）

%% 2. 初始化 AC 参数
acFile = 'AC_template.xlsx'; % 确保 AC_template1.xlsx 在路径中
fprintf('正在初始化空调参数...\n');
ACs = initializeACsFromExcel(acFile);
num_AC = length(ACs);

% 备份原始设置
for i = 1:num_AC
    ACs(i).Tset_original = ACs(i).Tset;
    ACs(i).Tmax_original = ACs(i).Tmax;
    ACs(i).Tmin_original = ACs(i).Tmin;
    if ~isfield(ACs(i), 'p_incentive')
        ACs(i).p_incentive = round(60*rand(), 1);
    end
end
fprintf('加载了 %d 台空调。\n', num_AC);

%% 3. 激励响应参数
p_min = 15;        
p_max = 50;        
p_min_prime = 10;  
p_max_prime = 40;  
T_set_max = 3;     

% 定义要仿真的价格点
p_incentive_range_full = linspace(0, 50, 25);
plot_price_indices = [1, 7, 13, 19, 25]; % 选择几个代表性价格

if length(p_incentive_range_full) < max(plot_price_indices)
    plot_price_indices = round(linspace(1, length(p_incentive_range_full), min(5, length(p_incentive_range_full))));
end
num_prices_to_simulate = length(plot_price_indices);
fprintf('将对 %d 个价格场景进行状态化仿真...\n', num_prices_to_simulate);

%% 4. 初始化数据存储 (用于绘图)
all_AC_Up = zeros(length(time_points), num_prices_to_simulate);
all_AC_Down = zeros(length(time_points), num_prices_to_simulate);

%% 5. 主循环：遍历选定的价格场景
for k = 1:num_prices_to_simulate
    p_idx = plot_price_indices(k); 
    current_p = p_incentive_range_full(p_idx); 
    fprintf('\n== 仿真价格场景 %d/%d (Price: %.1f 元) ==\n', k, num_prices_to_simulate, current_p);

    temp_ACs_for_price_scenario = ACs; 

    % --- 5.1: 预计算 (alpha, beta, gamma) 和 T_ja ---
    fprintf('  Step 5.1: 预计算单体参数...\n');
    parfor i = 1:num_AC
        % 重置
        temp_ACs_for_price_scenario(i).Tset = temp_ACs_for_price_scenario(i).Tset_original;
        temp_ACs_for_price_scenario(i).Tmax = temp_ACs_for_price_scenario(i).Tmax_original;
        temp_ACs_for_price_scenario(i).Tmin = temp_ACs_for_price_scenario(i).Tmin_original;

        % 计算参与度
        participation = calculateParticipation(current_p, base_price);
        [~, ~, deltaT_flex_magnitude] = incentiveTempAC(...
            current_p, p_min, p_max, p_min_prime, p_max_prime, T_set_max);
        temp_ACs_for_price_scenario(i).ptcp = (rand() < participation);

        % 调整温度范围
        if temp_ACs_for_price_scenario(i).ptcp
            temp_ACs_for_price_scenario(i).Tmax = temp_ACs_for_price_scenario(i).Tset_original + deltaT_flex_magnitude;
            temp_ACs_for_price_scenario(i).Tmin = temp_ACs_for_price_scenario(i).Tset_original - deltaT_flex_magnitude;
            % (可添加边界约束, e.g., max(18, ...), min(28, ...))
        end

        % 生成 T_ja
        base_ambient_temp = temp_ACs_for_price_scenario(i).Tset_original + 4*sin(2*pi*time_points/24);
        actual_temp_range = temp_ACs_for_price_scenario(i).Tmax - temp_ACs_for_price_scenario(i).Tmin;
        if abs(actual_temp_range) < 1e-6; actual_temp_range = 0.1; end
        noise = 0.2 * actual_temp_range * randn(size(time_points));
        temp_ACs_for_price_scenario(i).T_ja = min(max(base_ambient_temp + noise, temp_ACs_for_price_scenario(i).Tmin), temp_ACs_for_price_scenario(i).Tmax);

        % 计算 ABC (实现 式 2-12)
        [alpha, beta, gamma] = calculateACABC_single(...
            temp_ACs_for_price_scenario(i).R, temp_ACs_for_price_scenario(i).C, temp_ACs_for_price_scenario(i).eta,...
            temp_ACs_for_price_scenario(i).Tmax, temp_ACs_for_price_scenario(i).Tmin, temp_ACs_for_price_scenario(i).Tset, dt);
        temp_ACs_for_price_scenario(i).alpha = alpha;
        temp_ACs_for_price_scenario(i).beta = beta;
        temp_ACs_for_price_scenario(i).gamma = gamma;
    end
    fprintf('  Step 5.1: 完成。\n');

    % --- 5.1B: 【新增】计算聚合参数 (流程图 步骤3) ---
    fprintf('  Step 5.1B: 计算聚合模型参数 (A, B, C)...\n');
    temp_ACs = temp_ACs_for_price_scenario; % 获取 parfor 的结果
    ACs_participating = temp_ACs([temp_ACs.ptcp]);
    num_AC_participating = length(ACs_participating);
    
    if num_AC_participating == 0
        fprintf('  该价格下无空调参与，跳过此场景。\n');
        all_AC_Up(:, k) = 0;
        all_AC_Down(:, k) = 0;
        continue; % 跳到下一个价格场景
    end

    % 计算聚合参数 (实现 式 2-37)
    AggParams = calculateAggregatedACParams(ACs_participating);
    fprintf('  聚合参数: A=%.4f, B=%.4f, C=%.4f\n', AggParams.A, AggParams.B, AggParams.C);

    % --- 5.1C: 【新增】生成电网指令 & 初始化状态 ---
    fprintf('  Step 5.1C: 生成电网指令并初始化状态...\n');
    P_grid_command_series = generate_hourly_regulation_signal(T_steps_total, steps_per_hour, num_hours, num_AC_participating);
    
    % 初始化状态向量 (N_participating x 1)
    CURRENT_SOC_AC = [ACs_participating.SOC]'; 
    
    % 初始化本地存储
    local_AC_Up_total_per_price = zeros(T_steps_total, 1);
    local_AC_Down_total_per_price = zeros(T_steps_total, 1);

    % --- 5.2: 【重构】时间步进仿真 (状态化) ---
    fprintf('  Step 5.2: 开始 %d 步的状态化仿真...\n', T_steps_total);
    
    for t_idx = 1:T_steps_total
        
        % 1. 获取当前聚合状态 SOC(t)
        SOC_agg_t = mean(CURRENT_SOC_AC, 'omitnan');
        
        % 2. 获取电网指令 ΔP_S (流程图 步骤4)
        Delta_P_S_command = P_grid_command_series(t_idx);
        
        % 3. 预测目标聚合 SOC(t+1) (流程图 步骤5)
        % (实现 式 2-36)
        SOC_target_next = AggParams.A * SOC_agg_t + AggParams.B * Delta_P_S_command + AggParams.C;
        % 约束目标
        SOC_target_next = max(0, min(1, SOC_target_next));

        % 4. 临时变量 (用于 parfor)
        temp_AC_Up_agg = 0;
        temp_AC_Down_agg = 0;
        temp_SOC_for_next_step = zeros(num_AC_participating, 1);

        parfor i = 1:num_AC_participating
            
            ac_i = ACs_participating(i);
            soc_current_i = CURRENT_SOC_AC(i); % 获取 SOC(t)
            
            % A. 计算当前物理潜力 (用于记录和约束)
            [P_plus, P_minus] = calculateACAdjustmentPotentia(...
                0, 1e6, -1e6, ... % P_base, P_max, P_min (在SOC约束下意义不大，设为宽松)
                ac_i.alpha, ac_i.beta, ac_i.gamma,...
                soc_current_i, dt); % (!!! 使用 dt, 不是 t_adj !!!)
            
            temp_AC_Up_agg = temp_AC_Up_agg + P_plus;
            temp_AC_Down_agg = temp_AC_Down_agg + P_minus;
            
            % B. 反解理论功率 ΔP_j (流程图 步骤6)
            % ΔP_j = (SOC_target - α_j*SOC(t) - γ_j) / β_j
            delta_Pj_theory = 0;
            if abs(ac_i.beta) > 1e-9
                delta_Pj_theory = (SOC_target_next - ac_i.alpha * soc_current_i - ac_i.gamma) / ac_i.beta;
            end
            
            % C. 裁剪指令至物理可行域
            delta_Pj_clipped = max(P_minus, min(P_plus, delta_Pj_theory));

            % D. 更新状态 (实现 式 2-10)
            soc_next_i = updateACSOC_single(soc_current_i, delta_Pj_clipped, ...
                ac_i.alpha, ac_i.beta, ac_i.gamma);
                
            temp_SOC_for_next_step(i) = soc_next_i; % 存储 SOC(t+1)
        end
        
        % 5. 存储当前时间步 t 的聚合潜力
        local_AC_Up_total_per_price(t_idx) = temp_AC_Up_agg;
        local_AC_Down_total_per_price(t_idx) = temp_AC_Down_agg;
        
        % 6. 更新状态向量用于下一循环
        CURRENT_SOC_AC = temp_SOC_for_next_step; % CURRENT_SOC_AC 现在是 SOC(t+1)
        
    end % 结束 t_idx 循环
    
    fprintf('  Step 5.2: 仿真完成。\n');
    
    % 存储此价格场景的结果
    all_AC_Up(:, k) = local_AC_Up_total_per_price;
    all_AC_Down(:, k) = local_AC_Down_total_per_price;
    
end % 结束 k 价格循环

%% 6. 保存数据到Excel (与原 AC_main.m 一致)
excel_filename = 'AC_Plot_Data_Selected_Stateful.xlsx';
fprintf('\n正在将仿真结果保存到 %s ...\n', excel_filename);

dataTable = table(time_points', 'VariableNames', {'Time_Hours'});

for k = 1:num_prices_to_simulate
    p_idx = plot_price_indices(k);
    current_price = p_incentive_range_full(p_idx);
    up_col_name = matlab.lang.makeValidName(sprintf('Up_Price_%.1f', current_price));
    down_col_name = matlab.lang.makeValidName(sprintf('Down_Price_%.1f', current_price));
    dataTable.(up_col_name) = all_AC_Up(:, k);
    dataTable.(down_col_name) = all_AC_Down(:, k);
end

try
    writetable(dataTable, excel_filename);
    fprintf('数据成功保存到 %s\n', excel_filename);
catch ME
    fprintf('*** 保存Excel文件时出错: %s ***\n', ME.message);
end

%% 7. 绘图 (与原 AC_main.m 一致)
fprintf('正在生成绘图...\n');
colors = lines(num_prices_to_simulate); 

% --- 7.1 上调潜力图 ---
figure('Name', 'AC Cluster Up-Regulation Capacity (Stateful)', 'Position', [100 100 1000 600]);
hold on; grid on;
legend_entries_up = {};
for k = 1:num_prices_to_simulate
    p_idx = plot_price_indices(k);
    current_price = p_incentive_range_full(p_idx);
    plot(time_points, all_AC_Up(:, k), 'LineWidth', 1.5, 'Color', colors(k,:), ...
         'DisplayName', sprintf('电价: %.1f 分/℃', current_price));
    legend_entries_up{end+1} = sprintf('电价: %.1f 分/℃', current_price);
end
hold off;
xlabel('时间', 'FontSize', 18);
ylabel('集群上调潜力 (kW)', 'FontSize', 18);
if ~isempty(legend_entries_up)
    legend(legend_entries_up, 'Location', 'best', 'FontSize', 16);
end
set(gca, 'FontSize', 16);
xticks([0, 6, 12, 18, 24]);
xticklabels({'00:00', '06:00', '12:00', '18:00', '24:00'});
xlim([0, 24]);
print('-dpng', '-r400', 'AC_Cluster_Up_Regulation_Stateful.png');
fprintf('空调集群上调潜力图已保存为 AC_Cluster_Up_Regulation_Stateful.png\n');

% --- 7.2 下调潜力图 ---
figure('Name', 'AC Cluster Down-Regulation Capacity (Stateful)', 'Position', [100 750 1000 600]);
hold on; grid on;
legend_entries_down = {};
for k = 1:num_prices_to_simulate
    p_idx = plot_price_indices(k);
    current_price = p_incentive_range_full(p_idx);
    plot(time_points, all_AC_Down(:, k), 'LineWidth', 1.5, 'Color', colors(k,:), ...
         'DisplayName', sprintf('电价: %.1f 分/℃', current_price));
    legend_entries_down{end+1} = sprintf('电价: %.1f 分/℃', current_price);
end
hold off;
xlabel('时间', 'FontSize', 18);
ylabel('集群下调潜力 (kW)', 'FontSize', 18);
if ~isempty(legend_entries_down)
    legend(legend_entries_down, 'Location', 'best', 'FontSize', 16);
end
set(gca, 'FontSize', 16);
xticks([0, 6, 12, 18, 24]);
xticklabels({'00:00', '06:00', '12:00', '18:00', '24:00'});
xlim([0, 24]);
print('-dpng', '-r400', 'AC_Cluster_Down_Regulation_Stateful.png');
fprintf('空调集群下调潜力图已保存为 AC_Cluster_Down_Regulation_Stateful.png\n');

fprintf('所有 AC 状态化仿真和绘图任务完成。\n');
toc;
