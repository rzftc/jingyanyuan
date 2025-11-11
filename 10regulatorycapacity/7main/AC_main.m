% main_ac_only.m
clear; close all; clc;

%% 1. System Initialization
rng(2023, 'Threefry'); % Fixed random seed for reproducibility

T_total = 24; % Total simulation duration (hours)
dt = 5/60;    % Time resolution (hours), e.g., 5 minutes
time_points = 0:dt:T_total; % Simulation time points (from 0 to T_total)

base_price = 30; % Base electricity price (yuan/kWh)

%% 2. Initialize AC Parameters
acFile = 'AC_template.xlsx'; % AC data file path (assuming it's in the same directory or MATLAB path)
ACs = initializeACsFromExcel(acFile);
num_AC = length(ACs);

% Store original AC parameters for resetting in each price iteration
% This ensures each incentive scenario starts from the same baseline AC characteristics.
for i = 1:num_AC
    ACs(i).Tset_original = ACs(i).Tset;
    ACs(i).Tmax_original = ACs(i).Tmax;
    ACs(i).Tmin_original = ACs(i).Tmin;
    % Ensure p_incentive exists, if not from Excel, give a random default.
    % (Excel file 'AC_template.xlsx - Sheet1.csv' seems to already have p_incentive)
    if ~isfield(ACs(i), 'p_incentive')
        ACs(i).p_incentive = round(60*rand(), 1);
    end
end

%% 3. Incentive Response Parameters
% These parameters define the range and behavior of temperature setpoint adjustment
% based on incentive electricity prices.
p_min = 15;        % Lower price threshold for down-regulation elasticity (yuan)
p_max = 50;        % Upper price threshold for down-regulation elasticity (yuan)
p_min_prime = 10;  % Lower price threshold for up-regulation elasticity (yuan)
p_max_prime = 40;  % Upper price threshold for up-regulation elasticity (yuan)
T_set_max = 3;     % Maximum allowed temperature setpoint deviation (Celsius)

% Define the full range of potential incentive prices.
p_incentive_range_full = linspace(0, 50, 25);

% --- 修改：定义仅用于绘图和仿真的价格索引 ---
% Define a set of indices for incentive prices to be plotted AND simulated.
% These indices correspond to the 'p_incentive_range_full' array.
plot_price_indices = [1, 7, 13, 19, 25];

% Adjust plot_price_indices if the actual number of potential prices is less than expected.
if length(p_incentive_range_full) < max(plot_price_indices)
    plot_price_indices = round(linspace(1, length(p_incentive_range_full), min(5, length(p_incentive_range_full))));
    plot_price_indices = unique(plot_price_indices); % Ensure unique indices
end
num_prices_to_simulate = length(plot_price_indices);
fprintf('将仅对以下价格索引进行仿真: %s\n', num2str(plot_price_indices));
% --- 修改结束 ---

%% 4. Initialize Data Storage for Aggregated AC Response (Only for selected prices)
% Dimensions: [Number_of_Time_Points x Number_of_Simulated_Prices]
all_AC_Up = zeros(length(time_points), num_prices_to_simulate);
all_AC_Down = zeros(length(time_points), num_prices_to_simulate);

%% 5. Main Loop: Iterate ONLY through selected incentive price indices for simulation
% For each selected incentive price, simulate the entire 24-hour period for all AC units.
% --- 修改：循环仅遍历 plot_price_indices ---
for k = 1:num_prices_to_simulate % Loop counter from 1 to number of prices to simulate
    p_idx = plot_price_indices(k); % Get the actual price index from the selected list
    current_p = p_incentive_range_full(p_idx); % Get the current incentive price for this iteration
    fprintf('\n== Simulating for Selected Incentive Price Index %d (Price: %.1f Yuan, dt = %.2f hours) ==\n', p_idx, current_p, dt);

    % Create a temporary, independent copy of the ACs structure for this specific price scenario.
    temp_ACs_for_price_scenario = ACs; % Starts with original AC parameters.

    % --- Step 5.1: Apply Incentive-Based Temperature Adjustment and Pre-calculate AC parameters ---
    parfor i = 1:num_AC
        % Reset individual AC's settable temperature range to its original values
        temp_ACs_for_price_scenario(i).Tset = temp_ACs_for_price_scenario(i).Tset_original;
        temp_ACs_for_price_scenario(i).Tmax = temp_ACs_for_price_scenario(i).Tmax_original;
        temp_ACs_for_price_scenario(i).Tmin = temp_ACs_for_price_scenario(i).Tmin_original;

        % Calculate participation probability
        participation = calculateParticipation(current_p, base_price);

        % Calculate potential temperature adjustment
        [~, ~, deltaT_flex_magnitude] = incentiveTempAC(...
            current_p, p_min, p_max, p_min_prime, p_max_prime, T_set_max);

        % Decide participation
        temp_ACs_for_price_scenario(i).ptcp = (rand() < participation);

        % Adjust temperature range if participating
        if temp_ACs_for_price_scenario(i).ptcp
            temp_ACs_for_price_scenario(i).Tmax = temp_ACs_for_price_scenario(i).Tset_original + deltaT_flex_magnitude;
            temp_ACs_for_price_scenario(i).Tmin = temp_ACs_for_price_scenario(i).Tset_original - deltaT_flex_magnitude;
            temp_ACs_for_price_scenario(i).Tmax = min(temp_ACs_for_price_scenario(i).Tmax, 28);
            temp_ACs_for_price_scenario(i).Tmin = max(temp_ACs_for_price_scenario(i).Tmin, 18);
        end

        % Generate ambient temperature profile
        base_ambient_temp = temp_ACs_for_price_scenario(i).Tset_original + 4*sin(2*pi*time_points/24);
        actual_temp_range = temp_ACs_for_price_scenario(i).Tmax - temp_ACs_for_price_scenario(i).Tmin;
        noise_factor = 0.2;
        noise = noise_factor * actual_temp_range * randn(size(time_points));
        temp_ACs_for_price_scenario(i).T_ja = min(max(base_ambient_temp + noise, temp_ACs_for_price_scenario(i).Tmin), temp_ACs_for_price_scenario(i).Tmax);

        % Pre-calculate coefficients
        [alpha, beta, gamma] = calculateACABC_single(...
            temp_ACs_for_price_scenario(i).R, temp_ACs_for_price_scenario(i).C, temp_ACs_for_price_scenario(i).eta,...
            temp_ACs_for_price_scenario(i).Tmax, temp_ACs_for_price_scenario(i).Tmin, temp_ACs_for_price_scenario(i).Tset, dt);
        temp_ACs_for_price_scenario(i).alpha = alpha;
        temp_ACs_for_price_scenario(i).beta = beta;
        temp_ACs_for_price_scenario(i).gamma = gamma;
    end

    % Initialize local accumulators
    local_AC_Up_total_per_price = zeros(length(time_points), 1);
    local_AC_Down_total_per_price = zeros(length(time_points), 1);

    % --- Step 5.2: Time-Stepping Simulation ---
    for t_idx = 1:length(time_points)
        current_time_step_total_DeltaP_plus = 0;
        current_time_step_total_DeltaP_minus = 0;

        for i = 1:num_AC
            if temp_ACs_for_price_scenario(i).ptcp
                % Calculate baseline power
                P_base_val = ACbaseP_single(...
                    temp_ACs_for_price_scenario(i).T_ja(t_idx), temp_ACs_for_price_scenario(i).Tset, ...
                    temp_ACs_for_price_scenario(i).R, temp_ACs_for_price_scenario(i).eta);

                % Update SOC
                SOC_val = calculateACS_single(temp_ACs_for_price_scenario(i).T_ja(t_idx),...
                                            temp_ACs_for_price_scenario(i).Tmax, temp_ACs_for_price_scenario(i).Tmin);

                % Calculate adjustment potential
                [DeltaP_plus_t, DeltaP_minus_t] = calculateACAdjustmentPotentia(...
                    P_base_val, 2*abs(P_base_val), 0,...
                    temp_ACs_for_price_scenario(i).alpha, temp_ACs_for_price_scenario(i).beta, temp_ACs_for_price_scenario(i).gamma,...
                    SOC_val, dt);

                % Aggregate potentials
                current_time_step_total_DeltaP_plus = current_time_step_total_DeltaP_plus + DeltaP_plus_t;
                current_time_step_total_DeltaP_minus = current_time_step_total_DeltaP_minus + DeltaP_minus_t;
            end
        end
        local_AC_Up_total_per_price(t_idx) = current_time_step_total_DeltaP_plus;
        local_AC_Down_total_per_price(t_idx) = current_time_step_total_DeltaP_minus;
    end
    % --- 修改：使用循环计数器 k 作为列索引存储结果 ---
    all_AC_Up(:, k) = local_AC_Up_total_per_price;
    all_AC_Down(:, k) = local_AC_Down_total_per_price;
    % --- 修改结束 ---
end
% --- 修改：主仿真循环结束 ---

% --- 新增代码：在绘图前保存数据到Excel ---
excel_filename = 'AC_Plot_Data_Selected.xlsx'; % 文件名稍作区分
fprintf('\n正在将仿真和绘图所需的数据保存到 %s ...\n', excel_filename);

% 创建一个包含时间列的表
dataTable = table(time_points', 'VariableNames', {'Time_Hours'});

% 循环添加选中价格的上调和下调潜力数据
% --- 修改：循环遍历模拟的价格数量 ---
for k = 1:num_prices_to_simulate
    p_idx = plot_price_indices(k); % 获取对应的原始价格索引
    current_price = p_incentive_range_full(p_idx); % 获取实际价格值

    % 生成列名
    up_col_name = matlab.lang.makeValidName(sprintf('Up_Price_%.1f', current_price));
    down_col_name = matlab.lang.makeValidName(sprintf('Down_Price_%.1f', current_price));

    % 添加数据列到表中 (使用循环计数器 k 作为列索引)
    dataTable.(up_col_name) = all_AC_Up(:, k);
    dataTable.(down_col_name) = all_AC_Down(:, k);
end
% --- 修改结束 ---

% 写入Excel文件
try
    writetable(dataTable, excel_filename);
    fprintf('数据成功保存到 %s\n', excel_filename);
catch ME
    fprintf('*** 保存Excel文件时出错: %s ***\n', ME.message);
end
% --- 新增代码结束 ---


%% 6. Plotting Results: Cluster Up- and Down-Regulation Capacity vs. Time for Selected Incentive Prices
% This section generates plots to visualize how the aggregated AC flexibility
% changes over time, and how it is influenced by different incentive prices.

colors = lines(num_prices_to_simulate); % 使用模拟的价格数量

% --- 6.1 Plotting Aggregated AC Up-Regulation Capacity ---
figure('Name', 'AC Cluster Up-Regulation Capacity by Selected Incentive Price', 'Position', [100 100 1000 600]);
hold on; grid on;

legend_entries_up = {};
% --- 修改：循环和索引与Excel保存部分保持一致 ---
for k = 1:num_prices_to_simulate
    p_idx = plot_price_indices(k); % 获取原始索引
    current_price = p_incentive_range_full(p_idx); % 获取价格值

    % Plot (使用循环计数器 k 作为列索引)
    plot(time_points, all_AC_Up(:, k), 'LineWidth', 1.5, 'Color', colors(k,:), ...
         'DisplayName', sprintf('电价: %.1f 分/℃', current_price));

    legend_entries_up{end+1} = sprintf('电价: %.1f 分/℃', current_price);
end
% --- 修改结束 ---
hold off;

xlabel('时间', 'FontSize', 18); % <-- 修改
ylabel('集群上调潜力 (kW)', 'FontSize', 18);
if ~isempty(legend_entries_up)
    legend(legend_entries_up, 'Location', 'best', 'FontSize', 16);
end
set(gca, 'FontSize', 16);

% --- 新增修改：设置X轴刻度为HH:MM格式 ---
xticks([0, 6, 12, 18, 24]);
xticklabels({'00:00', '06:00', '12:00', '18:00', '24:00'});
xlim([0, 24]); % 确保X轴范围与刻度匹配
% --- 修改结束 ---

print('-dpng', '-r400', 'AC_Cluster_Up_Regulation_by_Selected_Price.png');
fprintf('空调集群上调潜力图已保存为 AC_Cluster_Up_Regulation_by_Selected_Price.png\n');

% --- 6.2 Plotting Aggregated AC Down-Regulation Capacity ---
figure('Name', 'AC Cluster Down-Regulation Capacity by Selected Incentive Price', 'Position', [100 750 1000 600]);
hold on; grid on;

legend_entries_down = {};
% --- 修改：循环和索引与Excel保存部分保持一致 ---
for k = 1:num_prices_to_simulate
    p_idx = plot_price_indices(k);
    current_price = p_incentive_range_full(p_idx);

    % Plot (使用循环计数器 k 作为列索引)
    plot(time_points, all_AC_Down(:, k), 'LineWidth', 1.5, 'Color', colors(k,:), ...
         'DisplayName', sprintf('电价: %.1f 分/℃', current_price));

    legend_entries_down{end+1} = sprintf('电价: %.1f 分/℃', current_price);
end
% --- 修改结束 ---
hold off;

xlabel('时间', 'FontSize', 18); % <-- 修改
ylabel('集群下调潜力 (kW)', 'FontSize', 18);
if ~isempty(legend_entries_down)
    legend(legend_entries_down, 'Location', 'best', 'FontSize', 16);
end
set(gca, 'FontSize', 16);

% --- 新增修改：设置X轴刻度为HH:MM格式 ---
xticks([0, 6, 12, 18, 24]);
xticklabels({'00:00', '06:00', '12:00', '18:00', '24:00'});
xlim([0, 24]); % 确保X轴范围与刻度匹配
% --- 修改结束 ---

print('-dpng', '-r400', 'AC_Cluster_Down_Regulation_by_Selected_Price.png');
fprintf('空调集群下调潜力图已保存为 AC_Cluster_Down_Regulation_by_Selected_Price.png\n');

fprintf('All AC-only simulation (selected prices) and plotting tasks completed.\n');