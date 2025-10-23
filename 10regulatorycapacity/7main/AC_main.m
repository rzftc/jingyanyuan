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

% Define the range of incentive prices to simulate and plot.
% This will create 25 distinct incentive price scenarios.
p_incentive_range = linspace(0, 50, 25); 

%% 4. Initialize Data Storage for Aggregated AC Response across all Incentive Prices
% all_AC_Up will store total up-regulation capacity [kW] for each time point and each incentive price.
% all_AC_Down will store total down-regulation capacity [kW] for each time point and each incentive price.
% Dimensions: [Number_of_Time_Points x Number_of_Incentive_Prices]
all_AC_Up = zeros(length(time_points), length(p_incentive_range));
all_AC_Down = zeros(length(time_points), length(p_incentive_range));

%% 5. Main Loop: Iterate through each incentive price scenario
% For each incentive price, simulate the entire 24-hour period for all AC units.
for p_idx = 1:length(p_incentive_range)
    current_p = p_incentive_range(p_idx); % Get the current incentive price for this iteration
    fprintf('\n== Simulating for Incentive Price: %.1f Yuan (dt = %.2f hours) ==\n', current_p, dt);

    % Create a temporary, independent copy of the ACs structure for this specific price scenario.
    % This is crucial to ensure that modifications to AC properties (like Tmax/Tmin, ptcp, alpha/beta/gamma)
    % for one price scenario do not affect subsequent price scenarios.
    % Also, using 'temp_ACs_for_price_scenario' in the parfor loop helps with MATLAB's transparency rules.
    temp_ACs_for_price_scenario = ACs; % Starts with original AC parameters.

    % --- Step 5.1: Apply Incentive-Based Temperature Adjustment and Pre-calculate AC parameters ---
    % This section determines each AC's participation and adjusted comfort range
    % for the current incentive price, and pre-calculates dynamic model coefficients.
    % Using parfor for device-level operations, assuming enough workers are available.
    parfor i = 1:num_AC 
        % Reset individual AC's settable temperature range to its original values
        % before applying the new incentive's effect.
        temp_ACs_for_price_scenario(i).Tset = temp_ACs_for_price_scenario(i).Tset_original;
        temp_ACs_for_price_scenario(i).Tmax = temp_ACs_for_price_scenario(i).Tmax_original;
        temp_ACs_for_price_scenario(i).Tmin = temp_ACs_for_price_scenario(i).Tmin_original;
        
        % Calculate participation probability for the current AC under current_p
        participation = calculateParticipation(current_p, base_price);
        
        % Calculate the potential temperature adjustment range (deltaT_flex_magnitude is the magnitude of change)
        % Note: incentiveTempAC returns T_set_up and T_set_down as positive values representing flexibility.
        % The 'deltaT' output is a random value within [T_set_down, T_set_up].
        [~, ~, deltaT_flex_magnitude] = incentiveTempAC(...
            current_p, p_min, p_max, p_min_prime, p_max_prime, T_set_max);

        % Randomly decide if this AC unit participates based on 'participation' probability.
        % This introduces user behavior uncertainty.
        temp_ACs_for_price_scenario(i).ptcp = (rand() < participation);

        % If the AC unit decides to participate, adjust its effective operating temperature range.
        % The 'deltaT_flex_magnitude' is used to define the half-width of the new comfort band
        % around the original setpoint (Tset_original).
        if temp_ACs_for_price_scenario(i).ptcp
            temp_ACs_for_price_scenario(i).Tmax = temp_ACs_for_price_scenario(i).Tset_original + deltaT_flex_magnitude;
            temp_ACs_for_price_scenario(i).Tmin = temp_ACs_for_price_scenario(i).Tset_original - deltaT_flex_magnitude;
            % Ensure Tmax/Tmin don't exceed hard physical limits if any, or general comfort bounds
            % (e.g., typically Tmax < 28C, Tmin > 16C, depending on context).
            temp_ACs_for_price_scenario(i).Tmax = min(temp_ACs_for_price_scenario(i).Tmax, 28); 
            temp_ACs_for_price_scenario(i).Tmin = max(temp_ACs_for_price_scenario(i).Tmin, 16); 
        end

        % Generate ambient temperature profile (T_ja) for the entire simulation period (all time points).
        % This profile includes a base sinusoidal pattern + random noise, bounded by Tmin/Tmax.
        % Tset_original is used for the base profile's center.
        base_ambient_temp = temp_ACs_for_price_scenario(i).Tset_original + 4*sin(2*pi*time_points/24);
        actual_temp_range = temp_ACs_for_price_scenario(i).Tmax - temp_ACs_for_price_scenario(i).Tmin;
        noise_factor = 0.2; % 20% of the comfort range
        noise = noise_factor * actual_temp_range * randn(size(time_points));
        % Apply new Tmax/Tmin to bound the ambient temperature, simulating control adherence.
        temp_ACs_for_price_scenario(i).T_ja = min(max(base_ambient_temp + noise, temp_ACs_for_price_scenario(i).Tmin), temp_ACs_for_price_scenario(i).Tmax);

        % Pre-calculate constant coefficients (alpha, beta, gamma) for the AC's thermal model.
        % These depend on the AC's physical properties (R, C, eta) and its *adjusted* Tmax/Tmin, Tset, dt.
        % The Tset here is still the original Tset, consistent with `ac_ev_relation_analys_fast.m`'s usage.
        [alpha, beta, gamma] = calculateACABC_single(...
            temp_ACs_for_price_scenario(i).R, temp_ACs_for_price_scenario(i).C, temp_ACs_for_price_scenario(i).eta,...
            temp_ACs_for_price_scenario(i).Tmax, temp_ACs_for_price_scenario(i).Tmin, temp_ACs_for_price_scenario(i).Tset, dt); 
        temp_ACs_for_price_scenario(i).alpha = alpha;
        temp_ACs_for_price_scenario(i).beta = beta;
        temp_ACs_for_price_scenario(i).gamma = gamma;
    end

    % Initialize local accumulators for total aggregated power for this incentive price.
    local_AC_Up_total_per_price = zeros(length(time_points), 1);
    local_AC_Down_total_per_price = zeros(length(time_points), 1);

    % --- Step 5.2: Time-Stepping Simulation for the current Incentive Price ---
    % This loop progresses through time, calculating aggregated power for all ACs at each time point.
    for t_idx = 1:length(time_points)
        % For progress tracking in console (optional, commented out for cleaner output by default)
        % if mod(t_idx-1, round(1/dt)) == 0 && t_idx > 1
        %      fprintf('  Processing time %.1f hours...\n', time_points(t_idx));
        % end

        current_time_step_total_DeltaP_plus = 0;
        current_time_step_total_DeltaP_minus = 0;

        % Calculate individual AC properties and sum them for the current time step.
        for i = 1:num_AC
            if temp_ACs_for_price_scenario(i).ptcp % Only consider participating ACs
                % Calculate baseline power for the current AC at the current time point.
                % Uses the adjusted Tset and the current ambient temperature (T_ja).
                P_base_val = ACbaseP_single(...
                    temp_ACs_for_price_scenario(i).T_ja(t_idx), temp_ACs_for_price_scenario(i).Tset, ... 
                    temp_ACs_for_price_scenario(i).R, temp_ACs_for_price_scenario(i).eta);

                % Update SOC (State of Charge) based on the current ambient temperature and adjusted Tmax/Tmin.
                % The SOC is a normalized representation of how "full" the AC's thermal storage is
                % within its current comfort band.
                SOC_val = calculateACS_single(temp_ACs_for_price_scenario(i).T_ja(t_idx),... 
                                             temp_ACs_for_price_scenario(i).Tmax, temp_ACs_for_price_scenario(i).Tmin);

                % Calculate the adjustment potential (up and down) for this individual AC.
                % Pmax is approximated as 2 * absolute_baseline_power for cooling (ability to cool more).
                % Pmin is 0 (ability to stop cooling).
                % dt (last parameter) is the time step for potential calculation.
                [DeltaP_plus_t, DeltaP_minus_t] = calculateACAdjustmentPotentia(...
                    P_base_val, 2*abs(P_base_val), 0,... 
                    temp_ACs_for_price_scenario(i).alpha, temp_ACs_for_price_scenario(i).beta, temp_ACs_for_price_scenario(i).gamma,...
                    SOC_val, dt); 

                % Aggregate the individual AC's potentials to the total cluster potential for this time step.
                current_time_step_total_DeltaP_plus = current_time_step_total_DeltaP_plus + DeltaP_plus_t;
                current_time_step_total_DeltaP_minus = current_time_step_total_DeltaP_minus + DeltaP_minus_t;
            end
        end
        % Store the aggregated potentials for the current time step and incentive price.
        local_AC_Up_total_per_price(t_idx) = current_time_step_total_DeltaP_plus;
        local_AC_Down_total_per_price(t_idx) = current_time_step_total_DeltaP_minus;
    end
    % Store the results for the current incentive price scenario into the main result arrays.
    all_AC_Up(:, p_idx) = local_AC_Up_total_per_price;
    all_AC_Down(:, p_idx) = local_AC_Down_total_per_price;
end

%% 6. Plotting Results: Cluster Up- and Down-Regulation Capacity vs. Time for Selected Incentive Prices
% This section generates plots to visualize how the aggregated AC flexibility
% changes over time, and how it is influenced by different incentive prices.

% Define a set of indices for incentive prices to be plotted.
% These indices correspond to the 'p_incentive_range' array.
% Choosing 5 representative prices across the range (e.g., low, medium-low, mid, medium-high, high incentive).
% Example indices for a 25-point p_incentive_range: 1 (min price), 7, 13 (mid price), 19, 25 (max price).
plot_price_indices = [1, 7, 13, 19, 25]; 

% Adjust plot_price_indices if the actual number of simulated prices is less than expected.
if length(p_incentive_range) < max(plot_price_indices)
    plot_price_indices = round(linspace(1, length(p_incentive_range), min(5, length(p_incentive_range))));
    plot_price_indices = unique(plot_price_indices); % Ensure unique indices
end

colors = lines(length(plot_price_indices)); % Get distinct colors for each plotted line

% --- 6.1 Plotting Aggregated AC Up-Regulation Capacity ---
figure('Name', 'AC Cluster Up-Regulation Capacity by Incentive Price', 'Position', [100 100 1000 600]);
hold on; % Keep the plot active to add multiple lines
grid on; % Add a grid for readability

legend_entries_up = {}; % Initialize cell array for legend labels
for k = 1:length(plot_price_indices)
    p_idx = plot_price_indices(k); % Current index from the selected set
    current_price = p_incentive_range(p_idx); % Get the actual incentive price
    
    % Plot the up-regulation capacity time series for this price
    plot(time_points, all_AC_Up(:, p_idx), 'LineWidth', 1.5, 'Color', colors(k,:), ...
         'DisplayName', sprintf('电价: %.1f 分/℃', current_price));
    
    legend_entries_up{end+1} = sprintf('电价: %.1f 分/℃', current_price);
end
hold off; % Release the plot

xlabel('时间 (小时)', 'FontSize', 16);
ylabel('集群上调潜力 (kW)', 'FontSize', 16);

% Add legend
if ~isempty(legend_entries_up)
    legend(legend_entries_up, 'Location', 'best', 'FontSize', 14);
end
set(gca, 'FontSize', 14); % Set font size for axes ticks

% Save the plot to a PNG file
print('-dpng', '-r400', 'AC_Cluster_Up_Regulation_by_Price.png');
fprintf('空调集群上调潜力图已保存为 AC_Cluster_Up_Regulation_by_Price.png\n');

% --- 6.2 Plotting Aggregated AC Down-Regulation Capacity ---
figure('Name', 'AC Cluster Down-Regulation Capacity by Incentive Price', 'Position', [100 750 1000 600]);
hold on; 
grid on;

legend_entries_down = {}; % Initialize cell array for legend labels
for k = 1:length(plot_price_indices)
    p_idx = plot_price_indices(k);
    current_price = p_incentive_range(p_idx);
    
    % Plot the down-regulation capacity time series for this price
    plot(time_points, all_AC_Down(:, p_idx), 'LineWidth', 1.5, 'Color', colors(k,:), ...
         'DisplayName', sprintf('电价: %.1f 分/℃', current_price));
    
    legend_entries_down{end+1} = sprintf('电价: %.1f 分/℃', current_price);
end
hold off;

xlabel('时间 (小时)', 'FontSize', 16);
ylabel('集群下调潜力 (kW)', 'FontSize', 16);

% Add legend
if ~isempty(legend_entries_down)
    legend(legend_entries_down, 'Location', 'best', 'FontSize', 14);
end
set(gca, 'FontSize', 14);

% Save the plot to a PNG file
print('-dpng', '-r400', 'AC_Cluster_Down_Regulation_by_Price.png');
fprintf('空调集群下调潜力图已保存为 AC_Cluster_Down_Regulation_by_Price.png\n');

fprintf('All AC-only simulation and plotting tasks completed.\n');