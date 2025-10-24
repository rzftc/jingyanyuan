%% 主程序：预测山东省电动汽车和空调规模至2040年 (使用估算数据 V5.1 - 显式参数传递)
clear; clc; close all;

% --- 预测目标年份 ---
target_years = [2030, 2035, 2040];
predict_end_year = max(target_years);

fprintf('开始预测山东省电动汽车和空调规模 (使用估算数据 V5.1)...\n');
fprintf('目标年份: %s\n', num2str(target_years));
fprintf('预测至: %d\n', predict_end_year);

%% --- 1. 电动汽车规模预测 (Bass 模型 - 参数与 V2/V3/V4/V5 相同) ---
fprintf('\n--- 电动汽车规模预测 (山东省 - 估算) ---\n');
start_hist_year_ev = 2023;
params_ev = struct();
params_ev.start_year = start_hist_year_ev;
params_ev.predict_years = predict_end_year - start_hist_year_ev + 1;
params_ev.N0 = 1.5e6; % 2023年初保有量估算

years_total_ev = params_ev.predict_years;
potential_growth_years_ev = 12; potential_start_ev = 3e6; potential_peak_ev = 6e6;
if years_total_ev <= potential_growth_years_ev
    params_ev.M = round(linspace(potential_start_ev, potential_peak_ev, years_total_ev));
else
    m_growth_ev = round(linspace(potential_start_ev, potential_peak_ev, potential_growth_years_ev));
    m_stable_ev = potential_peak_ev * ones(1, years_total_ev - potential_growth_years_ev);
    params_ev.M = [m_growth_ev, m_stable_ev];
end
fprintf('EV 潜在市场 M (估算): 从 %.1f M/年 增长到 %.1f M/年\n', params_ev.M(1)/1e6, params_ev.M(end)/1e6);

params_ev.p_high = 0.02; params_ev.p_low  = 0.008;
params_ev.q_high = 0.45; params_ev.q_low  = 0.35;
params_ev.r_ac = 0.07;

base_Sur_ev = [1, 0.96, 0.92, 0.88, 0.83, 0.78, 0.72, 0.65, 0.58, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05];
required_len_ev = params_ev.predict_years + length(base_Sur_ev);
extended_Sur_ev = zeros(1, required_len_ev); len_base_ev = length(base_Sur_ev);
extended_Sur_ev(1:len_base_ev) = base_Sur_ev;
params_ev.Sur = extended_Sur_ev;

try
    [ev_high_full, ev_low_full] = ev_scale_prediction(params_ev); % Assume ev_scale works
    predicted_years_ev = params_ev.start_year : (params_ev.start_year + params_ev.predict_years - 1);
    ev_results = table('Size', [length(target_years), 3], ...
                       'VariableTypes', {'int32', 'double', 'double'}, ...
                       'VariableNames', {'Year', 'High_Scenario_EV', 'Low_Scenario_EV'});
    results_idx = 1;
    for year_val = target_years
        idx = find(predicted_years_ev == year_val);
        if ~isempty(idx) && idx <= length(ev_high_full) % Add bounds check
            ev_results.Year(results_idx) = year_val;
            ev_results.High_Scenario_EV(results_idx) = round(ev_high_full(idx));
            ev_results.Low_Scenario_EV(results_idx) = round(ev_low_full(idx));
            results_idx = results_idx + 1;
        else; fprintf('警告: 无法在 EV 预测结果中找到年份 %d\n', year_val); end
    end
    fprintf('电动汽车规模预测结果 (单位：辆):\n'); disp(ev_results);
catch ME_ev
    fprintf('电动汽车规模预测出错: %s\n', ME_ev.message);
     fprintf('发生在文件 %s 的第 %d 行\n', ME_ev.stack(1).file, ME_ev.stack(1).line);
end

%% --- 2. 空调规模预测 (Logistic 模型 - 使用 V5.1 核心逻辑和调整后参数) ---
fprintf('\n--- 空调规模预测 (山东省 - 估算 V5.1) ---\n');

params_ac = struct();
params_ac.method = 'logistic';
params_ac.hist_years = 2018:2022;
params_ac.N_hist = [30e6, 33e6, 36e6, 38e6, 40e6]; % 历史保有量估算

start_predict_year_ac = params_ac.hist_years(end) + 1;
params_ac.predict_years = predict_end_year - start_predict_year_ac + 1;

% [保持 V3/V4 调整] Logistic 模型参数 K, r
params_ac.K_high = 85e6; params_ac.K_low  = 75e6;
params_ac.r_high = 0.08; % 增长率 (略微调高)
params_ac.r_low  = 0.06; % 增长率 (略微调高)

% [保持 V3/V4 调整] 偶然淘汰率 r_ac
params_ac.r_ac = 0.04; % 4%

% [保持 V3/V4 调整] 留存率 Sur (15年) - 确保 Sur(1)=1
base_Sur_ac = [1, 0.98, 0.95, 0.91, 0.86, 0.80, 0.73, 0.65, 0.56, 0.46, 0.36, 0.26, 0.17, 0.10, 0.05]; % Sur(k)=k年后留存
required_len_ac = length(params_ac.N_hist) + params_ac.predict_years + length(base_Sur_ac);
extended_Sur_ac = zeros(1, required_len_ac); len_base_ac = length(base_Sur_ac);
extended_Sur_ac(1:len_base_ac) = base_Sur_ac;
params_ac.Sur = extended_Sur_ac;

try
    % !!! 调用 V5.1 版本的 ac_scale_prediction !!!
    [ac_high_full, ac_low_full] = ac_scale_prediction(params_ac);

    % 提取目标年份的结果 (同前)
    predicted_years_ac = start_predict_year_ac : (start_predict_year_ac + params_ac.predict_years - 1);
    ac_results = table('Size', [length(target_years), 3], ...
                       'VariableTypes', {'int32', 'double', 'double'}, ...
                       'VariableNames', {'Year', 'High_Scenario_AC', 'Low_Scenario_AC'});
    results_idx = 1;
     for year_val = target_years
        idx = find(predicted_years_ac == year_val);
        if ~isempty(idx) && idx <= length(ac_high_full) % 增加长度检查
            ac_results.Year(results_idx) = year_val;
            ac_results.High_Scenario_AC(results_idx) = round(ac_high_full(idx));
            ac_results.Low_Scenario_AC(results_idx) = round(ac_low_full(idx));
            results_idx = results_idx + 1;
        else
             idx_hist = find(params_ac.hist_years == year_val);
             if ~isempty(idx_hist)
                 fprintf('注意: 空调年份 %d 属于历史数据。\n', year_val);
                 ac_results.Year(results_idx) = year_val;
                 ac_results.High_Scenario_AC(results_idx) = round(params_ac.N_hist(idx_hist));
                 ac_results.Low_Scenario_AC(results_idx) = round(params_ac.N_hist(idx_hist));
                 results_idx = results_idx + 1;
             else; fprintf('警告: 无法在 AC 预测结果或历史数据中找到年份 %d\n', year_val); end
        end
    end
    fprintf('空调规模预测结果 (单位：台):\n'); disp(ac_results);
catch ME_ac
    fprintf('空调规模预测出错: %s\n', ME_ac.message);
     fprintf('发生在文件 %s 的第 %d 行\n', ME_ac.stack(1).file, ME_ac.stack(1).line);
end

fprintf('\n预测完成 (请注意使用的是调整后的估算数据 V5.1)。\n');

% --- 合并结果显示 ---
final_results = table(target_years', 'VariableNames', {'Year'});
if exist('ev_results', 'var') && istable(ev_results) && ~isempty(ev_results)
    final_results = outerjoin(final_results, ev_results(:, {'Year', 'High_Scenario_EV', 'Low_Scenario_EV'}), 'Keys', 'Year', 'MergeKeys', true);
else
    final_results.High_Scenario_EV = NaN(height(final_results), 1);
    final_results.Low_Scenario_EV = NaN(height(final_results), 1);
end
if exist('ac_results', 'var') && istable(ac_results) && ~isempty(ac_results)
     final_results = outerjoin(final_results, ac_results(:, {'Year', 'High_Scenario_AC', 'Low_Scenario_AC'}), 'Keys', 'Year', 'MergeKeys', true);
else
    final_results.High_Scenario_AC = NaN(height(final_results), 1);
    final_results.Low_Scenario_AC = NaN(height(final_results), 1);
end

fprintf('\n--- 综合预测结果 (调整后估算 V5.1) ---\n');
disp(final_results);