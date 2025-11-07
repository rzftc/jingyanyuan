%% 预测特定 220kV 节点的 EV 和 AC 规模 (降尺度方法)
clear; clc; close all;

% =======================================================================
% TODO: 用户必须根据 220kV 节点的实际数据修改以下分配因子
% =======================================================================
% 示例：假设此 220kV 节点...
% 1. 覆盖了山东省 1.5% 的空调保有量/负荷/家庭户数
% 2. 覆盖了山东省 1.0% 的电动汽车保有量/人口/机动车数量
%
% *** 请将这两个值修改为您节点的真实比例 (0.0 到 1.0 之间) ***
node_ac_allocation_factor = 0.015;  % (1.5%) - 请修改此值
node_ev_allocation_factor = 0.010;  % (1.0%) - 请修改此值
% =======================================================================

fprintf('开始使用节点分配比例进行预测:\n');
fprintf('  - AC 分配因子: %.2f%%\n', node_ac_allocation_factor * 100);
fprintf('  - EV 分配因子: %.2f%%\n', node_ev_allocation_factor * 100);

% --- 预测目标年份 ---
target_years = [2030, 2035, 2040];
predict_end_year = max(target_years);

fprintf('预测至: %d\n', predict_end_year);

%% --- 1. 电动汽车省级规模预测 (Bass 模型) ---
% (此部分逻辑与 test2.m 相同，但变量名增加了 _province 后缀)
fprintf('\n--- 步骤 1: 正在计算【省级】EV 预测... ---\n');
start_hist_year_ev = 2023;
params_ev_province = struct();
params_ev_province.start_year = start_hist_year_ev;
params_ev_province.predict_years = predict_end_year - start_hist_year_ev + 1;
params_ev_province.N0 = 1.5e6; % 2023年省级初始保有量估算

years_total_ev = params_ev_province.predict_years;
potential_growth_years_ev = 12; potential_start_ev = 3e6; potential_peak_ev = 6e6;
if years_total_ev <= potential_growth_years_ev
    params_ev_province.M = round(linspace(potential_start_ev, potential_peak_ev, years_total_ev));
else
    m_growth_ev = round(linspace(potential_start_ev, potential_peak_ev, potential_growth_years_ev));
    m_stable_ev = potential_peak_ev * ones(1, years_total_ev - potential_growth_years_ev);
    params_ev_province.M = [m_growth_ev, m_stable_ev];
end
fprintf('省级 EV 潜在市场 M (估算): 从 %.1f M/年 增长到 %.1f M/年\n', params_ev_province.M(1)/1e6, params_ev_province.M(end)/1e6);

params_ev_province.p_high = 0.02; params_ev_province.p_low  = 0.008;
params_ev_province.q_high = 0.45; params_ev_province.q_low  = 0.35;
params_ev_province.r_ac = 0.07;

base_Sur_ev = [1, 0.96, 0.92, 0.88, 0.83, 0.78, 0.72, 0.65, 0.58, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05];
required_len_ev = params_ev_province.predict_years + length(base_Sur_ev);
extended_Sur_ev = zeros(1, required_len_ev); len_base_ev = length(base_Sur_ev);
extended_Sur_ev(1:len_base_ev) = base_Sur_ev;
params_ev_province.Sur = extended_Sur_ev;

try
    [ev_high_full_province, ev_low_full_province] = ev_scale_prediction(params_ev_province);
    predicted_years_ev = params_ev_province.start_year : (params_ev_province.start_year + params_ev_province.predict_years - 1);
catch ME_ev
    fprintf('省级电动汽车规模预测出错: %s\n', ME_ev.message);
    ev_high_full_province = nan(1, params_ev_province.predict_years);
    ev_low_full_province = nan(1, params_ev_province.predict_years);
    predicted_years_ev = params_ev_province.start_year : (params_ev_province.start_year + params_ev_province.predict_years - 1);
end

%% --- 2. 空调省级规模预测 (Logistic 模型) ---
% (此部分逻辑与 test2.m 相同，但变量名增加了 _province 后缀)
fprintf('\n--- 步骤 2: 正在计算【省级】AC 预测... ---\n');

params_ac_province = struct();
params_ac_province.method = 'logistic';
params_ac_province.hist_years = 2018:2022;
params_ac_province.N_hist = [30e6, 33e6, 36e6, 38e6, 40e6]; % 省级历史保有量估算

start_predict_year_ac = params_ac_province.hist_years(end) + 1;
params_ac_province.predict_years = predict_end_year - start_predict_year_ac + 1;

params_ac_province.K_high = 85e6; params_ac_province.K_low  = 75e6;
params_ac_province.r_high = 0.08; 
params_ac_province.r_low  = 0.06; 
params_ac_province.r_ac = 0.04; % 偶然淘汰率

base_Sur_ac = [1, 0.98, 0.95, 0.91, 0.86, 0.80, 0.73, 0.65, 0.56, 0.46, 0.36, 0.26, 0.17, 0.10, 0.05]; 
required_len_ac = length(params_ac_province.N_hist) + params_ac_province.predict_years + length(base_Sur_ac);
extended_Sur_ac = zeros(1, required_len_ac); len_base_ac = length(base_Sur_ac);
extended_Sur_ac(1:len_base_ac) = base_Sur_ac;
params_ac_province.Sur = extended_Sur_ac;

try
    [ac_high_full_province, ac_low_full_province] = ac_scale_prediction(params_ac_province);
    predicted_years_ac = start_predict_year_ac : (start_predict_year_ac + params_ac_province.predict_years - 1);
catch ME_ac
    fprintf('省级空调规模预测出错: %s\n', ME_ac.message);
    ac_high_full_province = nan(1, params_ac_province.predict_years);
    ac_low_full_province = nan(1, params_ac_province.predict_years);
    predicted_years_ac = start_predict_year_ac : (start_predict_year_ac + params_ac_province.predict_years - 1);
end

%% --- 3. 计算【节点级】预测结果 ---
fprintf('\n--- 步骤 3: 正在应用分配因子计算【节点级】结果... ---\n');

% --- 3a. EV 节点预测 ---
ev_n0_node = params_ev_province.N0 * node_ev_allocation_factor;
ev_high_full_node = ev_high_full_province * node_ev_allocation_factor;
ev_low_full_node  = ev_low_full_province * node_ev_allocation_factor;

% --- 3b. AC 节点预测 ---
% 历史数据也需要按比例缩小
ac_hist_node = params_ac_province.N_hist * node_ac_allocation_factor;
% 预测数据按比例缩小
ac_high_full_node = ac_high_full_province * node_ac_allocation_factor;
ac_low_full_node  = ac_low_full_province * node_ac_allocation_factor;


%% --- 4. 提取并显示目标年份的【节点级】结果 ---
node_results = table(target_years', 'VariableNames', {'Year'});

% --- 填充 EV 节点结果 ---
ev_high_node_target = zeros(size(target_years));
ev_low_node_target = zeros(size(target_years));
for i = 1:length(target_years)
    year_val = target_years(i);
    idx = find(predicted_years_ev == year_val);
    if ~isempty(idx) && idx <= length(ev_high_full_node)
        ev_high_node_target(i) = ev_high_full_node(idx);
        ev_low_node_target(i) = ev_low_full_node(idx);
    else
        fprintf('警告: 无法在 EV 预测结果中找到年份 %d\n', year_val);
    end
end
node_results.High_Scenario_EV_Node = round(ev_high_node_target');
node_results.Low_Scenario_EV_Node = round(ev_low_node_target');

% --- 填充 AC 节点结果 ---
ac_high_node_target = zeros(size(target_years));
ac_low_node_target = zeros(size(target_years));
for i = 1:length(target_years)
    year_val = target_years(i);
    idx_pred = find(predicted_years_ac == year_val); % 查找预测年份
    idx_hist = find(params_ac_province.hist_years == year_val); % 查找历史年份
    
    if ~isempty(idx_pred) && idx_pred <= length(ac_high_full_node)
        ac_high_node_target(i) = ac_high_full_node(idx_pred);
        ac_low_node_target(i) = ac_low_full_node(idx_pred);
    elseif ~isempty(idx_hist) && idx_hist <= length(ac_hist_node)
        % 如果目标年份在历史数据中
        ac_high_node_target(i) = ac_hist_node(idx_hist);
        ac_low_node_target(i) = ac_hist_node(idx_hist);
    else
        fprintf('警告: 无法在 AC 预测或历史结果中找到年份 %d\n', year_val);
    end
end
node_results.High_Scenario_AC_Node = round(ac_high_node_target');
node_results.Low_Scenario_AC_Node = round(ac_low_node_target');

% --- 显示最终表格 ---
fprintf('\n======================================================\n');
fprintf('--- 特定 220kV 节点预测结果 (估算) ---\n');
disp(node_results);
fprintf('======================================================\n');


%% --- 步骤 6: 多模型节点级预测对比 (集成随机森林) ---
fprintf('\n--- 步骤 6: 使用多种新模型进行节点级对比预测 ---\n');

% --- 6.0 准备历史数据和预测时间轴 ---
% AC 历史 (2018-2022)
ac_hist_node = params_ac_province.N_hist * node_ac_allocation_factor;
ac_hist_years = params_ac_province.hist_years;
years_to_predict_ac = predict_end_year - ac_hist_years(end);
future_years_axis_ac = (ac_hist_years(end)+1) : predict_end_year;

% EV 历史 (估算 2019-2023)
if ~isfield(params_ev_province, 'N_hist')
    ev_hist_prov_estimated = [0.5e6, 0.7e6, 0.95e6, 1.2e6, params_ev_province.N0]; 
    ev_hist_years = (params_ev_province.start_year - length(ev_hist_prov_estimated) + 1) : params_ev_province.start_year;
else
    ev_hist_prov_estimated = params_ev_province.N_hist;
    ev_hist_years = (params_ev_province.start_year - length(ev_hist_prov_estimated)) : (params_ev_province.start_year - 1);
end
ev_hist_node = ev_hist_prov_estimated * node_ev_allocation_factor;
years_to_predict_ev = predict_end_year - ev_hist_years(end);
future_years_axis_ev = (ev_hist_years(end)+1) : predict_end_year;


fprintf('节点级历史数据 (估算):\n');
fprintf('  AC (万台, %d-%d): %s\n', ac_hist_years(1), ac_hist_years(end), num2str(round(ac_hist_node/1e4, 2)));
fprintf('  EV (万辆, %d-%d): %s\n', ev_hist_years(1), ev_hist_years(end), num2str(round(ev_hist_node/1e4, 2)));

% --- 6.1 应用灰色预测 GM(1,1) ---
fprintf('  > 正在执行灰色预测 GM(1,1)...\n');
ac_pred_grey = grey_prediction_gm11(ac_hist_node, years_to_predict_ac);
ev_pred_grey = grey_prediction_gm11(ev_hist_node, years_to_predict_ev);

% --- 6.2 应用 ARIMA 模型 ---
fprintf('  > 正在执行 ARIMA 预测...\n');
ac_pred_arima = arima_prediction(ac_hist_node, years_to_predict_ac, [1, 1, 0]);
ev_pred_arima = arima_prediction(ev_hist_node, years_to_predict_ev, [1, 2, 0]); 

% --- 6.3 应用动态 Gompertz 模型 ---
fprintf('  > 正在执行动态 Gompertz 预测...\n');
K_ac_dynamic = linspace(ac_hist_node(end)*1.5, ac_hist_node(end)*2.5, length(ac_hist_node) + years_to_predict_ac);
ac_pred_gompertz = gompertz_dynamic_prediction(ac_hist_node, years_to_predict_ac, K_ac_dynamic);
K_ev_dynamic = linspace(ev_hist_node(end)*5, ev_hist_node(end)*15, length(ev_hist_node) + years_to_predict_ev);
ev_pred_gompertz = gompertz_dynamic_prediction(ev_hist_node, years_to_predict_ev, K_ev_dynamic);

% --- 6.4 [修改] 应用自回归随机森林 (Random Forest Autoregressive) ---
fprintf('  > 正在执行自回归随机森林预测 (算例)...\n');
% 为算例创建与多元回归兼容的(虚拟)影响因子
% AC (历史: 2018-2022, 预测: 2023-2040)
n_hist_ac = length(ac_hist_node);
n_pred_ac = years_to_predict_ac; % 2040-2022 = 18
% 历史特征需要与N_hist对齐 (5个点)
X_hist_ac = [linspace(8, 12, n_hist_ac)', [0; 0; 1; 1; 2]]; % 假设收入+政策 (5x2)
% 未来特征需要与预测年份对齐 (18个点)
X_pred_ac = [linspace(12.5, 25, n_pred_ac)', [repmat(2, 7, 1); repmat(3, n_pred_ac-7, 1)]]; % (18x2)
ac_pred_rf_ar = random_forest_prediction(ac_hist_node, X_hist_ac, X_pred_ac, 50); % 50 棵树

% EV (历史: 2019-2023, 预测: 2024-2040)
n_hist_ev = length(ev_hist_node);
n_pred_ev = years_to_predict_ev; % 2040-2023 = 17
% 历史特征 (5x2)
X_hist_ev = [linspace(9, 13, n_hist_ev)', [100; 150; 250; 400; 600]]; % 假设收入+充电桩
% 未来特征 (17x2)
X_pred_ev = [linspace(13.5, 25, n_pred_ev)', round(linspace(800, 5000, n_pred_ev))'];
ev_pred_rf_ar = random_forest_prediction(ev_hist_node, X_hist_ev, X_pred_ev, 50);

% --- 6.5 汇总对比结果 (提取目标年份 2030, 2035, 2040) ---
res_compare_ac = table(target_years', 'VariableNames', {'Year'});
res_compare_ev = table(target_years', 'VariableNames', {'Year'});

% 提取目标年份索引
[~, idx_target_ac] = ismember(target_years, future_years_axis_ac);
[~, idx_target_ev] = ismember(target_years, future_years_axis_ev);

% 填充 AC 对比表 (单位: 万台)
res_compare_ac.Base_High     = round(ac_high_node_target'/1e4, 1);
if ~all(isnan(ac_pred_grey)); res_compare_ac.Grey_GM11 = round(ac_pred_grey(idx_target_ac)'/1e4, 1); else; res_compare_ac.Grey_GM11 = NaN; end
if ~all(isnan(ac_pred_arima)); res_compare_ac.ARIMA = round(ac_pred_arima(idx_target_ac)'/1e4, 1); else; res_compare_ac.ARIMA = NaN; end
if ~all(isnan(ac_pred_gompertz)); res_compare_ac.Gompertz = round(ac_pred_gompertz(idx_target_ac)'/1e4, 1); else; res_compare_ac.Gompertz = NaN; end
if ~all(isnan(ac_pred_rf_ar)); res_compare_ac.RandForest_AR = round(ac_pred_rf_ar(idx_target_ac)/1e4, 1); else; res_compare_ac.RandForest_AR = NaN; end % 修正

% 填充 EV 对比表 (单位: 万辆)
res_compare_ev.Base_High     = round(ev_high_node_target'/1e4, 1);
if ~all(isnan(ev_pred_grey)); res_compare_ev.Grey_GM11 = round(ev_pred_grey(idx_target_ev)'/1e4, 1); else; res_compare_ev.Grey_GM11 = NaN; end
if ~all(isnan(ev_pred_arima)); res_compare_ev.ARIMA = round(ev_pred_arima(idx_target_ev)'/1e4, 1); else; res_compare_ev.ARIMA = NaN; end
if ~all(isnan(ev_pred_gompertz)); res_compare_ev.Gompertz = round(ev_pred_gompertz(idx_target_ev)'/1e4, 1); else; res_compare_ev.Gompertz = NaN; end
if ~all(isnan(ev_pred_rf_ar)); res_compare_ev.RandForest_AR = round(ev_pred_rf_ar(idx_target_ev)/1e4, 1); else; res_compare_ev.RandForest_AR = NaN; end % 修正

fprintf('\n=== 节点级 AC 规模多模型预测对比 (单位: 万台) ===\n');
disp(res_compare_ac);
fprintf('\n=== 节点级 EV 规模多模型预测对比 (单位: 万辆) ===\n');
disp(res_compare_ev);

% 绘图对比
figure('Name', '多模型节点级预测对比', 'Color', 'w', 'Position', [100, 100, 1000, 400]);
subplot(1,2,1); hold on; title('节点 AC 预测对比');
plot(target_years, res_compare_ac.Base_High, 'k-o', 'LineWidth', 2, 'DisplayName', '原方法(High)');
plot(target_years, res_compare_ac.Grey_GM11, 'b--s', 'LineWidth', 1.5, 'DisplayName', 'GM(1,1)');
plot(target_years, res_compare_ac.ARIMA, 'g-.^', 'LineWidth', 1.5, 'DisplayName', 'ARIMA');
plot(target_years, res_compare_ac.Gompertz, 'm:x', 'LineWidth', 2, 'DisplayName', 'Dyn-Gompertz');
plot(target_years, res_compare_ac.RandForest_AR, 'c-d', 'LineWidth', 1.5, 'DisplayName', 'RF(自回归)'); % 修改
legend('Location', 'best'); xlabel('年份'); ylabel('万台'); grid on; xlim([target_years(1)-2, target_years(end)+2]);

subplot(1,2,2); hold on; title('节点 EV 预测对比');
plot(target_years, res_compare_ev.Base_High, 'k-o', 'LineWidth', 2, 'DisplayName', '原方法(High)');
plot(target_years, res_compare_ev.Grey_GM11, 'b--s', 'LineWidth', 1.5, 'DisplayName', 'GM(1,1)');
plot(target_years, res_compare_ev.ARIMA, 'g-.^', 'LineWidth', 1.5, 'DisplayName', 'ARIMA');
plot(target_years, res_compare_ev.Gompertz, 'm:x', 'LineWidth', 2, 'DisplayName', 'Dyn-Gompertz');
plot(target_years, res_compare_ev.RandForest_AR, 'c-d', 'LineWidth', 1.5, 'DisplayName', 'RF(自回归)'); % 修改
legend('Location', 'best'); xlabel('年份'); ylabel('万辆'); grid on; xlim([target_years(1)-2, target_years(end)+2]);

fprintf('\n新增模型算例执行完毕。\n');
% --- 结束新增步骤 ---



%% --- 5. 绘制 "Nature" 风格【省级】和【节点级】预测曲线 (中文/放大/分图/无标题/高DPI版) ---
% (*** 此部分已按您的最新要求修改 ***)

fprintf('\n--- 步骤 5: 正在生成四张高DPI预测曲线 (中文版)... ---\n');

% --- 5a. 准备 AC 省级绘图数据 (合并历史与预测) ---
ac_plot_years_prov = [params_ac_province.hist_years, predicted_years_ac];
ac_plot_low_prov   = [params_ac_province.N_hist, ac_low_full_province];
ac_plot_high_prov  = [params_ac_province.N_hist, ac_high_full_province];
ac_plot_mean_prov  = (ac_plot_low_prov + ac_plot_high_prov) / 2;
ac_fill_x_prov = [ac_plot_years_prov, fliplr(ac_plot_years_prov)];
ac_fill_y_prov = [ac_plot_low_prov, fliplr(ac_plot_high_prov)];

% --- 5b. 准备 EV 省级绘图数据 (合并初始N0与预测) ---
ev_plot_years_prov = [params_ev_province.start_year - 1, predicted_years_ev];
ev_plot_low_prov   = [params_ev_province.N0, ev_low_full_province];
ev_plot_high_prov  = [params_ev_province.N0, ev_high_full_province];
ev_plot_mean_prov  = (ev_plot_low_prov + ev_plot_high_prov) / 2;
ev_fill_x_prov = [ev_plot_years_prov, fliplr(ev_plot_years_prov)];
ev_fill_y_prov = [ev_plot_low_prov, fliplr(ev_plot_high_prov)];

% --- 5c. 准备 AC 节点绘图数据 ---
ac_plot_years_node = [params_ac_province.hist_years, predicted_years_ac];
ac_plot_low_node   = [ac_hist_node, ac_low_full_node];
ac_plot_high_node  = [ac_hist_node, ac_high_full_node];
ac_plot_mean_node  = (ac_plot_low_node + ac_plot_high_node) / 2;
ac_fill_x_node = [ac_plot_years_node, fliplr(ac_plot_years_node)];
ac_fill_y_node = [ac_plot_low_node, fliplr(ac_plot_high_node)];

% --- 5d. 准备 EV 节点绘图数据 ---
ev_plot_years_node = [params_ev_province.start_year - 1, predicted_years_ev];
ev_plot_low_node   = [ev_n0_node, ev_low_full_node];
ev_plot_high_node  = [ev_n0_node, ev_high_full_node];
ev_plot_mean_node  = (ev_plot_low_node + ev_plot_high_node) / 2;
ev_fill_x_node = [ev_plot_years_node, fliplr(ev_plot_years_node)];
ev_fill_y_node = [ev_plot_low_node, fliplr(ev_plot_high_node)];

% --- 5e. 中文字体设置 ---
chinese_font = 'SimHei'; % 尝试使用 'SimHei' (黑体)
try
    set(0, 'DefaultAxesFontName', chinese_font);
    set(0, 'DefaultTextFontName', chinese_font);
catch
    warning('未能设置中文字体 "%s"，图像中的中文可能显示为方框。', chinese_font);
end

% --- 5f. 绘制【图 1: 省级空调 (AC)】 ---
figureHandle_AC_Prov = figure('Color', 'white', 'Position', [100, 100, 800, 450]);
ax1 = axes(figureHandle_AC_Prov);
hold(ax1, 'on');
fill(ax1, ac_fill_x_prov, ac_fill_y_prov / 1e6, [0.6 0.8 1.0], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.6, 'DisplayName', '高低情景范围');
plot(ax1, ac_plot_years_prov, ac_plot_mean_prov / 1e6, 'Color', [0.2 0.4 0.8], 'LineWidth', 2, ...
     'DisplayName', '均值预测');
ylabel(ax1, '空调规模 (百万台)', 'FontSize', 12);
xlabel(ax1, '年份', 'FontSize', 12);
% title(ax1, '山东省空调 (AC) 规模预测', 'FontSize', 14, 'FontWeight', 'bold'); % <<< 标题已移除
set(ax1, 'FontName', chinese_font, 'FontSize', 12, 'Box', 'off', 'TickDir', 'out');
grid(ax1, 'on');
legend(ax1, 'Location', 'northwest', 'Box', 'off', 'FontSize', 12);
set(ax1, 'XLim', [min(ac_plot_years_prov), max(ac_plot_years_prov)]);
hold(ax1, 'off');

% --- 5g. 绘制【图 2: 省级电动汽车 (EV)】 ---
figureHandle_EV_Prov = figure('Color', 'white', 'Position', [150, 150, 800, 450]);
ax2 = axes(figureHandle_EV_Prov);
hold(ax2, 'on');
fill(ax2, ev_fill_x_prov, ev_fill_y_prov / 1e6, [1.0 0.8 0.6], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.6, 'DisplayName', '高低情景范围');
plot(ax2, ev_plot_years_prov, ev_plot_mean_prov / 1e6, 'Color', [0.8 0.4 0.2], 'LineWidth', 2, ...
     'DisplayName', '均值预测');
ylabel(ax2, '电动汽车规模 (百万台)', 'FontSize', 12); 
xlabel(ax2, '年份', 'FontSize', 12);
% title(ax2, '山东省电动汽车 (EV) 规模预测', 'FontSize', 14, 'FontWeight', 'bold'); % <<< 标题已移除
set(ax2, 'FontName', chinese_font, 'FontSize', 12, 'Box', 'off', 'TickDir', 'out');
grid(ax2, 'on');
legend(ax2, 'Location', 'northwest', 'Box', 'off', 'FontSize', 12);
set(ax2, 'XLim', [min(ac_plot_years_prov), max(ac_plot_years_prov)]);
hold(ax2, 'off');

% --- 5h. 绘制【图 3: 节点空调 (AC)】 ---
figureHandle_AC_Node = figure('Color', 'white', 'Position', [200, 200, 800, 450]);
ax3 = axes(figureHandle_AC_Node);
hold(ax3, 'on');
% 绘制填充区域 (单位：万台)
fill(ax3, ac_fill_x_node, ac_fill_y_node / 1e4, [0.6 0.8 1.0], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.6, 'DisplayName', '高低情景范围');
% 绘制均值线
plot(ax3, ac_plot_years_node, ac_plot_mean_node / 1e4, 'Color', [0.2 0.4 0.8], 'LineWidth', 2, ...
     'DisplayName', '均值预测');
% 风格设置
ylabel(ax3, '空调规模 (万台)', 'FontSize', 12);
xlabel(ax3, '年份', 'FontSize', 12);
% title(ax3, '220kV 节点空调 (AC) 规模预测', 'FontSize', 14, 'FontWeight', 'bold'); % <<< 标题已移除
set(ax3, 'FontName', chinese_font, 'FontSize', 12, 'Box', 'off', 'TickDir', 'out');
grid(ax3, 'on');
legend(ax3, 'Location', 'northwest', 'Box', 'off', 'FontSize', 12);
set(ax3, 'XLim', [min(ac_plot_years_node), max(ac_plot_years_node)]);
hold(ax3, 'off');

% --- 5i. 绘制【图 4: 节点电动汽车 (EV)】 ---
figureHandle_EV_Node = figure('Color', 'white', 'Position', [250, 250, 800, 450]);
ax4 = axes(figureHandle_EV_Node);
hold(ax4, 'on');
% 绘制填充区域 (单位：万台)
fill(ax4, ev_fill_x_node, ev_fill_y_node / 1e4, [1.0 0.8 0.6], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.6, 'DisplayName', '高低情景范围');
% 绘制均值线
plot(ax4, ev_plot_years_node, ev_plot_mean_node / 1e4, 'Color', [0.8 0.4 0.2], 'LineWidth', 2, ...
     'DisplayName', '均值预测');
% 风格设置
ylabel(ax4, '电动汽车规模 (万台)', 'FontSize', 12);
xlabel(ax4, '年份', 'FontSize', 12);
% title(ax4, '220kV 节点电动汽车 (EV) 规模预测', 'FontSize', 14, 'FontWeight', 'bold'); % <<< 标题已移除
set(ax4, 'FontName', chinese_font, 'FontSize', 12, 'Box', 'off', 'TickDir', 'out');
grid(ax4, 'on');
legend(ax4, 'Location', 'northwest', 'Box', 'off', 'FontSize', 12);
set(ax4, 'XLim', [min(ac_plot_years_node), max(ac_plot_years_node)]);
hold(ax4, 'off');

% --- 5j. 保存高DPI图像 ---
try
    dpi = 300; % 设置高分辨率 (300 DPI)
    
    output_filename_ac_prov = 'province_ac_scale_high_dpi.png';
    print(figureHandle_AC_Prov, output_filename_ac_prov, '-dpng', sprintf('-r%d', dpi));
    fprintf('【省级】空调高DPI预测曲线图已保存为: %s\n', output_filename_ac_prov);
    
    output_filename_ev_prov = 'province_ev_scale_high_dpi.png';
    print(figureHandle_EV_Prov, output_filename_ev_prov, '-dpng', sprintf('-r%d', dpi));
    fprintf('【省级】EV高DPI预测曲线图已保存为: %s\n', output_filename_ev_prov);
    
    output_filename_ac_node = 'node_ac_scale_high_dpi.png';
    print(figureHandle_AC_Node, output_filename_ac_node, '-dpng', sprintf('-r%d', dpi));
    fprintf('【节点级】空调高DPI预测曲线图已保存为: %s\n', output_filename_ac_node);
    
    output_filename_ev_node = 'node_ev_scale_high_dpi.png';
    print(figureHandle_EV_Node, output_filename_ev_node, '-dpng', sprintf('-r%d', dpi));
    fprintf('【节点级】EV高DPI预测曲线图已保存为: %s\n', output_filename_ev_node);
catch ME_save
    fprintf('保存图像失败: %s\n', ME_save.message);
end

fprintf('--- 绘图完成 ---\n');

% 恢复默认字体（可选，良好实践）
set(0, 'DefaultAxesFontName', 'default');
set(0, 'DefaultTextFontName', 'default');