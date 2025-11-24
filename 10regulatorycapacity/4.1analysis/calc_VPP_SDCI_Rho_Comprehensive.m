%% calc_VPP_SDCI_Rho_Comprehensive.m
% 功能：
% 1. 计算不同调节时长 (5min, 15min, 60min) 下 AC 和 EV 集群的互补性 (SDCI) 和相关性 (Rho)。
% 2. 计算不同激励电价下 AC 和 EV 集群的互补性 (SDCI) 和相关性 (Rho)。
% 3. 绘制共 8 张分析图表 (4张时长对比 + 4张电价趋势)。

clear; close all; clc;

% 设置默认字体和绘图参数
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultLineLineWidth', 1.5);

%% ========================================================================
%  第一部分：不同调节时长 (dt) 下的指标计算与绘图
% ========================================================================
fprintf('==========================================================\n');
fprintf('  开始执行任务 1：不同调节时长 (dt) 下的指标计算\n');
fprintf('==========================================================\n');

% 1. 定义文件名和参数
dt_labels = {'5min', '15min', '60min'};
ev_files_dt = {'main_potential_5min.mat', 'main_potential_15min.mat', 'main_potential_60min.mat'};
ac_files_dt = {'AC_Stateful_Simulation_Results_5min.mat', 'AC_Stateful_Simulation_Results_15min.mat', 'AC_Stateful_Simulation_Results_60min.mat'};

num_dt = length(dt_labels);

% 2. 初始化结果存储数组
res_dt_sdci_up   = zeros(1, num_dt);
res_dt_rho_up    = zeros(1, num_dt);
res_dt_sdci_down = zeros(1, num_dt);
res_dt_rho_down  = zeros(1, num_dt);

% 3. 循环处理每个 dt
for i = 1:num_dt
    fprintf('\n--- 正在处理 dt = %s ---\n', dt_labels{i});
    
    % 加载文件
    [ac_up, ac_down] = load_ac_data(ac_files_dt{i});
    [ev_up, ev_down] = load_ev_data(ev_files_dt{i});
    
    % 数据对齐与预处理
    [ac_up, ev_up] = align_data(ac_up, ev_up);
    [ac_down, ev_down] = align_data(abs(ac_down), abs(ev_down)); % 下调取绝对值
    
    % 计算指标
    n_dummy = ones(length(ac_up), 1); % 聚合模型视为 n=1
    
    % 上调
    res_dt_sdci_up(i) = local_calSDCI(n_dummy, n_dummy, ac_up, ev_up);
    res_dt_rho_up(i)  = local_calRho(n_dummy, ac_up, n_dummy, ev_up);
    
    % 下调
    res_dt_sdci_down(i) = local_calSDCI(n_dummy, n_dummy, ac_down, ev_down);
    res_dt_rho_down(i)  = local_calRho(n_dummy, ac_down, n_dummy, ev_down);
    
    fprintf('  SDCI(Up)=%.4f, Rho(Up)=%.4f\n', res_dt_sdci_up(i), res_dt_rho_up(i));
    fprintf('  SDCI(Down)=%.4f, Rho(Down)=%.4f\n', res_dt_sdci_down(i), res_dt_rho_down(i));
end

% 4. 绘制图表 (4张图)
plot_bar_chart(dt_labels, res_dt_sdci_up, '不同时长 SDCI (上调)', 'SDCI 值', 'SDCI_Up_Duration');
plot_bar_chart(dt_labels, res_dt_rho_up,  '不同时长 Rho (上调)',  'Rho 值',  'Rho_Up_Duration');
plot_bar_chart(dt_labels, res_dt_sdci_down, '不同时长 SDCI (下调)', 'SDCI 值', 'SDCI_Down_Duration');
plot_bar_chart(dt_labels, res_dt_rho_down,  '不同时长 Rho (下调)',  'Rho 值',  'Rho_Down_Duration');


%% ========================================================================
%  第二部分：不同激励电价下的指标计算与绘图
% ========================================================================
fprintf('\n==========================================================\n');
fprintf('  开始执行任务 2：不同激励电价下的指标计算\n');
fprintf('==========================================================\n');

% 1. 定义价格列表 (0 到 50，共10个点)
price_list = linspace(0, 50, 10);
num_prices = length(price_list);

% 2. 初始化结果存储数组
res_price_sdci_up   = nan(1, num_prices);
res_price_rho_up    = nan(1, num_prices);
res_price_sdci_down = nan(1, num_prices);
res_price_rho_down  = nan(1, num_prices);

% 3. 循环处理每个价格
for k = 1:num_prices
    current_p = price_list(k);
    
    % 构造文件名
    % AC 文件名格式: AC_Stateful_Simulation_Results_Price_5.6.mat
    ac_file_p = sprintf('AC_Stateful_Simulation_Results_Price_%.1f.mat', current_p);
    
    % EV 文件名格式: results_incentive_5.56.mat
    ev_file_p = sprintf('results_incentive_%.2f.mat', current_p);
    
    fprintf('\n--- 正在处理价格 P = %.2f ---\n', current_p);
    fprintf('  寻找 AC 文件: %s\n', ac_file_p);
    fprintf('  寻找 EV 文件: %s\n', ev_file_p);
    
    % 检查文件是否存在
    if ~exist(ac_file_p, 'file') || ~exist(ev_file_p, 'file')
        fprintf('  [警告] 文件缺失，跳过此价格点。\n');
        continue;
    end
    
    try
        % 加载文件
        [ac_up, ac_down] = load_ac_data(ac_file_p);
        [ev_up, ev_down] = load_ev_data(ev_file_p);
        
        % 数据对齐与预处理
        [ac_up, ev_up] = align_data(ac_up, ev_up);
        [ac_down, ev_down] = align_data(abs(ac_down), abs(ev_down)); % 下调取绝对值
        
        % 计算指标
        n_dummy = ones(length(ac_up), 1);
        
        % 上调
        res_price_sdci_up(k) = local_calSDCI(n_dummy, n_dummy, ac_up, ev_up);
        res_price_rho_up(k)  = local_calRho(n_dummy, ac_up, n_dummy, ev_up);
        
        % 下调
        res_price_sdci_down(k) = local_calSDCI(n_dummy, n_dummy, ac_down, ev_down);
        res_price_rho_down(k)  = local_calRho(n_dummy, ac_down, n_dummy, ev_down);
        
        fprintf('  SDCI(Up)=%.4f, Rho(Up)=%.4f\n', res_price_sdci_up(k), res_price_rho_up(k));
    catch ME
        fprintf('  [错误] 计算出错: %s\n', ME.message);
    end
end

% 4. 绘制图表 (4张折线图)
plot_line_chart(price_list, res_price_sdci_up, '激励电价 (元)', 'SDCI (上调) 随电价变化', 'SDCI 值', [0 0.4470 0.7410], 'Price_SDCI_Up');
plot_line_chart(price_list, res_price_rho_up,  '激励电价 (元)', 'Rho (上调) 随电价变化',  'Rho 值',  [0.8500 0.3250 0.0980], 'Price_Rho_Up');
plot_line_chart(price_list, res_price_sdci_down,'激励电价 (元)', 'SDCI (下调) 随电价变化', 'SDCI 值', [0.9290 0.6940 0.1250], 'Price_SDCI_Down');
plot_line_chart(price_list, res_price_rho_down, '激励电价 (元)', 'Rho (下调) 随电价变化',  'Rho 值',  [0.4940 0.1840 0.5560], 'Price_Rho_Down');

fprintf('\n全部计算与绘图完成。\n');


%% ========================================================================
%  辅助函数库
% ========================================================================

% --- 加载 AC 数据 ---
function [ac_up, ac_down] = load_ac_data(filename)
    if ~exist(filename, 'file'), error('AC文件不存在: %s', filename); end
    data = load(filename);
    if ~isfield(data, 'results'), error('%s 中无 results', filename); end
    
    % 兼容不同版本的字段名
    if isfield(data.results, 'Agg_P_Potential_Up_History')
        ac_up = data.results.Agg_P_Potential_Up_History;
    elseif isfield(data.results, 'AC_Up')
        ac_up = data.results.AC_Up;
    else
        error('无法找到 AC Up 数据');
    end
    
    if isfield(data.results, 'Agg_P_Potential_Down_History')
        ac_down = data.results.Agg_P_Potential_Down_History;
    elseif isfield(data.results, 'AC_Down')
        ac_down = data.results.AC_Down;
    else
        error('无法找到 AC Down 数据');
    end
    
    ac_up = ac_up(:); ac_down = ac_down(:);
end

% --- 加载 EV 数据 ---
function [ev_up, ev_down] = load_ev_data(filename)
    if ~exist(filename, 'file'), error('EV文件不存在: %s', filename); end
    data = load(filename);
    if ~isfield(data, 'results'), error('%s 中无 results', filename); end
    
    if isfield(data.results, 'EV_Up')
        ev_up = data.results.EV_Up;
    else
        error('无法找到 EV Up 数据');
    end
    
    if isfield(data.results, 'EV_Down')
        ev_down = data.results.EV_Down;
    else
        error('无法找到 EV Down 数据');
    end
    
    ev_up = ev_up(:); ev_down = ev_down(:);
end

% --- 数据对齐 ---
function [d1, d2] = align_data(d1, d2)
    len = min(length(d1), length(d2));
    d1 = d1(1:len);
    d2 = d2(1:len);
end

% --- 绘图：柱状图 (用于时长对比) ---
function plot_bar_chart(x_labels, y_data, title_str, ylabel_str, fig_name)
    figure('Name', fig_name, 'Color', 'w', 'Position', [200, 200, 500, 400]);
    b = bar(y_data);
    b.FaceColor = 'flat';
    b.CData(1,:) = [0.2 0.6 0.8]; % 5min 蓝
    b.CData(2,:) = [0.2 0.7 0.5]; % 15min 绿
    b.CData(3,:) = [0.8 0.4 0.2]; % 60min 橙
    
    set(gca, 'XTickLabel', x_labels);
    ylabel(ylabel_str);
    title(title_str);
    grid on;
    
    % 添加数值标签
    xtips = b.XEndPoints;
    ytips = b.YEndPoints;
    labels = string(round(b.YData, 4));
    text(xtips, ytips, labels, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 10);
    
    % 简单的边界美化
    if min(y_data) >= 0
        ylim([0, max(y_data)*1.2]);
    else
        ylim([min(y_data)*1.2, max(y_data)*1.2]);
    end
end

% --- 绘图：折线图 (用于电价趋势) ---
function plot_line_chart(x_data, y_data, xlabel_str, title_str, ylabel_str, line_color, fig_name)
    figure('Name', fig_name, 'Color', 'w', 'Position', [300, 300, 600, 400]);
    
    % 过滤掉 NaN 数据 (针对可能缺失的文件)
    valid_idx = ~isnan(y_data);
    
    plot(x_data(valid_idx), y_data(valid_idx), '-o', ...
        'LineWidth', 2, 'MarkerSize', 6, ...
        'Color', line_color, 'MarkerFaceColor', line_color);
    
    xlabel(xlabel_str);
    ylabel(ylabel_str);
    title(title_str);
    grid on;
    
    % 简单的边界美化
    xlim([min(x_data)-2, max(x_data)+2]);
end

% --- 指标计算函数 (SDCI) ---
function [SDCI] = local_calSDCI(n_AC, n_EV, deltaP_AC, deltaP_EV)
    AC = n_AC .* deltaP_AC;
    EV = n_EV .* deltaP_EV;
    
    min_vals = min(AC, EV);
    max_vals = max(AC, EV);
    
    total_max = sum(max_vals);
    if total_max == 0
        SDCI = 0;
    else
        SDCI = sum(min_vals) / total_max;
    end
end

% --- 指标计算函数 (Rho) ---
function [rho] = local_calRho(n_AC, deltaP_AC, n_EV, deltaP_EV)
    AC_series = n_AC .* deltaP_AC;
    EV_series = n_EV .* deltaP_EV;
    
    T = length(AC_series);
    if T < 2
        rho = NaN;
        return;
    end
    
    try
        rho = corr(AC_series, EV_series, 'Type', 'Spearman');
    catch
        % 手动计算
        rank_AC = local_tiedrank(AC_series);
        rank_EV = local_tiedrank(EV_series);
        d_sq = (rank_AC - rank_EV).^2;
        sum_d_sq = sum(d_sq);
        rho = 1 - (6 * sum_d_sq) / (T * (T^2 - 1));
    end
    rho = max(min(rho, 1), -1);
end

function r = local_tiedrank(x)
    [~, p] = sort(x);
    r = zeros(size(x));
    r(p) = 1:length(x);
end