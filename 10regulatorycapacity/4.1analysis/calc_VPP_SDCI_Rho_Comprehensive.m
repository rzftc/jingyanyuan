%% calc_VPP_SDCI_Rho_Comprehensive.m
% 功能：
% 1. 计算不同调节时长 (5min, 15min, 60min) 下 AC 和 EV 集群的互补性 (SDCI) 和相关性 (Rho)。
% 2. 计算不同激励电价下 AC 和 EV 集群的互补性 (SDCI) 和相关性 (Rho)。
% 3. 绘制共 8 张分析图表 (4张时长对比 + 4张电价趋势)。
%    [修改] 保存为高DPI PNG，无标题，中文文件名和图例。
%    [修改] 统一时间轴：将 AC 数据 (0:24) 转化为 (6:30)，与 EV 数据对齐。

clear; close all; clc;

% 设置默认字体和绘图参数 (推荐使用支持中文的字体)
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesFontName', 'Microsoft YaHei'); 
set(0, 'DefaultTextFontName', 'Microsoft YaHei');

%% ========================================================================
%  第一部分：不同调节时长 (dt) 下的指标计算与绘图
% ========================================================================
fprintf('==========================================================\n');
fprintf('  开始执行任务 1：不同调节时长 (dt) 下的指标计算\n');
fprintf('==========================================================\n');

% 1. 定义文件名和参数 (标签改为中文)
dt_labels = {'5分钟', '15分钟', '60分钟'};
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
    try
        [ac_up, ac_down] = load_ac_data(ac_files_dt{i});
        [ev_up, ev_down] = load_ev_data(ev_files_dt{i});
        
        % ==================== [修改] 时间轴对齐逻辑 ====================
        % AC 原数据: 00:00 -> 24:00 (24h)
        % EV 原数据: 06:00 -> 30:00 (24h)
        % 操作: 将 AC 数据移位，使其变为 06:00 -> 30:00
        % 方法: 新 AC = [原AC(06:00-24:00); 原AC(00:00-06:00)]
        
        len_ac = length(ac_up);
        % 计算对应 6 小时的数据点数 (假设 len_ac 对应 24 小时)
        idx_6h = round(len_ac * (6/24));
        
        if idx_6h > 0 && idx_6h < len_ac
            % 拼接数据
            ac_up_shifted = [ac_up(idx_6h+1:end); ac_up(1:idx_6h)];
            ac_down_shifted = [ac_down(idx_6h+1:end); ac_down(1:idx_6h)];
            
            ac_up = ac_up_shifted;
            ac_down = ac_down_shifted;
            fprintf('  已执行时间轴对齐：AC 数据从 0:24 调整为 6:30 (循环移位 %d 点)\n', idx_6h);
        else
            warning('  AC 数据长度异常，跳过时间轴对齐。');
        end
        % ===============================================================
        
        % 数据长度截取 (取两者最小值，防止微小误差)
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
    catch ME
        fprintf('  [错误] 处理 dt=%s 时出错: %s\n', dt_labels{i}, ME.message);
    end
end

% 4. 绘制图表 (4张图) - 传入中文文件名，移除标题参数
plot_bar_chart(dt_labels, res_dt_sdci_up, 'SDCI 值', '不同时长_SDCI_上调');
plot_bar_chart(dt_labels, res_dt_rho_up,  'Rho 值',  '不同时长_Rho_上调');
plot_bar_chart(dt_labels, res_dt_sdci_down, 'SDCI 值', '不同时长_SDCI_下调');
plot_bar_chart(dt_labels, res_dt_rho_down,  'Rho 值',  '不同时长_Rho_下调');


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
    ac_file_p = sprintf('AC_Stateful_Simulation_Results_Price_%.1f.mat', current_p);
    
    % EV 文件名格式 (根据您的实际文件名调整)
    ev_file_p = sprintf('results_incentive_%.2f.mat', current_p); 
    
    fprintf('\n--- 正在处理价格 P = %.2f ---\n', current_p);
    
    % 检查文件是否存在
    if ~exist(ac_file_p, 'file') 
        fprintf('  [警告] AC文件缺失: %s，跳过此价格点。\n', ac_file_p);
        continue;
    end
    if ~exist(ev_file_p, 'file')
        % 尝试另一种常见的命名格式
        ev_file_p_alt = sprintf('results_chunk_ev_price_%.1f.mat', current_p);
        if exist(ev_file_p_alt, 'file')
            ev_file_p = ev_file_p_alt;
        else
             fprintf('  [警告] EV文件缺失: %s，跳过此价格点。\n', ev_file_p);
             continue;
        end
    end
    
    try
        % 加载文件
        [ac_up, ac_down] = load_ac_data(ac_file_p);
        [ev_up, ev_down] = load_ev_data(ev_file_p);
        
        % ==================== [修改] 时间轴对齐逻辑 ====================
        % 同样应用到价格循环中
        len_ac = length(ac_up);
        idx_6h = round(len_ac * (6/24)); 
        
        if idx_6h > 0 && idx_6h < len_ac
            ac_up_shifted = [ac_up(idx_6h+1:end); ac_up(1:idx_6h)];
            ac_down_shifted = [ac_down(idx_6h+1:end); ac_down(1:idx_6h)];
            
            ac_up = ac_up_shifted;
            ac_down = ac_down_shifted;
            % fprintf('  已对齐 AC 时间轴 (6:30)。\n');
        end
        % ===============================================================
        
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

% 4. 绘制图表 (4张折线图) - 传入中文文件名，移除标题参数
plot_line_chart(price_list, res_price_sdci_up, '激励电价 (元)', 'SDCI 值', [0 0.4470 0.7410], '电价敏感度_SDCI_上调');
plot_line_chart(price_list, res_price_rho_up,  '激励电价 (元)', 'Rho 值',  [0.8500 0.3250 0.0980], '电价敏感度_Rho_上调');
plot_line_chart(price_list, res_price_sdci_down,'激励电价 (元)', 'SDCI 值', [0.9290 0.6940 0.1250], '电价敏感度_SDCI_下调');
plot_line_chart(price_list, res_price_rho_down, '激励电价 (元)', 'Rho 值',  [0.4940 0.1840 0.5560], '电价敏感度_Rho_下调');

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
function plot_bar_chart(x_labels, y_data, ylabel_str, file_name_cn)
    figure('Name', file_name_cn, 'Color', 'w', 'Position', [200, 200, 600, 450]);
    b = bar(y_data);
    b.FaceColor = 'flat';
    % 注意：这里假设有3个柱子
    if length(y_data) >= 3
        b.CData(1,:) = [0.2 0.6 0.8]; % 5min 蓝
        b.CData(2,:) = [0.2 0.7 0.5]; % 15min 绿
        b.CData(3,:) = [0.8 0.4 0.2]; % 60min 橙
    else
        b.FaceColor = [0.2 0.6 0.8];
    end
    
    set(gca, 'XTickLabel', x_labels, 'FontSize', 12, 'FontName', 'Microsoft YaHei');
    ylabel(ylabel_str, 'FontName', 'Microsoft YaHei', 'FontSize', 12);
    
    grid on;
    
    % --- 动态添加数值标签 ---
    xtips = b.XEndPoints;
    ytips = b.YEndPoints;
    labels = string(round(b.YData, 4));
    
    for k = 1:length(ytips)
        if ytips(k) >= 0
            va = 'bottom'; v_offset = 0.01 * max(abs(y_data));
        else
            va = 'top'; v_offset = -0.01 * max(abs(y_data));
        end
        text(xtips(k), ytips(k) + v_offset, labels(k), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', va, ...
            'FontSize', 10, 'FontName', 'Arial'); 
    end

    % --- 边界美化 ---
    y_max = max(y_data); y_min = min(y_data);
    padding = 0.15; abs_max = max(abs(y_data)); if abs_max == 0, abs_max = 1; end
    margin = abs_max * padding;
    
    if y_min >= 0
        ylim([0, y_max + margin]);
    elseif y_max <= 0
        ylim([y_min - margin, 0]);
    else
        ylim([y_min - margin, y_max + margin]);
    end
    
    yline(0, 'k-', 'LineWidth', 1.0);
    legend({'计算数值'}, 'Location', 'best', 'FontName', 'Microsoft YaHei');

    print(gcf, [file_name_cn '.png'], '-dpng', '-r300');
    fprintf('  已保存图片: %s.png\n', file_name_cn);
end

% --- 绘图：折线图 (用于电价趋势) ---
function plot_line_chart(x_data, y_data, xlabel_str, ylabel_str, line_color, file_name_cn)
    figure('Name', file_name_cn, 'Color', 'w', 'Position', [300, 300, 600, 450]);
    valid_idx = ~isnan(y_data);
    
    plot(x_data(valid_idx), y_data(valid_idx), '-o', ...
        'LineWidth', 2, 'MarkerSize', 6, ...
        'Color', line_color, 'MarkerFaceColor', line_color, ...
        'DisplayName', '指标趋势'); 
    
    xlabel(xlabel_str, 'FontName', 'Microsoft YaHei', 'FontSize', 12);
    ylabel(ylabel_str, 'FontName', 'Microsoft YaHei', 'FontSize', 12);
    
    grid on;
    set(gca, 'FontSize', 12, 'FontName', 'Microsoft YaHei');
    legend('show', 'Location', 'best', 'FontName', 'Microsoft YaHei');
    xlim([min(x_data)-2, max(x_data)+2]);
    
    print(gcf, [file_name_cn '.png'], '-dpng', '-r300');
    fprintf('  已保存图片: %s.png\n', file_name_cn);
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