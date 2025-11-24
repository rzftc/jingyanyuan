%% calc_SDCI_Rho_5min_Specific.m
% 功能：基于指定的 5min 仿真结果文件计算 AC 和 EV 的 SDCI 和 Rho 指标
% 输入文件：
%   1. AC_Stateful_Simulation_Results_5min.mat (AC 数据)
%   2. main_potential_5min.mat (EV 数据)

clear; close all; clc;

%% 1. 文件路径定义
ac_file = 'AC_Stateful_Simulation_Results_5min.mat';
ev_file = 'main_potential_5min.mat';

%% 2. 加载 AC 数据
fprintf('正在加载 AC 数据: %s ...\n', ac_file);
if ~exist(ac_file, 'file')
    error('未找到文件 %s，请确保该文件在当前路径下。', ac_file);
end
ac_data_struct = load(ac_file);
if ~isfield(ac_data_struct, 'results')
    error('AC 文件中未找到 "results" 结构体。');
end

% 提取 AC 潜力序列
% 注意：AC 代码中 AC_Up 通常为正，AC_Down 为负
if isfield(ac_data_struct.results, 'AC_Up')
    AC_Up_Raw = ac_data_struct.results.AC_Up;
elseif isfield(ac_data_struct.results, 'Agg_P_Potential_Up_History')
     % 兼容可能的字段名变化
    AC_Up_Raw = ac_data_struct.results.Agg_P_Potential_Up_History;
else
    error('无法在 AC 结果中找到上调潜力数据 (AC_Up 或 Agg_P_Potential_Up_History)');
end

if isfield(ac_data_struct.results, 'AC_Down')
    AC_Down_Raw = ac_data_struct.results.AC_Down;
elseif isfield(ac_data_struct.results, 'Agg_P_Potential_Down_History')
    AC_Down_Raw = ac_data_struct.results.Agg_P_Potential_Down_History;
else
    error('无法在 AC 结果中找到下调潜力数据 (AC_Down 或 Agg_P_Potential_Down_History)');
end

%% 3. 加载 EV 数据
fprintf('正在加载 EV 数据: %s ...\n', ev_file);
if ~exist(ev_file, 'file')
    error('未找到文件 %s，请确保该文件在当前路径下。', ev_file);
end
ev_data_struct = load(ev_file);
if ~isfield(ev_data_struct, 'results')
    error('EV 文件中未找到 "results" 结构体。');
end

% 提取 EV 潜力序列
% EV 代码中 EV_Up 为正，EV_Down 为负
if isfield(ev_data_struct.results, 'EV_Up')
    EV_Up_Raw = ev_data_struct.results.EV_Up;
else
    error('无法在 EV 结果中找到 EV_Up 数据');
end

if isfield(ev_data_struct.results, 'EV_Down')
    EV_Down_Raw = ev_data_struct.results.EV_Down;
else
    error('无法在 EV 结果中找到 EV_Down 数据');
end

%% 4. 数据预处理与对齐
% 确保数据是列向量
AC_Up_Raw = AC_Up_Raw(:); AC_Down_Raw = AC_Down_Raw(:);
EV_Up_Raw = EV_Up_Raw(:); EV_Down_Raw = EV_Down_Raw(:);

% 获取最小长度以对齐数据
len = min([length(AC_Up_Raw), length(EV_Up_Raw)]);
fprintf('数据对齐: 将数据截取至 %d 个时间步。\n', len);

% 截取数据
AC_Up = AC_Up_Raw(1:len);
EV_Up = EV_Up_Raw(1:len);

% 对于下调潜力，取绝对值进行指标计算 (表示调节能力的大小)
AC_Down_Abs = abs(AC_Down_Raw(1:len));
EV_Down_Abs = abs(EV_Down_Raw(1:len));

%% 5. 指标计算
% 构造辅助输入：
% 由于我们要计算的是“集群”层面的指标，我们将集群视为一个整体。
% 相当于 n=1, deltaP = 聚合功率。
n_dummy = ones(len, 1); 

% --- 计算上调 (Up) 指标 ---
SDCI_Up = local_calSDCI(n_dummy, n_dummy, AC_Up, EV_Up);
Rho_Up  = local_calRho(n_dummy, AC_Up, n_dummy, EV_Up);

% --- 计算下调 (Down) 指标 ---
SDCI_Down = local_calSDCI(n_dummy, n_dummy, AC_Down_Abs, EV_Down_Abs);
Rho_Down  = local_calRho(n_dummy, AC_Down_Abs, n_dummy, EV_Down_Abs);

%% 6. 结果显示与绘图
fprintf('\n==============================================\n');
fprintf('        AC 与 EV 集群互补性分析结果        \n');
fprintf('==============================================\n');
fprintf('【上调能力 (Up-Regulation)】\n');
fprintf('  SDCI (互补性指数):   %.4f\n', SDCI_Up);
fprintf('  Rho  (相关系数):     %.4f\n', Rho_Up);
fprintf('\n');
fprintf('【下调能力 (Down-Regulation)】\n');
fprintf('  SDCI (互补性指数):   %.4f\n', SDCI_Down);
fprintf('  Rho  (相关系数):     %.4f\n', Rho_Down);
fprintf('==============================================\n');

% --- 绘图 ---
figure('Name', 'AC-EV 5min 互补性指标', 'Position', [300, 300, 600, 400], 'Color', 'w');
bar_data = [SDCI_Up, SDCI_Down; Rho_Up, Rho_Down];
b = bar(bar_data);
b(1).FaceColor = [0.2 0.6 0.8]; % 蓝色 (Up)
b(2).FaceColor = [0.8 0.4 0.2]; % 橙色 (Down)

set(gca, 'XTickLabel', {'SDCI (互补性)', 'Rho (相关性)'}, 'FontSize', 12);
ylabel('指标值', 'FontSize', 12);
legend({'上调 (Up)', '下调 (Down)'}, 'Location', 'best');
grid on;
ylim([-1.1, 1.1]); % Rho 范围是 -1 到 1，SDCI 是 0 到 1

% 添加数值标签
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(round(b(1).YData, 3));
text(xtips1, ytips1, labels1, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(round(b(2).YData, 3));
text(xtips2, ytips2, labels2, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');


%% ---------------------------------------------------------
%  本地辅助函数 (确保脚本独立运行，不依赖外部文件)
% ---------------------------------------------------------

function [SDCI] = local_calSDCI(n_AC, n_EV, deltaP_AC, deltaP_EV)
% 计算同方向互补性指数 (SDCI)
    % 计算总调节能力
    AC = n_AC .* deltaP_AC;
    EV = n_EV .* deltaP_EV;
    
    % 计算逐时段min/max
    min_vals = min(AC, EV);
    max_vals = max(AC, EV);
    
    % 处理分母为零的情况
    total_max = sum(max_vals);
    if total_max == 0
        SDCI = 0;
    else
        SDCI = sum(min_vals) / total_max;
    end
end

function [rho] = local_calRho(n_AC, deltaP_AC, n_EV, deltaP_EV)
% 计算斯皮尔曼秩相关系数
    AC_series = n_AC .* deltaP_AC;
    EV_series = n_EV .* deltaP_EV;
    
    T = length(AC_series);
    if T < 2
        rho = NaN;
        return;
    end
    
    % 使用 MATLAB 内置 corr 函数计算 Spearman 相关系数
    % 如果没有统计工具箱，使用简化公式（不处理并列秩）或原有逻辑
    try
        rho = corr(AC_series, EV_series, 'Type', 'Spearman');
    catch
        % 手动计算 (calRho.m 的逻辑)
        rank_AC = tiedrank_local(AC_series);
        rank_EV = tiedrank_local(EV_series);
        d_sq = (rank_AC - rank_EV).^2;
        sum_d_sq = sum(d_sq);
        rho = 1 - (6 * sum_d_sq) / (T * (T^2 - 1));
    end
    
    % 数值边界处理
    rho = max(min(rho, 1), -1); 
end

function r = tiedrank_local(x)
    % 简化的 tiedrank 实现，防止缺少统计工具箱
    [~, p] = sort(x);
    r = zeros(size(x));
    r(p) = 1:length(x);
    % 处理并列值(简化版：取平均秩)
    % 为保证代码简洁，此处仅做基础排序，如需精确处理并列可使用 MATLAB 内置 tiedrank
end