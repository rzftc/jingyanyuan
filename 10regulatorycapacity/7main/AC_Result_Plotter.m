%% AC_Result_Plotter.m
% 功能：读取 AC_Stateful_Simulation_Results.mat 并复现所有分析图表
% 依赖：必须先运行修改后的 AC_main_Stateful_Sim_potential.m 以生成完整数据

clear; close all; clc;

% --- 1. 加载数据 ---
mat_filename = 'AC_Stateful_Simulation_Results.mat';
fprintf('正在加载数据文件: %s ...\n', mat_filename);

if ~exist(mat_filename, 'file')
    error('错误: 未找到文件 %s。\n请先运行仿真脚本 AC_main_Stateful_Sim_potential.m 生成数据。', mat_filename);
end

load(mat_filename);

if ~exist('results', 'var')
    error('错误: MAT文件中未找到 "results" 结构体。');
end

% --- 2. 数据提取与恢复 ---
fprintf('正在提取数据...\n');

% 检查并提取时间轴
if isfield(results, 'time_points')
    time_points = results.time_points;
elseif isfield(results, 'time_points_absolute')
    time_points = results.time_points_absolute;
else
    error('数据中缺失时间轴 (time_points)。');
end

% 提取 SOC 并恢复维度 (T x N)
if isfield(results, 'Individual_SOC_History_Transposed')
    Individual_SOC_History = results.Individual_SOC_History_Transposed';
elseif isfield(results, 'SOC_AC') % 兼容旧版命名
    Individual_SOC_History = results.SOC_AC';
else
    warning('数据中缺失 SOC 历史数据，图2, 3 将无法绘制。');
    Individual_SOC_History = [];
end

% 计算聚合 SOC
if ~isempty(Individual_SOC_History)
    Agg_SOC_History = mean(Individual_SOC_History, 2, 'omitnan');
    num_AC_participating = size(Individual_SOC_History, 2);
else
    Agg_SOC_History = [];
    num_AC_participating = 0;
end

% 提取温度数据并恢复维度 (T x N)
if isfield(results, 'Individual_Temp_History_Transposed')
    Individual_Temp_History = results.Individual_Temp_History_Transposed';
else
    warning('数据中缺失温度历史 (Individual_Temp_History)，图3 将无法绘制。');
    Individual_Temp_History = [];
end

% 提取其他关键变量
fields_to_load = {
    'Agg_P_Command_History', 'Agg_P_Achieved_History', ...
    'Individual_Power_History', 'Total_Power_History', ...
    'Agg_Baseline_Power', 'Agg_Total_Power', ...
    'Agg_P_Potential_Up_History', 'Agg_P_Potential_Down_History', ...
    'Agg_Model_Potential_Up_History', 'Agg_Model_Potential_Down_History'
};

% 动态加载变量到工作区，如果不存在则设为空
for i = 1:length(fields_to_load)
    fname = fields_to_load{i};
    if isfield(results, fname)
        eval([fname ' = results.' fname ';']);
    elseif isfield(results, 'AC_Up') && strcmp(fname, 'Agg_P_Potential_Up_History')
         Agg_P_Potential_Up_History = results.AC_Up; % 兼容旧名
    elseif isfield(results, 'AC_Down') && strcmp(fname, 'Agg_P_Potential_Down_History')
         Agg_P_Potential_Down_History = results.AC_Down; % 兼容旧名
    else
        warning('数据中缺失字段: %s，相关图表可能不完整。', fname);
        eval([fname ' = [];']);
    end
end

P_standby = 0.05; % 最小待机功率 (用于图5参考线)

%% 3. 开始绘图
fprintf('正在生成图表...\n');

% --- 图 1: 功率跟踪对比 ---
if ~isempty(Agg_P_Command_History) && ~isempty(Agg_P_Achieved_History)
    figure('Name', '功率跟踪对比 (理论分解)', 'Position', [100 100 1000 450]);
    ax1 = axes;
    hold(ax1, 'on');
    plot(ax1, time_points, Agg_P_Command_History, 'k:', 'LineWidth', 2.5, ...
        'DisplayName', '电网调节指令 (ΔP_s)');
    plot(ax1, time_points, Agg_P_Achieved_History, 'r-', 'LineWidth', 1.5, ...
        'DisplayName', '各空调的响应功率之和 (ΣΔP_j)');
    hold(ax1, 'off');
    xlabel(ax1, '时间 (小时)', 'FontSize', 12);
    ylabel(ax1, '聚合功率 (kW)', 'FontSize', 12);
    title(ax1, '图1：空调聚合响应功率 vs 电网指令', 'FontSize', 14);
    legend(ax1, 'show', 'Location', 'best');
    set(ax1, 'FontSize', 11);
    xticks(ax1, [0, 6, 12, 18, 24]);
    xticklabels(ax1, {'00:00', '06:00', '12:00', '18:00', '24:00'});
    xlim(ax1, [0, 24]);
    grid(ax1, 'on');
end

% --- 图 2: SOC状态对比 (多曲线) ---
if ~isempty(Individual_SOC_History) && ~isempty(Agg_SOC_History)
    figure('Name', 'SOC状态对比 (多曲线)', 'Position', [100 550 1000 450]);
    ax2 = axes;
    hold(ax2, 'on');
    h_individual = plot(ax2, time_points, Individual_SOC_History, 'LineWidth', 0.5);
    h_agg = plot(ax2, time_points, Agg_SOC_History, 'k--', 'LineWidth', 3, ...
        'DisplayName', '空调聚合模型的SOC (均值)');
    
    % 优化图例显示 (只显示一条单体曲线)
    if num_AC_participating > 0
        set(h_individual(1), 'DisplayName', '单体空调的SOC');
        if num_AC_participating > 1
            set(h_individual(2:end), 'HandleVisibility', 'off');
        end
        legend(ax2, [h_agg, h_individual(1)], 'Location', 'best', 'FontSize', 11);
    else
        legend(ax2, h_agg, 'Location', 'best', 'FontSize', 11);
    end
    
    hold(ax2, 'off');
    xlabel(ax2, '时间 (小时)', 'FontSize', 12);
    ylabel(ax2, '空调SOC', 'FontSize', 12);
    title(ax2, '图2：单体空调SOC 与 聚合模型SOC 对比 (多曲线)', 'FontSize', 14);
    set(ax2, 'FontSize', 11);
    xticks(ax2, [0, 6, 12, 18, 24]);
    xticklabels(ax2, {'00:00', '06:00', '12:00', '18:00', '24:00'});
    xlim(ax2, [0, 24]);
    ylim(ax2, [-0.1, 1.1]);
    grid(ax2, 'on');
end

% --- 图 3: 室内温度变化 (反推) ---
if ~isempty(Individual_Temp_History)
    figure('Name', '室内温度变化 (反推)', 'Position', [100 300 1000 450]);
    ax3 = axes;
    hold(ax3, 'on');
    plot(ax3, time_points, Individual_Temp_History, 'LineWidth', 0.5);
    hold(ax3, 'off');
    xlabel(ax3, '时间 (小时)', 'FontSize', 12);
    ylabel(ax3, '室内温度 (°C)', 'FontSize', 12);
    title(ax3, '图3：单体空调室内温度变化 (基于SOC反推)', 'FontSize', 14);
    xlim(ax3, [0, 24]); grid(ax3, 'on');
    set(ax3, 'FontSize', 11);
    xticks(ax3, [0, 6, 12, 18, 24]);
    xticklabels(ax3, {'00:00', '06:00', '12:00', '18:00', '24:00'});
end

% --- 图 4: 单体空调(调节)功率 ---
if ~isempty(Individual_Power_History)
    figure('Name', '单体空调调节功率曲线', 'Position', [100 400 1000 450]);
    ax4 = axes;
    hold(ax4, 'on');
    plot(ax4, time_points, Individual_Power_History, 'LineWidth', 0.5); 
    yline(ax4, 0, 'k--', 'LineWidth', 2, 'DisplayName', '0 kW 参考线');
    hold(ax4, 'off');
    xlabel(ax4, '时间 (小时)', 'FontSize', 12);
    ylabel(ax4, '单体调节功率 (kW)', 'FontSize', 12);
    title(ax4, '图4：单体空调实际调节功率 (ΔP_j) 变化', 'FontSize', 14);
    xlim(ax4, [0, 24]); grid(ax4, 'on');
    set(ax4, 'FontSize', 11);
    xticks(ax4, [0, 6, 12, 18, 24]);
    xticklabels(ax4, {'00:00', '06:00', '12:00', '18:00', '24:00'});
end

% --- 图 5: 单体空调总制冷功率 ---
if ~isempty(Total_Power_History)
    figure('Name', '单体空调总制冷功率', 'Position', [100 500 1000 450]);
    ax5 = axes;
    hold(ax5, 'on');
    plot(ax5, time_points, Total_Power_History, 'LineWidth', 0.5);
    yline(ax5, P_standby, 'k--', 'LineWidth', 2, 'DisplayName', '最小待机功率');
    hold(ax5, 'off');
    xlabel(ax5, '时间 (小时)', 'FontSize', 12);
    ylabel(ax5, '单体总制冷功率 (kW)', 'FontSize', 12);
    title(ax5, '图5：单体空调总制冷功率 (P_{base} + ΔP_j) 变化', 'FontSize', 14);
    xlim(ax5, [0, 24]); grid(ax5, 'on');
    set(ax5, 'FontSize', 11);
    xticks(ax5, [0, 6, 12, 18, 24]);
    xticklabels(ax5, {'00:00', '06:00', '12:00', '18:00', '24:00'});
end

% --- 图 6: 聚合功率对比 ---
if ~isempty(Agg_Baseline_Power) && ~isempty(Agg_Total_Power)
    figure('Name', '聚合功率对比', 'Position', [100 600 1000 450]);
    ax6 = axes;
    hold(ax6, 'on');
    plot(ax6, time_points, Agg_Baseline_Power, 'b--', 'LineWidth', 2, ...
        'DisplayName', '所有AC的基线功率之和 (ΣP_{base})');
    plot(ax6, time_points, Agg_Total_Power, 'r-', 'LineWidth', 2, ...
        'DisplayName', '所有AC的总制冷功率之和 (ΣP_{total})');
    hold(ax6, 'off');
    xlabel(ax6, '时间 (小时)', 'FontSize', 12);
    ylabel(ax6, '聚合功率 (kW)', 'FontSize', 12);
    title(ax6, '图6：聚合基线功率 vs 聚合总制冷功率', 'FontSize', 14);
    legend(ax6, 'show', 'Location', 'best', 'FontSize', 11);
    set(ax6, 'FontSize', 11);
    xticks(ax6, [0, 6, 12, 18, 24]);
    xticklabels(ax6, {'00:00', '06:00', '12:00', '18:00', '24:00'});
    xlim(ax6, [0, 24]); grid(ax6, 'on');
end

% --- 图 7: 聚合潜力对比 ---
if ~isempty(Agg_P_Potential_Up_History) && ~isempty(Agg_Model_Potential_Up_History)
    figure('Name', '聚合潜力对比', 'Position', [100 700 1000 450]);
    ax7 = axes;
    hold(ax7, 'on');
    plot(ax7, time_points, Agg_P_Potential_Up_History, 'b-', 'LineWidth', 1.5, ...
        'DisplayName', '单体累加上调 (ΣP_{+,j})');
    plot(ax7, time_points, Agg_P_Potential_Down_History, 'r-', 'LineWidth', 1.5, ...
        'DisplayName', '单体累加下调 (ΣP_{-,j})');
    plot(ax7, time_points, Agg_Model_Potential_Up_History, 'b--', 'LineWidth', 2, ...
        'DisplayName', '聚合模型上调 (P_{+,agg})');
    plot(ax7, time_points, Agg_Model_Potential_Down_History, 'r--', 'LineWidth', 2, ...
        'DisplayName', '聚合模型下调 (P_{-,agg})');
    hold(ax7, 'off');
    xlabel(ax7, '时间 (小时)', 'FontSize', 12);
    ylabel(ax7, '调节潜力 (kW)', 'FontSize', 12);
    title(ax7, '图7：调节潜力对比 (单体累加 vs 聚合模型计算)', 'FontSize', 14);
    legend(ax7, 'show', 'Location', 'best', 'FontSize', 11);
    set(ax7, 'FontSize', 11);
    xticks(ax7, [0, 6, 12, 18, 24]);
    xticklabels(ax7, {'00:00', '06:00', '12:00', '18:00', '24:00'});
    xlim(ax7, [0, 24]); grid(ax7, 'on');
end

fprintf('绘图程序执行完毕。\n');