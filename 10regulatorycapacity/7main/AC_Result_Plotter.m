%% AC_Result_Plotter.m
% 功能：读取仿真结果并绘制高DPI、无标题、中文图例的分析图表
%       1. 读取 AC_Stateful_Simulation_Results.mat (单次仿真)
%       2. 读取 results_AC 文件夹下的批量结果 (不同电价)
% 依赖：AC_main_Stateful_Sim_potential_diff_inc.m 生成的数据

clear; close all; clc;

% 设置默认字体为支持中文的字体 (防止乱码)
set(0, 'DefaultAxesFontName', 'Microsoft YaHei'); 
set(0, 'DefaultTextFontName', 'Microsoft YaHei');

%% === 第一部分：单次仿真结果绘图 ===

mat_filename = 'AC_Stateful_Simulation_Results.mat';
fprintf('检查单次仿真数据文件: %s ...\n', mat_filename);

has_single_result = false;
if exist(mat_filename, 'file')
    load(mat_filename);
    if exist('results', 'var')
        has_single_result = true;
    else
        warning('文件存在但未找到 "results" 结构体。');
    end
else
    fprintf('提示: 未找到单次仿真文件。跳过单体详细分析，直接进行多场景对比。\n');
end

if has_single_result
    % --- 1. 数据提取 ---
    fprintf('正在提取单次仿真数据...\n');
    
    if isfield(results, 'time_points')
        time_points = results.time_points;
    elseif isfield(results, 'time_points_absolute')
        time_points = results.time_points_absolute;
    else
        error('数据中缺失时间轴。');
    end
    
    % 提取数据 (兼容不同版本的字段名)
    if isfield(results, 'Individual_SOC_History_Transposed')
        Individual_SOC_History = results.Individual_SOC_History_Transposed';
    elseif isfield(results, 'SOC_AC')
        Individual_SOC_History = results.SOC_AC';
    else
        Individual_SOC_History = [];
    end
    
    if ~isempty(Individual_SOC_History)
        Agg_SOC_History = mean(Individual_SOC_History, 2, 'omitnan');
        num_AC_participating = size(Individual_SOC_History, 2);
    else
        Agg_SOC_History = [];
        num_AC_participating = 0;
    end

    if isfield(results, 'Individual_Temp_History_Transposed')
        Individual_Temp_History = results.Individual_Temp_History_Transposed';
    else
        Individual_Temp_History = [];
    end
    
    % 动态提取其他变量
    fields_map = {
        'Agg_P_Command_History', 'Agg_P_Command_History';
        'Agg_P_Achieved_History', 'Agg_P_Achieved_History';
        'Individual_Power_History', 'Individual_Power_History';
        'Total_Power_History', 'Total_Power_History';
        'Agg_Baseline_Power', 'Agg_Baseline_Power';
        'Agg_Total_Power', 'Agg_Total_Power';
        'Agg_P_Potential_Up_History', 'AC_Up';
        'Agg_P_Potential_Down_History', 'AC_Down';
        'Agg_Model_Potential_Up_History', 'Agg_Model_Potential_Up_History';
        'Agg_Model_Potential_Down_History', 'Agg_Model_Potential_Down_History'
    };

    for i = 1:size(fields_map, 1)
        target_var = fields_map{i, 1};
        source_field = fields_map{i, 2};
        if isfield(results, source_field)
            eval([target_var ' = results.' source_field ';']);
        elseif isfield(results, target_var) % 尝试直接匹配目标名
            eval([target_var ' = results.' target_var ';']);
        else
            eval([target_var ' = [];']);
        end
    end
    
    P_standby = 0.05; 

    % --- 2. 开始绘图 (图1-7) ---
    fprintf('正在生成单次仿真图表...\n');
    
    % 图 1: 功率跟踪对比
    if ~isempty(Agg_P_Command_History) && ~isempty(Agg_P_Achieved_History)
        figure('Position', [100 100 1000 450]);
        plot(time_points, Agg_P_Command_History, 'k:', 'LineWidth', 2.5, 'DisplayName', '电网调节指令');
        hold on;
        plot(time_points, Agg_P_Achieved_History, 'r-', 'LineWidth', 1.5, 'DisplayName', '聚合响应功率');
        hold off;
        xlabel('时间 (小时)', 'FontSize', 12);
        ylabel('功率 (kW)', 'FontSize', 12);
        legend('show', 'Location', 'best', 'FontSize', 11);
        set(gca, 'FontSize', 11);
        xlim([0, 24]); grid on;
        % 无标题
        print(gcf, '图1_功率跟踪对比.png', '-dpng', '-r300');
    end
    
    % 图 2: SOC状态对比
    if ~isempty(Individual_SOC_History) && ~isempty(Agg_SOC_History)
        figure('Position', [100 550 1000 450]);
        hold on;
        h_ind = plot(time_points, Individual_SOC_History, 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
        h_agg = plot(time_points, Agg_SOC_History, 'k--', 'LineWidth', 2.5, 'DisplayName', '聚合SOC (均值)');
        % 只给第一条单体曲线加图例，避免图例过长
        if num_AC_participating > 0
            set(h_ind(1), 'DisplayName', '单体空调SOC');
            if num_AC_participating > 1, set(h_ind(2:end), 'HandleVisibility', 'off'); end
            legend([h_agg, h_ind(1)], 'Location', 'best', 'FontSize', 11);
        else
            legend(h_agg, 'Location', 'best', 'FontSize', 11);
        end
        hold off;
        xlabel('时间 (小时)', 'FontSize', 12);
        ylabel('SOC', 'FontSize', 12);
        set(gca, 'FontSize', 11);
        xlim([0, 24]); ylim([-0.1, 1.1]); grid on;
        print(gcf, '图2_SOC状态对比.png', '-dpng', '-r300');
    end
    
    % 图 3: 室内温度变化
    if ~isempty(Individual_Temp_History)
        figure('Position', [100 300 1000 450]);
        plot(time_points, Individual_Temp_History, 'LineWidth', 0.5);
        xlabel('时间 (小时)', 'FontSize', 12);
        ylabel('温度 (°C)', 'FontSize', 12);
        xlim([0, 24]); grid on;
        set(gca, 'FontSize', 11);
        print(gcf, '图3_室内温度变化.png', '-dpng', '-r300');
    end
    
    % 图 4: 单体调节功率
    if ~isempty(Individual_Power_History)
        figure('Position', [100 400 1000 450]);
        plot(time_points, Individual_Power_History, 'LineWidth', 0.5);
        yline(0, 'k--', 'LineWidth', 1.5);
        xlabel('时间 (小时)', 'FontSize', 12);
        ylabel('功率 (kW)', 'FontSize', 12);
        xlim([0, 24]); grid on;
        set(gca, 'FontSize', 11);
        print(gcf, '图4_单体调节功率.png', '-dpng', '-r300');
    end
    
    % 图 5: 单体总制冷功率
    if ~isempty(Total_Power_History)
        figure('Position', [100 500 1000 450]);
        plot(time_points, Total_Power_History, 'LineWidth', 0.5);
        yline(P_standby, 'k--', 'LineWidth', 1.5, 'DisplayName', '待机功率');
        xlabel('时间 (小时)', 'FontSize', 12);
        ylabel('功率 (kW)', 'FontSize', 12);
        xlim([0, 24]); grid on;
        set(gca, 'FontSize', 11);
        print(gcf, '图5_单体总制冷功率.png', '-dpng', '-r300');
    end
    
    % 图 6: 聚合功率对比
    if ~isempty(Agg_Baseline_Power) && ~isempty(Agg_Total_Power)
        figure('Position', [100 600 1000 450]);
        hold on;
        plot(time_points, Agg_Baseline_Power, 'b--', 'LineWidth', 2, 'DisplayName', '聚合基线功率');
        plot(time_points, Agg_Total_Power, 'r-', 'LineWidth', 2, 'DisplayName', '聚合总制冷功率');
        hold off;
        xlabel('时间 (小时)', 'FontSize', 12);
        ylabel('功率 (kW)', 'FontSize', 12);
        legend('show', 'Location', 'best', 'FontSize', 11);
        set(gca, 'FontSize', 11);
        xlim([0, 24]); grid on;
        print(gcf, '图6_聚合功率对比.png', '-dpng', '-r300');
    end
    
    % 图 7: 聚合潜力对比
    if ~isempty(Agg_P_Potential_Up_History)
        figure('Position', [100 700 1000 450]);
        hold on;
        plot(time_points, Agg_P_Potential_Up_History, 'b-', 'LineWidth', 1.5, 'DisplayName', '单体累加上调潜力');
        plot(time_points, Agg_P_Potential_Down_History, 'r-', 'LineWidth', 1.5, 'DisplayName', '单体累加下调潜力');
        if ~isempty(Agg_Model_Potential_Up_History)
            plot(time_points, Agg_Model_Potential_Up_History, 'b--', 'LineWidth', 2, 'DisplayName', '聚合模型上调潜力');
            plot(time_points, Agg_Model_Potential_Down_History, 'r--', 'LineWidth', 2, 'DisplayName', '聚合模型下调潜力');
        end
        hold off;
        xlabel('时间 (小时)', 'FontSize', 12);
        ylabel('功率 (kW)', 'FontSize', 12);
        legend('show', 'Location', 'best', 'FontSize', 11);
        set(gca, 'FontSize', 11);
        xlim([0, 24]); grid on;
        print(gcf, '图7_聚合潜力对比.png', '-dpng', '-r300');
    end
    
    fprintf('单次仿真绘图完成。\n');
end

%% === 第二部分：批量结果绘图 (不同电价下的聚合制冷功率) ===

results_dir = 'data/results_AC';
fprintf('\n------------------------------------------------------\n');
fprintf('检查批量仿真文件夹 "%s" ...\n', results_dir);

if exist(results_dir, 'dir')
    
    % 搜索所有文件
    file_pattern = fullfile(results_dir, 'AC_Stateful_Simulation_Results_Price_*.mat');
    mat_files = dir(file_pattern);
    
    if isempty(mat_files)
        fprintf('  警告: 在 "%s" 中未找到结果文件。\n', results_dir);
    else
        fprintf('  发现 %d 个结果文件，开始读取...\n', length(mat_files));
        
        % 数据容器
        data_list = struct('price', {}, 'total_power', {}, 'time_points', {});
        
        for i = 1:length(mat_files)
            full_path = fullfile(mat_files(i).folder, mat_files(i).name);
            try
                temp_data = load(full_path);
                if isfield(temp_data, 'results')
                    res = temp_data.results;
                    
                    % 获取价格 (优先从结构体)
                    p = NaN;
                    if isfield(res, 'current_price')
                        p = res.current_price;
                    else
                        % 尝试从文件名解析
                        tokens = regexp(mat_files(i).name, 'Price_([\d\.]+).mat', 'tokens');
                        if ~isempty(tokens), p = str2double(tokens{1}{1}); end
                    end
                    
                    % 获取聚合总功率
                    if ~isnan(p) && isfield(res, 'Agg_Total_Power') && isfield(res, 'time_points')
                        data_list(end+1).price = p;
                        data_list(end).total_power = res.Agg_Total_Power;
                        
                        % 兼容时间轴名称
                        if isfield(res, 'time_points')
                            data_list(end).time_points = res.time_points;
                        else
                            data_list(end).time_points = res.time_points_absolute;
                        end
                    end
                end
            catch
                fprintf('  读取文件 %s 失败。\n', mat_files(i).name);
            end
        end
        
        if ~isempty(data_list)
            % 按价格排序
            [~, sort_idx] = sort([data_list.price]);
            data_list = data_list(sort_idx);
            
            % 绘图: 聚合总制冷功率对比
            figure('Position', [150 150 1000 600]);
            hold on;
            
            % 颜色映射 (蓝色 -> 红色)
            colors_multi = jet(length(data_list));
            
            legend_str = {};
            for i = 1:length(data_list)
                plot(data_list(i).time_points, data_list(i).total_power, ...
                    'LineWidth', 1.5, 'Color', colors_multi(i,:));
                legend_str{end+1} = sprintf('激励电价 = %.1f 分/kW', data_list(i).price);
            end
            
            hold off;
            xlabel('时间 (小时)', 'FontSize', 12);
            ylabel('聚合总制冷功率 (kW)', 'FontSize', 12);
            % 无标题
            legend(legend_str, 'Location', 'bestoutside', 'FontSize', 10);
            grid on;
            set(gca, 'FontSize', 11);
            xlim([0, 24]);
            
            % 保存
            print(gcf, '图8_不同电价下聚合制冷功率.png', '-dpng', '-r300');
            fprintf('  批量结果绘图完成: 图8_不同电价下聚合制冷功率.png\n');
        else
            fprintf('  未提取到有效的批量数据。\n');
        end
    end
else
    fprintf('提示: 未找到文件夹 "%s"，跳过批量绘图。\n', results_dir);
end

fprintf('\n所有绘图程序执行完毕。\n');