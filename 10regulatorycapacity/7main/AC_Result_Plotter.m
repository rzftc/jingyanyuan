%% AC_Result_Plotter.m
% 功能：读取仿真结果并绘制高DPI、无标题、中文图例的分析图表
%       1. 读取 AC_Stateful_Simulation_Results.mat (单次仿真)
%       2. 读取 results_AC 文件夹下的批量结果 (不同电价)
%       3. [新增] 读取不同 dt (5min, 15min, 60min) 的结果进行对比
%       4. [新增] 绘制激励电价 vs 聚合整体功率特性曲线 (验证死区/饱和区) - 强制包含(0,0)点
% 依赖：AC_main_Stateful_Sim_potential_diff_inc.m 生成的数据

clear; close all; clc;

% 设置默认字体为支持中文的字体 (防止乱码)
set(0, 'DefaultAxesFontName', 'Microsoft YaHei'); 
set(0, 'DefaultTextFontName', 'Microsoft YaHei');

%% === 第一部分：单次仿真结果绘图 ===

mat_filename = 'AC_Stateful_Simulation_Results_5min_pi.mat';
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
        'Agg_Model_Total_Power', 'Agg_Model_Total_Power'; % [新增] 聚合模型功率
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
        
        % --- [修改开始] 筛选第5小时后功率始终小于6kW的空调 ---
        % 1. 找到第5小时之后的时间索引
        idx_after_5h = time_points > 5;
        
        % 2. 找到符合条件的空调列索引 (在idx_after_5h时间段内，最大功率 < 6)
        % Total_Power_History 的每一列代表一台空调
        valid_ac_mask = max(Total_Power_History(idx_after_5h, :), [], 1) < 6;
        
        % 3. 提取需要绘制的数据
        Data_to_Plot = Total_Power_History(:, valid_ac_mask);
        
        fprintf('  图5筛选: 共 %d 台空调，其中 %d 台满足"5小时后功率<6kW"的条件。\n', ...
            size(Total_Power_History, 2), sum(valid_ac_mask));
            
        plot(time_points, Data_to_Plot, 'LineWidth', 0.5);
        % --- [修改结束] ---
        
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
    if ~isempty(Agg_Total_Power)
        figure('Position', [100 600 1000 450]);
        hold on;
        
        % 1. 单体累加的制冷功率 (红色实线)
        plot(time_points, Agg_Total_Power, 'r-', 'LineWidth', 2, 'DisplayName', '聚合总制冷功率 (单体累加)');
        
        % 2. 聚合模型制冷功率 (绿色虚线)
        if ~isempty(Agg_Model_Total_Power)
             plot(time_points, Agg_Model_Total_Power, 'g:', 'LineWidth', 2.5, 'DisplayName', '聚合模型制冷功率');
        end
        
        hold off;
        xlabel('时间 (小时)', 'FontSize', 12);
        ylabel('功率 (kW)', 'FontSize', 12);
        legend('show', 'Location', 'best', 'FontSize', 11);
        set(gca, 'FontSize', 11);
        xlim([0, 24]); grid on;
        print(gcf, '图8_聚合功率对比.png', '-dpng', '-r300');
    end

    fprintf('单次仿真绘图完成。\n');
end

%% === 第二部分：批量结果绘图 (不同电价下的结果对比) ===

results_dir = 'data/results_AC';
fprintf('\n------------------------------------------------------\n');
fprintf('检查批量仿真文件夹 "%s" ...\n', results_dir);

if exist(results_dir, 'dir')
    
    % 搜索所有文件
    file_pattern = fullfile(results_dir, 'AC_Stateful_Simulation_Results_Price_*_pi.mat');
    mat_files = dir(file_pattern);
    
    if isempty(mat_files)
        fprintf('  警告: 在 "%s" 中未找到结果文件。\n', results_dir);
    else
        fprintf('  发现 %d 个结果文件，开始读取...\n', length(mat_files));
        
        % 数据容器 (新增 up_potential, down_potential)
        data_list = struct('price', {}, 'total_power', {}, 'up_potential', {}, 'down_potential', {}, 'time_points', {});
        
        for i = 1:length(mat_files)
            full_path = fullfile(mat_files(i).folder, mat_files(i).name);
            try
                temp_data = load(full_path);
                if isfield(temp_data, 'results')
                    res = temp_data.results;
                    
                    % 获取价格
                    p = NaN;
                    if isfield(res, 'current_price')
                        p = res.current_price;
                    else
                        tokens = regexp(mat_files(i).name, 'Price_([\d\.]+).mat', 'tokens');
                        if ~isempty(tokens), p = str2double(tokens{1}{1}); end
                    end
                    
                    % 获取聚合数据 (总功率 + 上下调节潜力)
                    if ~isnan(p) && isfield(res, 'time_points')
                        data_list(end+1).price = p;
                        
                        % 1. 总功率 (Agg_Total_Power)
                        if isfield(res, 'Agg_Total_Power')
                            data_list(end).total_power = res.Agg_Total_Power;
                        else
                            data_list(end).total_power = [];
                        end
                        
                        % 2. 上调潜力 (Agg_P_Potential_Up_History)
                        if isfield(res, 'Agg_P_Potential_Up_History')
                            data_list(end).up_potential = res.Agg_P_Potential_Up_History;
                        else
                            data_list(end).up_potential = [];
                        end
                        
                        % 3. 下调潜力 (Agg_P_Potential_Down_History)
                        if isfield(res, 'Agg_P_Potential_Down_History')
                            data_list(end).down_potential = res.Agg_P_Potential_Down_History;
                        else
                            data_list(end).down_potential = [];
                        end
                        
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
            
            % 颜色映射 (蓝色 -> 红色)
            colors_multi = jet(length(data_list));
            
            % --- [保留原图] 图 8: 聚合总制冷功率对比 ---
            if ~isempty(data_list(1).total_power)
                figure('Position', [150 150 1000 600]);
                hold on;
                legend_str = {};
                for i = 1:length(data_list)
                    if ~isempty(data_list(i).total_power)
                        plot(data_list(i).time_points, data_list(i).total_power, ...
                            'LineWidth', 1.5, 'Color', colors_multi(i,:));
                        legend_str{end+1} = sprintf('激励电价 = %.1f 分/kW', data_list(i).price);
                    end
                end
                hold off;
                xlabel('时间 (小时)', 'FontSize', 12);
                ylabel('聚合总制冷功率 (kW)', 'FontSize', 12);
                legend(legend_str, 'Location', 'best', 'FontSize', 10);
                grid on; set(gca, 'FontSize', 11); xlim([0, 24]);
                print(gcf, '图8_不同电价下聚合制冷功率.png', '-dpng', '-r300');
                fprintf('  已保存: 图8_不同电价下聚合制冷功率.png\n');
            end

            % --- [新增] 图 9: 聚合上调潜力对比 ---
            if ~isempty(data_list(1).up_potential)
                figure('Position', [200 200 1000 600]);
                hold on;
                legend_str_up = {};
                for i = 1:length(data_list)
                    if ~isempty(data_list(i).up_potential)
                        plot(data_list(i).time_points, data_list(i).up_potential, ...
                            'LineWidth', 1.5, 'Color', colors_multi(i,:));
                        legend_str_up{end+1} = sprintf('激励电价 = %.1f 分/kW', data_list(i).price);
                    end
                end
                hold off;
                xlabel('时间 (小时)', 'FontSize', 12);
                ylabel('AC集群上调潜力 (kW)', 'FontSize', 12);

                legend(legend_str_up, 'Location', 'best', 'FontSize', 10);
                grid on; set(gca, 'FontSize', 11); xlim([0, 24]);
                print(gcf, '图9_不同电价下AC上调能力对比.png', '-dpng', '-r300');
                fprintf('  已保存: 图9_不同电价下AC上调能力对比.png\n');
            end

            % --- [新增] 图 10: 聚合下调潜力对比 ---
            if ~isempty(data_list(1).down_potential)
                figure('Position', [250 250 1000 600]);
                hold on;
                legend_str_down = {};
                for i = 1:length(data_list)
                    if ~isempty(data_list(i).down_potential)
                        plot(data_list(i).time_points, data_list(i).down_potential, ...
                            'LineWidth', 1.5, 'Color', colors_multi(i,:));
                        legend_str_down{end+1} = sprintf('激励电价 = %.1f 分/kW', data_list(i).price);
                    end
                end
                hold off;
                xlabel('时间 (小时)', 'FontSize', 12);
                ylabel('AC集群下调潜力 (kW)', 'FontSize', 12);
   
                legend(legend_str_down, 'Location', 'best', 'FontSize', 10);
                grid on; set(gca, 'FontSize', 11); xlim([0, 24]);
                print(gcf, '图10_不同电价下AC下调能力对比.png', '-dpng', '-r300');
                fprintf('  已保存: 图10_不同电价下AC下调能力对比.png\n');
            end

            % --- [新增] 图 11: 激励价格 vs 聚合整体功率特性曲线 (验证死区/饱和区) ---
            % 说明：提取不同电价下聚合总功率的最大值，验证价格响应特性
            fprintf('  正在绘制图 11 (激励价格-聚合整体功率特性)...\n');
            
            prices_for_curve = [data_list.price];
            max_powers_for_curve = zeros(size(prices_for_curve));
            
            for k = 1:length(data_list)
                if ~isempty(data_list(k).total_power)
                    % 提取该价格下的最大聚合功率 (表征容量)
                    max_powers_for_curve(k) = max(data_list(k).total_power);
                end
            end
            
            % === [修改] 强制添加 (0,0) 点，防止低价跳变 ===
            if ~any(prices_for_curve == 0)
                prices_for_curve = [0, prices_for_curve];
                max_powers_for_curve = [0, max_powers_for_curve];
            end
            
            % 重新排序确保曲线连贯
            [prices_for_curve, sort_idx] = sort(prices_for_curve);
            max_powers_for_curve = max_powers_for_curve(sort_idx);
            % ==================================================
            
            figure('Position', [300 300 800 500]);
            plot(prices_for_curve, max_powers_for_curve, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
            xlabel('激励电价 (分/kW)', 'FontSize', 12);
            ylabel('聚合整体功率峰值 (kW)', 'FontSize', 12);
            grid on; set(gca, 'FontSize', 11);
            
            print(gcf, '图11_AC激励价格-聚合整体功率特性.png', '-dpng', '-r300');
            fprintf('  已保存: 图11_AC激励价格-聚合整体功率特性.png\n');

        else
            fprintf('  未提取到有效的批量数据。\n');
        end
    end
else
    fprintf('提示: 未找到文件夹 "%s"，跳过批量绘图。\n', results_dir);
end
%% === 第三部分：不同时间步长 (dt) 对比绘图 (新增) ===
% 功能：加载 5min, 15min, 60min 的仿真结果并对比上下调节能力

fprintf('\n------------------------------------------------------\n');
fprintf('=== 开始执行不同时间步长 (dt) 对比绘图 ===\n');

% 1. 定义文件和标签
dt_files = {
    'AC_Stateful_Simulation_Results_5min.mat', ...
    'AC_Stateful_Simulation_Results_15min.mat', ...
    'AC_Stateful_Simulation_Results_60min.mat'
};
dt_labels = {'5 min', '15 min', '60 min'};
line_styles = {'-', '--', '-.'}; % 使用不同线型区分
colors_dt = lines(3); % 使用 distinct 颜色

% 2. 数据加载
data_dt_list = struct('label', {}, 'up', {}, 'down', {}, 'time', {});
valid_dt_count = 0;

for i = 1:length(dt_files)
    fname = dt_files{i};
    if exist(fname, 'file')
        fprintf('  加载文件: %s\n', fname);
        try
            tmp = load(fname);
            if isfield(tmp, 'results')
                res = tmp.results;
                
                % 检查必要字段
                if isfield(res, 'Agg_P_Potential_Up_History') && isfield(res, 'Agg_P_Potential_Down_History')
                    valid_dt_count = valid_dt_count + 1;
                    data_dt_list(valid_dt_count).label = dt_labels{i};
                    data_dt_list(valid_dt_count).up = res.Agg_P_Potential_Up_History;
                    data_dt_list(valid_dt_count).down = res.Agg_P_Potential_Down_History;
                    
                    % 获取时间轴
                    if isfield(res, 'time_points')
                        data_dt_list(valid_dt_count).time = res.time_points;
                    elseif isfield(res, 'time_points_absolute')
                        data_dt_list(valid_dt_count).time = res.time_points_absolute;
                    else
                        % 如果没有时间轴，根据数据长度生成默认的
                        len = length(res.Agg_P_Potential_Up_History);
                        data_dt_list(valid_dt_count).time = linspace(0, 24, len);
                    end
                end
            else
                warning('文件 %s 中未找到 "results" 结构体。', fname);
            end
        catch
            warning('读取文件 %s 失败。', fname);
        end
    else
        fprintf('  [警告] 文件不存在: %s，已跳过。\n', fname);
    end
end

if valid_dt_count > 0
    % 3. 绘图：上调能力对比
    figure('Name', 'AC集群不同dt上调能力对比', 'Position', [200 200 1000 500]);
    hold on; grid on;
    for i = 1:valid_dt_count
        plot(data_dt_list(i).time, data_dt_list(i).up, ...
            'LineWidth', 2.0, ...
            'LineStyle', line_styles{i}, ...
            'Color', colors_dt(i,:), ...
            'DisplayName', ['dt = ' data_dt_list(i).label]);
    end
    hold off;
    xlabel('时间 (小时)', 'FontSize', 14);
    ylabel('AC集群上调潜力 (kW)', 'FontSize', 14);
    
    legend('show', 'Location', 'best', 'FontSize', 12);
    set(gca, 'FontSize', 12);
    xlim([0, 24]); 
    print(gcf, '图9_AC_Cluster_Up_Comparison_dt.png', '-dpng', '-r300');
    fprintf('  上调对比图已保存为: 图9_AC_Cluster_Up_Comparison_dt.png\n');
    
    % 4. 绘图：下调能力对比
    figure('Name', 'AC集群不同dt下调能力对比', 'Position', [200 750 1000 500]);
    hold on; grid on;
    for i = 1:valid_dt_count
        plot(data_dt_list(i).time, data_dt_list(i).down, ...
            'LineWidth', 2.0, ...
            'LineStyle', line_styles{i}, ...
            'Color', colors_dt(i,:), ...
            'DisplayName', ['dt = ' data_dt_list(i).label]);
    end
    hold off;
    xlabel('时间 (小时)', 'FontSize', 14);
    ylabel('AC集群下调潜力 (kW)', 'FontSize', 14);
   
    legend('show', 'Location', 'best', 'FontSize', 12);
    set(gca, 'FontSize', 12);
    xlim([0, 24]);
    print(gcf, '图10_AC_Cluster_Down_Comparison_dt.png', '-dpng', '-r300');
    fprintf('  下调对比图已保存为: 图10_AC_Cluster_Down_Comparison_dt.png\n');
else
    fprintf('  未加载到任何有效数据，无法绘制dt对比图。\n');
end

fprintf('\n所有绘图程序执行完毕。\n');
close all;