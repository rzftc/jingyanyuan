

clear; close all; clc;

% [已修改] 恢复默认字体设置
% set(0, 'DefaultAxesFontName', 'Microsoft YaHei'); 
% set(0, 'DefaultTextFontName', 'Microsoft YaHei');

%% === 第一部分：单次仿真结果绘图 ===

mat_filename = 'AC_Stateful_Simulation_Results_5min_pi_8am.mat';
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
    
    % [修改] 确保时间轴是 8 到 32，如果数据直接是 8-32，无需处理
    if isrow(time_points), time_points = time_points'; end

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
        'Agg_Model_Total_Power', 'Agg_Model_Total_Power'; 
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

    % [修改] 已移除数据拼接逻辑，直接使用读取的数据

    % --- 2. 开始绘图 (图1-7) ---
    fprintf('正在生成单次仿真图表 (8:00 - 次日 8:00)...\n');
    
    % [修改] 定义通用 X 轴设置函数：范围 8-32，但标签显示为时刻，避免显示 "32"
    set_xaxis_custom = @() set(gca, ...
        'XLim', [8, 32], ...
        'XTick', 8:4:32, ...
        'XTickLabel', {'08:00','12:00','16:00','20:00','D2:00:00','D2:04:00','D2:08:00'}, ...
        'FontSize', 16);
    
    % 图 1: 功率跟踪对比
    if ~isempty(Agg_P_Command_History) && ~isempty(Agg_P_Achieved_History)
        figure('Position', [100 100 1000 450]);
        plot(time_points, Agg_P_Command_History, 'k:', 'LineWidth', 2.5, 'DisplayName', '电网调节指令');
        hold on;
        plot(time_points, Agg_P_Achieved_History, 'r-', 'LineWidth', 1.5, 'DisplayName', '聚合响应功率');
        hold off;
        % 字体放大
        xlabel('时间', 'FontSize', 20);
        ylabel('功率 (kW)', 'FontSize', 20);
        legend('show', 'Location', 'best', 'FontSize', 16);
        set_xaxis_custom(); grid on;
        print(gcf, '图1_功率跟踪对比.png', '-dpng', '-r300');
    end
    
    % 图 2: SOC状态对比
    if ~isempty(Individual_SOC_History) && ~isempty(Agg_SOC_History)
        figure('Position', [100 550 1000 450]);
        hold on;
        h_ind = plot(time_points, Individual_SOC_History, 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
        h_agg = plot(time_points, Agg_SOC_History, 'k--', 'LineWidth', 2.5, 'DisplayName', '聚合SOC (均值)');
        if num_AC_participating > 0
            set(h_ind(1), 'DisplayName', '单体空调SOC');
            if num_AC_participating > 1, set(h_ind(2:end), 'HandleVisibility', 'off'); end
            legend([h_agg, h_ind(1)], 'Location', 'best', 'FontSize', 16); % 字体放大
        else
            legend(h_agg, 'Location', 'best', 'FontSize', 16); % 字体放大
        end
        hold off;
        % 字体放大
        xlabel('时间', 'FontSize', 20);
        ylabel('SOC', 'FontSize', 20);
        set_xaxis_custom();
        ylim([-0.1, 1.1]); grid on;
        print(gcf, '图2_SOC状态对比.png', '-dpng', '-r300');
    end
    
    % 图 3: 室内温度变化
    if ~isempty(Individual_Temp_History)
        figure('Position', [100 300 1000 450]);
        plot(time_points, Individual_Temp_History, 'LineWidth', 0.5);
        % 字体放大
        xlabel('时间', 'FontSize', 20);
        ylabel('温度 (°C)', 'FontSize', 20);
        set_xaxis_custom(); grid on;
        print(gcf, '图3_室内温度变化.png', '-dpng', '-r300');
    end
    
    % 图 4: 单体调节功率
    if ~isempty(Individual_Power_History)
        figure('Position', [100 400 1000 450]);
        plot(time_points, Individual_Power_History, 'LineWidth', 0.5);
        yline(0, 'k--', 'LineWidth', 1.5);
        % 字体放大
        xlabel('时间', 'FontSize', 20);
        ylabel('功率 (kW)', 'FontSize', 20);
        set_xaxis_custom(); grid on;
        print(gcf, '图4_单体调节功率.png', '-dpng', '-r300');
    end
    
   % 图 5: 单体总制冷功率
    if ~isempty(Total_Power_History)
        figure('Position', [100 500 1000 450]);
        plot(time_points, Total_Power_History, 'LineWidth', 0.5);
        yline(P_standby, 'k--', 'LineWidth', 1.5, 'DisplayName', '待机功率');
        % 字体放大
        xlabel('时间', 'FontSize', 20);
        ylabel('功率 (kW)', 'FontSize', 20);
        set_xaxis_custom(); grid on;
        print(gcf, '图5_单体总制冷功率.png', '-dpng', '-r300');
    end
    
    % 图 6: 聚合功率对比
    if ~isempty(Agg_Baseline_Power) && ~isempty(Agg_Total_Power)
        figure('Position', [100 600 1000 450]);
        hold on;
        plot(time_points, Agg_Baseline_Power, 'b--', 'LineWidth', 2, 'DisplayName', '聚合基线功率');
        plot(time_points, Agg_Total_Power, 'r-', 'LineWidth', 2, 'DisplayName', '聚合总制冷功率');
        hold off;
        % 字体放大
        xlabel('时间', 'FontSize', 20);
        ylabel('功率 (kW)', 'FontSize', 20);
        legend('show', 'Location', 'best', 'FontSize', 16);
        set_xaxis_custom(); grid on;
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
        % 字体放大
        xlabel('时间', 'FontSize', 20);
        ylabel('功率 (kW)', 'FontSize', 20);
        legend('show', 'Location', 'best', 'FontSize', 16);
        set_xaxis_custom(); grid on;
        print(gcf, '图7_聚合潜力对比.png', '-dpng', '-r300');
    end

    % =====================================================================
    % 图 8: 聚合功率对比 (包含局部放大子图 - EMF矢量版)
    % =====================================================================
    if ~isempty(Agg_Total_Power)
        
        if ~isempty(Agg_Model_Total_Power)
            % --- 1. 计算偏差并定位放大区域 ---
            Agg_Total_Vec = Agg_Total_Power(:);
            Agg_Model_Vec = Agg_Model_Total_Power(:);
            Time_Vec = time_points(:);
            
            diff_curve = abs(Agg_Total_Vec - Agg_Model_Vec);
            [max_diff, idx_max] = max(diff_curve);
            
            % 定义放大窗口参数
            t_center = Time_Vec(idx_max);
            window_width = 0.2; % 小时
            % [修改] 边界限制调整为 8 到 32
            t_start_zoom = max(8, t_center - window_width/2);
            t_end_zoom = min(32, t_center + window_width/2);
            
            % 获取窗口内的数据索引
            mask_zoom = (Time_Vec >= t_start_zoom) & (Time_Vec <= t_end_zoom);
            
            % 计算Y轴显示范围
            data_in_window = [Agg_Total_Vec(mask_zoom); Agg_Model_Vec(mask_zoom)];
            if isempty(data_in_window)
                 % 容错处理
                 y_rect_min = 0; y_rect_max = 1; y_rect_h = 1;
            else
                y_min_zoom = min(data_in_window);
                y_max_zoom = max(data_in_window);
                y_margin = (y_max_zoom - y_min_zoom) * 0.1; 
                if y_margin == 0, y_margin = 1; end
                y_rect_min = y_min_zoom - y_margin;
                y_rect_max = y_max_zoom + y_margin;
                y_rect_h = y_rect_max - y_rect_min;
            end
            
            % --- 2. 绘制母图 (白色背景，保存为PNG) ---
            figure('Position', [100 600 1000 450]);
            hold on;
            plot(time_points, Agg_Total_Power, 'r-', 'LineWidth', 2, 'DisplayName', '聚合总制冷功率 (单体累加)');
            plot(time_points, Agg_Model_Total_Power, 'g:', 'LineWidth', 2.5, 'DisplayName', '聚合模型制冷功率');
            
            % 绘制标注矩形框
            rectangle('Position', [t_start_zoom, y_rect_min, window_width, y_rect_h], ...
                      'EdgeColor', 'k', 'LineWidth', 1.5, 'LineStyle', '--');
            
            hold off;
            % 字体放大
            xlabel('时间', 'FontSize', 20);
            ylabel('功率 (kW)', 'FontSize', 20);
            legend('show', 'Location', 'best', 'FontSize', 16);
            set_xaxis_custom(); grid on;
            
            print(gcf, '图8_聚合功率对比_母图.png', '-dpng', '-r300');
            fprintf('  已保存: 图8_聚合功率对比_母图.png\n');
            
            % --- 3. 绘制子图 (EMF矢量格式，支持透明) ---
            figure('Position', [150 650 500 350]); % 尺寸稍小
            
            % 在屏幕上保持白色背景，方便查看
            set(gcf, 'Color', 'w');
            
            hold on;
            if any(mask_zoom)
                plot(time_points(mask_zoom), Agg_Total_Power(mask_zoom), 'r-', 'LineWidth', 2);
                plot(time_points(mask_zoom), Agg_Model_Total_Power(mask_zoom), 'g:', 'LineWidth', 2.5);
            end
            hold off;
            
            % 设置坐标轴范围
            xlim([t_start_zoom, t_end_zoom]);
            ylim([y_rect_min, y_rect_max]);
            grid off; % <--- 移除子图网格线，保留坐标轴
            
            % 子图字体特大
            set(gca, 'FontSize', 24); 
            
            % 去掉标题、XY轴标签、图例
            xlabel(''); ylabel('');
            legend('off');
            
            % === 关键修改：保存为 EMF 矢量图 ===
            img_filename = '图8_聚合功率对比_子图.emf';
            
            % 设置属性以确保 EMF 转换时不带背景色
            set(gcf, 'Color', 'none'); 
            set(gca, 'Color', 'none');
            set(gcf, 'InvertHardcopy', 'off'); 
            
            % 使用 -dmeta 指令保存为增强型图元文件 (EMF)
            print(gcf, img_filename, '-dmeta');
            
            % 恢复白色背景
            set(gcf, 'Color', 'w');
            set(gca, 'Color', 'w');
            
            fprintf('  已保存: %s (EMF矢量格式，支持透明)\n', img_filename);
        end
    end

    fprintf('单次仿真绘图完成。\n');
end

%% === 第二部分：批量结果绘图 (不同电价下的结果对比) ===

results_dir = 'data/results_AC';
fprintf('\n------------------------------------------------------\n');
fprintf('检查批量仿真文件夹 "%s" ...\n', results_dir);

if exist(results_dir, 'dir')
    file_pattern = fullfile(results_dir, 'AC_Stateful_Simulation_Results_Price_*_pi_8am.mat');
    mat_files = dir(file_pattern);
    
    if isempty(mat_files)
        fprintf('  警告: 在 "%s" 中未找到结果文件。\n', results_dir);
    else
        fprintf('  发现 %d 个结果文件，开始读取...\n', length(mat_files));
        
        data_list = struct('price', {}, 'total_power', {}, 'up_potential', {}, 'down_potential', {}, 'time_points', {});
        
        for i = 1:length(mat_files)
            full_path = fullfile(mat_files(i).folder, mat_files(i).name);
            try
                temp_data = load(full_path);
                if isfield(temp_data, 'results')
                    res = temp_data.results;
                    p = NaN;
                    if isfield(res, 'current_price')
                        p = res.current_price;
                    else
                        tokens = regexp(mat_files(i).name, 'Price_([\d\.]+).mat', 'tokens');
                        if ~isempty(tokens), p = str2double(tokens{1}{1}); end
                    end
                    
                    if ~isnan(p) && isfield(res, 'time_points')
                        % [修改] 批量数据也假定已经是 8-32 的格式
                        tp = res.time_points;
                        if isfield(res, 'time_points_absolute'), tp = res.time_points_absolute; end
                        if isrow(tp), tp = tp'; end
                        
                        data_list(end+1).price = p;
                        data_list(end).time_points = tp;
                        
                        if isfield(res, 'Agg_Total_Power'), data_list(end).total_power = res.Agg_Total_Power; else, data_list(end).total_power = []; end
                        if isfield(res, 'Agg_P_Potential_Up_History'), data_list(end).up_potential = res.Agg_P_Potential_Up_History; else, data_list(end).up_potential = []; end
                        if isfield(res, 'Agg_P_Potential_Down_History'), data_list(end).down_potential = res.Agg_P_Potential_Down_History; else, data_list(end).down_potential = []; end
                    end
                end
            catch
                fprintf('  读取文件 %s 失败。\n', mat_files(i).name);
            end
        end
        
        if ~isempty(data_list)
            [~, sort_idx] = sort([data_list.price]);
            data_list = data_list(sort_idx);
            colors_multi = jet(length(data_list));
            
            % [修改] 重新定义 X 轴设置函数
            set_xaxis_custom = @() set(gca, ...
                'XLim', [8, 32], ...
                'XTick', 8:4:32, ...
                'XTickLabel', {'08:00','12:00','16:00','20:00','D2:00:00','D2:04:00','D2:08:00'}, ...
                'FontSize', 16);

            % 图 8: 聚合总制冷功率对比
            if ~isempty(data_list(1).total_power)
                figure('Position', [150 150 1000 600]);
                hold on;
                legend_str = {};
                for i = 1:length(data_list)
                    if ~isempty(data_list(i).total_power)
                        plot(data_list(i).time_points, data_list(i).total_power, 'LineWidth', 1.5, 'Color', colors_multi(i,:));
                        legend_str{end+1} = sprintf('%.2f 元/kW', data_list(i).price / 100);
                    end
                end
                hold off;
                % 字体放大
                xlabel('时间', 'FontSize', 20);
                ylabel('聚合总制冷功率 (kW)', 'FontSize', 20);
                legend(legend_str, 'Location', 'best', 'FontSize', 16);
                grid on; set_xaxis_custom();
                print(gcf, '图8_不同电价下聚合制冷功率.png', '-dpng', '-r300');
            end

            % 图 9: 聚合上调潜力对比
            if ~isempty(data_list(1).up_potential)
                figure('Position', [200 200 1000 600]);
                hold on;
                legend_str_up = {};
                for i = 1:length(data_list)
                    if ~isempty(data_list(i).up_potential)
                        plot(data_list(i).time_points, data_list(i).up_potential, 'LineWidth', 1.5, 'Color', colors_multi(i,:));
                        legend_str_up{end+1} = sprintf('%.2f 元/kW', data_list(i).price / 100);
                    end
                end
                hold off;
                % 字体放大
                xlabel('时间', 'FontSize', 20);
                ylabel('AC集群上调潜力 (kW)', 'FontSize', 20);
                legend(legend_str_up, 'Location', 'best', 'FontSize', 16);
                grid on; set_xaxis_custom();
                print(gcf, '图9_不同电价下AC上调能力对比.png', '-dpng', '-r300');
            end

            % 图 10: 聚合下调潜力对比
            if ~isempty(data_list(1).down_potential)
                figure('Position', [250 250 1000 600]);
                hold on;
                legend_str_down = {};
                for i = 1:length(data_list)
                    if ~isempty(data_list(i).down_potential)
                        plot(data_list(i).time_points, data_list(i).down_potential, 'LineWidth', 1.5, 'Color', colors_multi(i,:));
                        legend_str_down{end+1} = sprintf('%.2f 元/kW', data_list(i).price / 100);
                    end
                end
                hold off;
                % 字体放大
                xlabel('时间', 'FontSize', 20);
                ylabel('AC集群下调潜力 (kW)', 'FontSize', 20);
                legend(legend_str_down, 'Location', 'best', 'FontSize', 16);
                grid on; set_xaxis_custom();
                print(gcf, '图10_不同电价下AC下调能力对比.png', '-dpng', '-r300');
            end

            % 图 11: 激励价格 vs 聚合整体功率特性曲线
            % 注意：此图 X 轴是价格，无需调整
            fprintf('  正在绘制图 11 (激励价格-聚合整体功率特性)...\n');
            prices_for_curve = [data_list.price];
            max_powers_for_curve = zeros(size(prices_for_curve));
            for k = 1:length(data_list)
                if ~isempty(data_list(k).total_power)
                    max_powers_for_curve(k) = max(data_list(k).total_power);
                end
            end
            if ~any(prices_for_curve == 0)
                prices_for_curve = [0, prices_for_curve];
                max_powers_for_curve = [0, max_powers_for_curve];
            end
            [prices_for_curve, sort_idx] = sort(prices_for_curve);
            max_powers_for_curve = max_powers_for_curve(sort_idx);
            
            figure('Position', [300 300 800 500]);
            plot(prices_for_curve / 100, max_powers_for_curve, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
            % 字体放大
            xlabel('激励电价 (元/kW)', 'FontSize', 20);
            ylabel('聚合整体功率峰值 (kW)', 'FontSize', 20);
            grid on; set(gca, 'FontSize', 16);
            print(gcf, '图11_AC激励价格-聚合整体功率特性.png', '-dpng', '-r300');
        else
            fprintf('  未提取到有效的批量数据。\n');
        end
    end
else
    fprintf('提示: 未找到文件夹 "%s"，跳过批量绘图。\n', results_dir);
end

%% === 第三部分：不同时间步长 (dt) 对比绘图 ===

fprintf('\n------------------------------------------------------\n');
fprintf('=== 开始执行不同时间步长 (dt) 对比绘图 (8-32h) ===\n');

dt_files = {
    'AC_Stateful_Simulation_Results_5min_pi_8am.mat', ...
    'AC_Stateful_Simulation_Results_15min_pi_8am.mat', ...
    'AC_Stateful_Simulation_Results_60min_pi_8am.mat'
};
dt_labels = {'5 min', '15 min', '60 min'};
line_styles = {'-', '--', '-.'};
colors_dt = lines(3);

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
                if isfield(res, 'Agg_P_Potential_Up_History') && isfield(res, 'Agg_P_Potential_Down_History')
                    valid_dt_count = valid_dt_count + 1;
                    
                    % [修改] 假定数据已经是 8-32 格式，不再拼接
                    if isfield(res, 'time_points')
                        tp = res.time_points;
                    elseif isfield(res, 'time_points_absolute')
                        tp = res.time_points_absolute;
                    else
                        len = length(res.Agg_P_Potential_Up_History);
                        tp = linspace(8, 32, len)'; % 假定也是8-32
                    end
                    if isrow(tp), tp = tp'; end
                    
                    data_dt_list(valid_dt_count).label = dt_labels{i};
                    data_dt_list(valid_dt_count).up = res.Agg_P_Potential_Up_History;
                    data_dt_list(valid_dt_count).down = res.Agg_P_Potential_Down_History;
                    data_dt_list(valid_dt_count).time = tp;
                end
            end
        catch
            warning('读取文件 %s 失败。', fname);
        end
    end
end

if valid_dt_count > 0
    set_xaxis_custom = @() set(gca, ...
        'XLim', [8, 32], ...
        'XTick', 8:4:32, ...
        'XTickLabel', {'08:00','12:00','16:00','20:00','D2:00:00','D2:04:00','D2:08:00'}, ...
        'FontSize', 16);
    
    figure('Name', 'AC集群不同dt上调能力对比', 'Position', [200 200 1000 500]);
    hold on; grid on;
    for i = 1:valid_dt_count
        plot(data_dt_list(i).time, data_dt_list(i).up, 'LineWidth', 2.0, 'LineStyle', line_styles{i}, 'Color', colors_dt(i,:), 'DisplayName', ['dt = ' data_dt_list(i).label]);
    end
    hold off;
    % 字体放大
    xlabel('时间', 'FontSize', 20);
    ylabel('AC集群上调潜力 (kW)', 'FontSize', 20);
    legend('show', 'Location', 'best', 'FontSize', 16);
    set_xaxis_custom();
    print(gcf, '图9_AC_Cluster_Up_Comparison_dt.png', '-dpng', '-r300');
    
    figure('Name', 'AC集群不同dt下调能力对比', 'Position', [200 750 1000 500]);
    hold on; grid on;
    for i = 1:valid_dt_count
        plot(data_dt_list(i).time, data_dt_list(i).down, 'LineWidth', 2.0, 'LineStyle', line_styles{i}, 'Color', colors_dt(i,:), 'DisplayName', ['dt = ' data_dt_list(i).label]);
    end
    hold off;
    % 字体放大
    xlabel('时间', 'FontSize', 20);
    ylabel('AC集群下调潜力 (kW)', 'FontSize', 20);
    legend('show', 'Location', 'best', 'FontSize', 16);
    set_xaxis_custom();
    print(gcf, '图10_AC_Cluster_Down_Comparison_dt.png', '-dpng', '-r300');
else
    fprintf('  未加载到任何有效数据，无法绘制dt对比图。\n');
end

fprintf('\n所有绘图程序执行完毕。\n');
close all;