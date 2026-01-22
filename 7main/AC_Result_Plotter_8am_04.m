%% AC_Result_Plotter_8am.m
% 功能：读取仿真结果并绘制高DPI、无标题、中文图例的分析图表 (8:00起)
%       1. 读取 AC_Stateful_Simulation_Results_..._8am.mat (单次仿真)
%       2. 读取 results_AC 文件夹下的批量结果 (不同电价)
%       3. 读取不同 dt (5min, 15min, 60min) 的结果进行对比
%       4. 绘制激励电价 vs 聚合整体功率特性曲线 (验证死区/饱和区)
%       5. [修复] 图8子图保存为 EMF 矢量格式
%       6. [修改 V4] 图2 SOC 矩形框高度自适应数据范围，仅包裹分散区域
% 依赖：AC_main_1_inc_pi_8am.m 生成的数据

clear; close all; clc;

%% === 第一部分：单次仿真结果绘图 ===

mat_filename = 'AC_Stateful_Simulation_Results_5min_pi_8am_04.mat';
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
    
    % 确保时间轴是 8 到 32
    if isrow(time_points), time_points = time_points'; end

    % 提取数据
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

    % --- 2. 开始绘图 (图1-7) ---
    fprintf('正在生成单次仿真图表 (8:00 - 次日 8:00)...\n');
    
    % 定义通用 X 轴设置函数
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
        xlabel('时间', 'FontSize', 20);
        ylabel('功率 (kW)', 'FontSize', 20);
        legend('show', 'Location', 'best', 'FontSize', 16);
        set_xaxis_custom(); grid on;
        print(gcf, '图1_功率跟踪对比.png', '-dpng', '-r300');
    end
    
    % 图 2: SOC状态对比 (含彩色线条、自适应矩形框及子图)
    if ~isempty(Individual_SOC_History) && ~isempty(Agg_SOC_History)
        fig2 = figure('Position', [100 550 1000 450]);
        ax2 = gca;
        hold(ax2, 'on');
        
        % --- 彩色线条绘制 ---
        h_ind = plot(time_points, Individual_SOC_History, 'LineWidth', 0.5);
        num_lines = length(h_ind);
        line_colors = jet(num_lines); 
        for k = 1:num_lines
            set(h_ind(k), 'Color', line_colors(k, :));
        end

        h_agg = plot(time_points, Agg_SOC_History, 'k--', 'LineWidth', 2.5, 'DisplayName', '聚合SOC (均值)');
        
        if num_AC_participating > 0
            set(h_ind(1), 'DisplayName', '单体空调SOC');
            if num_AC_participating > 1, set(h_ind(2:end), 'HandleVisibility', 'off'); end
            legend([h_agg, h_ind(1)], 'Location', 'best', 'FontSize', 16);
        else
            legend(h_agg, 'Location', 'best', 'FontSize', 16);
        end
        
        % 基础格式设置
        xlabel('时间', 'FontSize', 20);
        ylabel('SOC', 'FontSize', 20);
        set_xaxis_custom();
        ylim([-0.1, 1.1]); grid on;

        % --- [修改]：寻找SOC最分散区域并绘制自适应高度的矩形 ---
        % 1. 计算每一时刻的SOC标准差
        soc_std = std(Individual_SOC_History, 0, 2);
        
        % 2. 找到标准差最大的时刻
        [max_std, idx_max_spread] = max(soc_std);
        t_center = time_points(idx_max_spread);
        
        % 3. 定义放大窗口 (30分钟)
        zoom_window_hours = 0.5; 
        t_zoom_start = max(8, t_center - zoom_window_hours/2);
        t_zoom_end = min(32, t_center + zoom_window_hours/2);
        
        % 4. 计算该窗口内的SOC最大/最小值，确定矩形高度
        mask_zoom = (time_points >= t_zoom_start) & (time_points <= t_zoom_end);
        if any(mask_zoom)
            data_in_window = Individual_SOC_History(mask_zoom, :);
            soc_min_zoom = min(data_in_window, [], 'all');
            soc_max_zoom = max(data_in_window, [], 'all');
            
            % 添加一点视觉边距 (例如上下各加 0.05 SOC)
            rect_y = max(-0.1, soc_min_zoom - 0.05);
            rect_top = min(1.1, soc_max_zoom + 0.05);
            rect_height = rect_top - rect_y;
        else
            % 兜底逻辑
            rect_y = 0; rect_height = 1;
        end

        fprintf('  [图2] SOC最分散时刻: %.2f (Std: %.4f), 窗口: %.2f-%.2f, SOC范围: %.2f-%.2f\n', ...
            t_center, max_std, t_zoom_start, t_zoom_end, rect_y, rect_y+rect_height);

        % 5. 在母图上绘制黑色虚线框 (高度紧凑)
        rect_width = t_zoom_end - t_zoom_start;
        rectangle(ax2, 'Position', [t_zoom_start, rect_y, rect_width, rect_height], ...
                  'EdgeColor', 'k', 'LineWidth', 2, 'LineStyle', '--');
       
        hold(ax2, 'off');
        print(fig2, '图2_SOC状态对比.png', '-dpng', '-r300');
        
        % --- [修改]：绘制局部放大子图 (去边框/网格/图例) ---
        fig2_sub = figure('Name', 'SOC子图', 'Color', 'w', 'Position', [150 600 600 400]);
        ax2_sub = axes(fig2_sub);
        hold(ax2_sub, 'on');
        
        t_sub = time_points(mask_zoom);
        data_sub = Individual_SOC_History(mask_zoom, :);
        agg_sub = Agg_SOC_History(mask_zoom);
        
        % 绘制彩色线条
        for k = 1:num_lines
            plot(ax2_sub, t_sub, data_sub(:, k), 'LineWidth', 1.0, 'Color', line_colors(k, :));
        end
        plot(ax2_sub, t_sub, agg_sub, 'k--', 'LineWidth', 3.0); % 加粗聚合线
        
        % 子图格式设置
        xlim(ax2_sub, [t_zoom_start, t_zoom_end]);
        % Y轴范围略大于数据范围
        ylim(ax2_sub, [rect_y, rect_y + rect_height]);
        % 
        % xlabel(ax2_sub, '时间', 'FontSize', 24, 'FontWeight', 'bold');
        % ylabel(ax2_sub, 'SOC', 'FontSize', 24, 'FontWeight', 'bold');
        
        % 去掉上边框和右边框，刻度朝外
        set(ax2_sub, 'FontSize', 26, 'Box', 'off', 'LineWidth', 1.5, 'TickDir', 'out');
        
        % 去掉网格
        grid(ax2_sub, 'off');
        
        % 确保无图例
        legend(ax2_sub, 'off');
        
        % 保存子图为 EMF (透明背景)
        set(fig2_sub, 'Color', 'none'); 
        set(ax2_sub, 'Color', 'none');
        set(fig2_sub, 'InvertHardcopy', 'off'); 
        
        print(fig2_sub, '图2_SOC状态对比_子图.emf', '-dmeta');
        fprintf('  已保存: 图2_SOC状态对比.png (母图) 和 图2_SOC状态对比_子图.emf (紧凑矢量子图)\n');
        
        set(fig2_sub, 'Color', 'w'); % 恢复白色
    end
    
    % 图 3: 室内温度变化
    if ~isempty(Individual_Temp_History)
        figure('Position', [100 300 1000 450]);
        plot(time_points, Individual_Temp_History, 'LineWidth', 0.5);
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
            t_start_zoom = max(8, t_center - window_width/2);
            t_end_zoom = min(32, t_center + window_width/2);
            
            % 获取窗口内的数据索引
            mask_zoom = (Time_Vec >= t_start_zoom) & (Time_Vec <= t_end_zoom);
            
            % 计算Y轴显示范围
            data_in_window = [Agg_Total_Vec(mask_zoom); Agg_Model_Vec(mask_zoom)];
            if isempty(data_in_window)
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
            
            % --- 2. 绘制母图 ---
            figure('Position', [100 600 1000 450]);
            hold on;
            plot(time_points, Agg_Total_Power, 'r-', 'LineWidth', 2, 'DisplayName', '聚合总制冷功率 (单体累加)');
            plot(time_points, Agg_Model_Total_Power, 'g:', 'LineWidth', 2.5, 'DisplayName', '聚合模型制冷功率');
            
            rectangle('Position', [t_start_zoom, y_rect_min, window_width, y_rect_h], ...
                      'EdgeColor', 'k', 'LineWidth', 1.5, 'LineStyle', '--');
            
            hold off;
            xlabel('时间', 'FontSize', 20);
            ylabel('功率 (kW)', 'FontSize', 20);
            legend('show', 'Location', 'best', 'FontSize', 16);
            set_xaxis_custom(); grid on;
            
            print(gcf, '图8_聚合功率对比_母图.png', '-dpng', '-r300');
            
            % --- 3. 绘制子图 ---
            figure('Position', [150 650 500 350]); 
            set(gcf, 'Color', 'w');
            
            hold on;
            if any(mask_zoom)
                plot(time_points(mask_zoom), Agg_Total_Power(mask_zoom), 'r-', 'LineWidth', 2);
                plot(time_points(mask_zoom), Agg_Model_Total_Power(mask_zoom), 'g:', 'LineWidth', 2.5);
            end
            hold off;
            
            xlim([t_start_zoom, t_end_zoom]);
            ylim([y_rect_min, y_rect_max]);
            grid off; 
            
            set(gca, 'FontSize', 24); 
            xlabel(''); ylabel('');
            legend('off');
            
            img_filename = '图8_聚合功率对比_子图.emf';
            set(gcf, 'Color', 'none'); 
            set(gca, 'Color', 'none');
            set(gcf, 'InvertHardcopy', 'off'); 
            
            print(gcf, img_filename, '-dmeta');
            set(gcf, 'Color', 'w'); set(gca, 'Color', 'w');
            fprintf('  已保存: %s\n', img_filename);
        end
    end

    fprintf('单次仿真绘图完成。\n');
end

%% === 第二部分：批量结果绘图 (不同电价下的结果对比) ===

results_dir = 'data/results_AC';
fprintf('\n------------------------------------------------------\n');
fprintf('检查批量仿真文件夹 "%s" ...\n', results_dir);

if exist(results_dir, 'dir')
    file_pattern = fullfile(results_dir, 'AC_Stateful_Simulation_Results_Price_*_pi_8am_04.mat');
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
                        legend_str{end+1} = sprintf('%.2f 元/kWh', data_list(i).price / 100);
                    end
                end
                hold off;
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
                        legend_str_up{end+1} = sprintf('%.2f 元/kWh', data_list(i).price / 100);
                    end
                end
                hold off;
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
                        legend_str_down{end+1} = sprintf('%.2f 元/kWh', data_list(i).price / 100);
                    end
                end
                hold off;
                xlabel('时间', 'FontSize', 20);
                ylabel('AC集群下调潜力 (kW)', 'FontSize', 20);
                legend(legend_str_down, 'Location', 'best', 'FontSize', 16);
                grid on; set_xaxis_custom();
                print(gcf, '图10_不同电价下AC下调能力对比.png', '-dpng', '-r300');
            end

            % 图 11: 激励价格 vs 聚合整体功率特性曲线 (样式修改版)
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
            
            % [修改]: 样式调整为与 EV 代码一致
            fig11 = figure('Name', 'AC激励价格-聚合整体功率特性', 'Position', [300 300 800 500], 'NumberTitle', 'off');
            plot(prices_for_curve / 100, max_powers_for_curve, 'bo-', 'LineWidth', 2.5, 'MarkerSize', 12, 'MarkerFaceColor', 'b');
            xlabel('激励电价 (元/kWh)', 'FontSize', 23);
            ylabel('聚合整体功率峰值 (kW)', 'FontSize', 23);
            set(gca, 'FontSize', 19);
            grid on;
            set(fig11, 'Renderer', 'painters');
            print(fig11, '图11_AC激励价格-聚合整体功率特性.png', '-dpng', '-r600');
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
    'AC_Stateful_Simulation_Results_5min_pi_8am_04.mat', ...
    'AC_Stateful_Simulation_Results_15min_pi_8am_04.mat', ...
    'AC_Stateful_Simulation_Results_60min_pi_8am_04.mat'
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
                    
                    if isfield(res, 'time_points')
                        tp = res.time_points;
                    elseif isfield(res, 'time_points_absolute')
                        tp = res.time_points_absolute;
                    else
                        len = length(res.Agg_P_Potential_Up_History);
                        tp = linspace(8, 32, len)'; 
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