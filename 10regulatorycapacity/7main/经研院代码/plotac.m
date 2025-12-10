% plot_AC_Data_From_Excel.m
clear; close all; clc;

%% 1. 参数设置
excel_filename = 'AC_Plot_Data_Selected.xlsx'; % 指定要读取的Excel文件名
output_png_up = 'AC_Cluster_Up_Regulation_from_Excel.png'; % 上调图输出文件名
output_png_down = 'AC_Cluster_Down_Regulation_from_Excel.png'; % 下调图输出文件名

%% 2. 从Excel文件读取数据
fprintf('正在从 %s 读取数据...\n', excel_filename);
try
    dataTable = readtable(excel_filename);
catch ME
    error('无法读取Excel文件 "%s"。请确保文件存在于当前目录或MATLAB路径中。\n错误信息: %s', excel_filename, ME.message);
end

% 检查数据是否成功读取
if isempty(dataTable)
    error('从 "%s" 读取的数据为空。', excel_filename);
end

% 提取时间列
if ~ismember('Time_Hours', dataTable.Properties.VariableNames)
    error('Excel文件中缺少 "Time_Hours" 列。');
end
time_points = dataTable.Time_Hours;

% 动态识别数据列 (上调和下调)
variableNames = dataTable.Properties.VariableNames;
up_columns = {};
down_columns = {};
prices_found = []; % 用于存储从列名中提取的价格

for i = 1:length(variableNames)
    colName = variableNames{i};
    % 查找上调列 (以 "Up_Price_" 开头)
    if startsWith(colName, 'Up_Price_')
        up_columns{end+1} = colName;
        % 尝试从列名提取价格值
        price_str = erase(colName, 'Up_Price_');
        price_val = str2double(strrep(price_str, '_', '.')); % 替换可能存在的下划线
        if ~isnan(price_val)
            prices_found(end+1) = price_val;
        end
    % 查找下调列 (以 "Down_Price_" 开头)
    elseif startsWith(colName, 'Down_Price_')
        down_columns{end+1} = colName;
        % 无需重复提取价格
    end
end

% 对找到的价格进行排序，以便图例顺序一致
prices_found = unique(prices_found);
prices_found = sort(prices_found);

if isempty(up_columns) || isempty(down_columns) || isempty(prices_found)
    error('未能从Excel文件中找到符合命名规则 (Up_Price_X 或 Down_Price_X) 的数据列。');
end

num_plots = length(prices_found);
fprintf('从Excel文件中识别出 %d 个价格水平的数据。\n', num_plots);

%% 3. 绘制图形
colors = lines(num_plots); % 为每个价格曲线获取颜色

% --- 3.1 绘制上调潜力图 ---
figure('Name', 'AC Cluster Up-Regulation Capacity from Excel', 'Position', [100 100 1000 600]);
hold on; grid on;
legend_entries_up = {};

for k = 1:num_plots
    current_price = prices_found(k);
    % 构造对应的列名
    up_col_name = matlab.lang.makeValidName(sprintf('Up_Price_%.1f', current_price));

    % 检查列是否存在于表中
    if ismember(up_col_name, dataTable.Properties.VariableNames)
        % 绘制数据
        plot(time_points, dataTable.(up_col_name), 'LineWidth', 1.5, 'Color', colors(k,:), ...
             'DisplayName', sprintf('电价: %.1f 分/℃', current_price));
        legend_entries_up{end+1} = sprintf('电价: %.1f 分/℃', current_price);
    else
        warning('在Excel数据中未找到列 "%s"，跳过绘制。', up_col_name);
    end
end
hold off;

xlabel('时间 (小时)', 'FontSize', 16);
ylabel('集群上调潜力 (kW)', 'FontSize', 16);
if ~isempty(legend_entries_up)
    legend(legend_entries_up, 'Location', 'best', 'FontSize', 14);
end
set(gca, 'FontSize', 14);

% 保存图像
print('-dpng', '-r400', output_png_up);
fprintf('空调集群上调潜力图已保存为 %s\n', output_png_up);

% --- 3.2 绘制下调潜力图 ---
figure('Name', 'AC Cluster Down-Regulation Capacity from Excel', 'Position', [100 750 1000 600]);
hold on; grid on;
legend_entries_down = {};

for k = 1:num_plots
    current_price = prices_found(k);
    % 构造对应的列名
    down_col_name = matlab.lang.makeValidName(sprintf('Down_Price_%.1f', current_price));

    % 检查列是否存在于表中
    if ismember(down_col_name, dataTable.Properties.VariableNames)
        % 绘制数据
        plot(time_points, dataTable.(down_col_name), 'LineWidth', 1.5, 'Color', colors(k,:), ...
             'DisplayName', sprintf('电价: %.1f 分/℃', current_price));
        legend_entries_down{end+1} = sprintf('电价: %.1f 分/℃', current_price);
    else
        warning('在Excel数据中未找到列 "%s"，跳过绘制。', down_col_name);
    end
end
hold off;

xlabel('时间 (小时)', 'FontSize', 16);
ylabel('集群下调潜力 (kW)', 'FontSize', 16);
if ~isempty(legend_entries_down)
    legend(legend_entries_down, 'Location', 'best', 'FontSize', 14);
end
set(gca, 'FontSize', 14);

% 保存图像
print('-dpng', '-r400', output_png_down);
fprintf('空调集群下调潜力图已保存为 %s\n', output_png_down);

fprintf('所有绘图任务完成。\n');