% loadSpatialAnalysisInputs.m
function [ev_data, ac_data, zone_definitions] = loadSpatialAnalysisInputs(ev_filename, ac_filename, zones_filename)
% loadSpatialAnalysisInputs 从Excel文件加载用于空间分析的输入数据。
%
% 输入:
%   ev_filename:    EV数据Excel文件名。
%   ac_filename:    AC数据Excel文件名。
%   zones_filename: 功能区定义Excel文件名。
%
% 输出:
%   ev_data:        包含EV数据的MATLAB table。
%   ac_data:        包含AC数据的MATLAB table。
%   zone_definitions: 包含功能区定义的MATLAB cell array。

    disp('开始从Excel文件加载输入数据...');
    ev_data = [];
    ac_data = [];
    zone_definitions = {};

    % --- 1. 加载 EV 数据 ---
    if exist(ev_filename, 'file')
        try
            ev_data = readtable(ev_filename, 'Sheet', 'EVData');
            % 确保 ChargingMode 是 categorical 类型，如果它是从 cellstr 读取的话
            if ismember('ChargingMode', ev_data.Properties.VariableNames) && iscellstr(ev_data.ChargingMode)
                ev_data.ChargingMode = categorical(ev_data.ChargingMode);
            end
            disp([ev_filename, ' 数据加载成功。']);
        catch ME
            warning(['加载 ', ev_filename, ' 失败: ', ME.message]);
        end
    else
        warning([ev_filename, ' 未找到。']);
    end

    % --- 2. 加载 AC 数据 ---
    if exist(ac_filename, 'file')
        try
            ac_data = readtable(ac_filename, 'Sheet', 'ACData');
            disp([ac_filename, ' 数据加载成功。']);
        catch ME
            warning(['加载 ', ac_filename, ' 失败: ', ME.message]);
        end
    else
        warning([ac_filename, ' 未找到。']);
    end

    % --- 3. 加载并处理功能区定义 ---
    if exist(zones_filename, 'file')
        try
            zones_table = readtable(zones_filename, 'Sheet', 'ZoneDefinitions');
            zone_definitions = cell(height(zones_table), 1);
            for i = 1:height(zones_table)
                zone = struct();
                zone.name = zones_table.ZoneName{i};
                
                % 解析逗号分隔的坐标字符串
                x_coords_str = strsplit(zones_table.PolygonX_CSV{i}, ',');
                y_coords_str = strsplit(zones_table.PolygonY_CSV{i}, ',');
                
                zone.polygon_x = str2double(x_coords_str);
                zone.polygon_y = str2double(y_coords_str);
                
                % 检查解析是否成功 (例如，是否有非数字导致NaN)
                if any(isnan(zone.polygon_x)) || any(isnan(zone.polygon_y))
                    warning('功能区 "%s" 的坐标解析包含NaN值，请检查Excel文件中的格式。', zone.name);
                end
                
                zone_definitions{i} = zone;
            end
            disp([zones_filename, ' 功能区定义加载并处理成功。']);
        catch ME
            warning(['加载或处理 ', zones_filename, ' 失败: ', ME.message]);
        end
    else
        warning([zones_filename, ' 未找到。']);
    end
    
    if isempty(ev_data)
        warning('EV数据为空，请检查对应的Excel文件和路径。');
    end
    if isempty(ac_data)
        warning('AC数据为空，请检查对应的Excel文件和路径。');
    end
    if isempty(zone_definitions)
        warning('功能区定义为空，地点分类将无法进行或不准确。请检查对应的Excel文件和路径。');
    end

    disp('输入数据加载过程结束。');
end