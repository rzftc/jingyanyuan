% generateSpatialAnalysisInputFiles.m
function generateSpatialAnalysisInputFiles(ev_filename, ac_filename, zones_filename, overwrite_existing)
% generateSpatialAnalysisInputFiles 生成用于空间分析的示例Excel输入文件。
%
% 输入:
%   ev_filename:        EV数据Excel文件名 (例如 'ev_data_example.xlsx')。
%   ac_filename:        AC数据Excel文件名 (例如 'ac_data_example.xlsx')。
%   zones_filename:     功能区定义Excel文件名 (例如 'zone_definitions_example.xlsx')。
%   overwrite_existing: 布尔值，如果为true，则覆盖已存在的文件；否则，如果文件已存在则跳过生成。
%                       默认为 false。

    if nargin < 4
        overwrite_existing = false;
    end

    disp('开始生成示例Excel输入文件...');

    % --- 1. 生成 EV 数据 Excel 文件 ---
    if ~exist(ev_filename, 'file') || overwrite_existing
        disp(['正在生成: ', ev_filename]);
        num_ev_samples = 100;
        Latitude = 39.88 + (40.02 - 39.88) * rand(num_ev_samples, 1);
        Longitude = 116.28 + (116.42 - 116.28) * rand(num_ev_samples, 1);
        ChargingAmount = round(10 + 40 * rand(num_ev_samples, 1), 2); % kWh
        InitialSOC = round(0.1 + 0.7 * rand(num_ev_samples, 1), 2);    % 10% - 80%
        ChargingModeOptions = {'Fast'; 'Slow'};
        ChargingMode = ChargingModeOptions(randi(length(ChargingModeOptions), num_ev_samples, 1));
        ChargingDuration = round(0.5 + 7.5 * rand(num_ev_samples, 1), 1); % hours
        StartTime = round(mod(randn(num_ev_samples,1)*4 + 18, 24), 1); % 数值型小时 (0-23.9)

        ev_table = table(Latitude, Longitude, ChargingAmount, InitialSOC, ...
                         cellstr(ChargingMode), ChargingDuration, StartTime, ... % 存储为cellstr以兼容性
                         'VariableNames', {'Latitude', 'Longitude', 'ChargingAmount', 'InitialSOC', ...
                                           'ChargingMode', 'ChargingDuration', 'StartTime'});
        try
            writetable(ev_table, ev_filename, 'Sheet', 'EVData');
            disp([ev_filename, ' 已成功生成。']);
        catch ME
            warning(['生成 ', ev_filename, ' 失败: ', ME.message]);
        end
    else
        disp([ev_filename, ' 已存在，跳过生成。']);
    end

    % --- 2. 生成 AC 数据 Excel 文件 ---
    if ~exist(ac_filename, 'file') || overwrite_existing
        disp(['正在生成: ', ac_filename]);
        num_ac_samples = 150;
        Latitude = 39.88 + (40.02 - 39.88) * rand(num_ac_samples, 1);
        Longitude = 116.28 + (116.42 - 116.28) * rand(num_ac_samples, 1);
        TsetMax = round(24 + 6 * rand(num_ac_samples, 1)); % 24-30 度
        TsetMin = round(18 + 5 * rand(num_ac_samples, 1)); % 18-23 度
        % 确保 Tmin <= Tmax
        idx_swap = TsetMin > TsetMax;
        temp_min = TsetMin(idx_swap);
        TsetMin(idx_swap) = TsetMax(idx_swap);
        TsetMax(idx_swap) = temp_min;
        % 确保不完全相同
        idx_equal = TsetMin == TsetMax;
        TsetMax(idx_equal & TsetMax < 30) = TsetMax(idx_equal & TsetMax < 30) + 1;
        TsetMin(idx_equal & TsetMin > 18) = TsetMin(idx_equal & TsetMin > 18) - 1;
        
        RatedPower = round(0.8 + (3.5 - 0.8) * rand(num_ac_samples, 1), 2); % kW

        ac_table = table(Latitude, Longitude, TsetMax, TsetMin, RatedPower, ...
                         'VariableNames', {'Latitude', 'Longitude', 'TsetMax', 'TsetMin', 'RatedPower'});
        try
            writetable(ac_table, ac_filename, 'Sheet', 'ACData');
            disp([ac_filename, ' 已成功生成。']);
        catch ME
            warning(['生成 ', ac_filename, ' 失败: ', ME.message]);
        end
    else
        disp([ac_filename, ' 已存在，跳过生成。']);
    end

    % --- 3. 生成功能区定义 Excel 文件 ---
    if ~exist(zones_filename, 'file') || overwrite_existing
        disp(['正在生成: ', zones_filename]);
        ZoneName = {'Residential'; 'Commercial'; 'Office'};
        % 坐标用逗号分隔的字符串形式存储
        PolygonX_CSV = {
            '116.30,116.35,116.35,116.30,116.30'; % Residential X
            '116.36,116.40,116.40,116.36,116.36'; % Commercial X
            '116.30,116.35,116.35,116.30,116.30'  % Office X
        };
        PolygonY_CSV = {
            '39.90,39.90,39.95,39.95,39.90';     % Residential Y
            '39.90,39.90,39.95,39.95,39.90';     % Commercial Y
            '39.96,39.96,40.00,40.00,39.96'      % Office Y
        };
        
        zones_table = table(ZoneName, PolygonX_CSV, PolygonY_CSV);
        try
            writetable(zones_table, zones_filename, 'Sheet', 'ZoneDefinitions');
            disp([zones_filename, ' 已成功生成。']);
        catch ME
            warning(['生成 ', zones_filename, ' 失败: ', ME.message]);
        end
    else
        disp([zones_filename, ' 已存在，跳过生成。']);
    end
    
    disp('示例Excel输入文件生成过程结束。');
end