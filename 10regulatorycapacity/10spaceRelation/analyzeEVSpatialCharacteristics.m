% analyzeEVSpatialCharacteristics.m
function ev_spatial_stats = analyzeEVSpatialCharacteristics(ev_data, zone_definitions)
% analyzeEVSpatialCharacteristics 分析电动汽车充电的空间特性。
%
% 输入:
%   ev_data: 一个表 (table) 或结构体数组，包含EV充电数据。
%            必须包含 'Latitude' 和 'Longitude' 列/字段。
%            建议包含以下指标作为列/字段 (根据您的数据调整):
%            - 'ChargingAmount' (单次充电电量 kWh)
%            - 'InitialSOC' (充电前剩余电量，0-1范围)
%            - 'ChargingMode' (充电方式, 例如 'Fast', 'Slow', 或具体功率值)
%            - 'ChargingDuration' (充电时长, 小时或分钟)
%            - 'StartTime' (起始充电时刻, datetime 或数值小时)
%   zone_definitions: 用于 classifyLocationByFunctionalArea 的功能区定义。
%
% 输出:
%   ev_spatial_stats: 一个结构体，其字段名是各功能区的名称 (以及一个 'Overall' 字段表示总体统计)。
%                     每个功能区字段下包含该区域的EV充电统计信息。

    if ~(istable(ev_data) || isstruct(ev_data))
        error('ev_data 必须是一个表或结构体数组。');
    end
    
    if istable(ev_data)
        if ~all(ismember({'Latitude', 'Longitude'}, ev_data.Properties.VariableNames))
            error('ev_data 表中必须包含 "Latitude" 和 "Longitude" 列。');
        end
        num_evs = height(ev_data);
        % 先将表转换为结构体数组，便于统一处理字段缺失的情况
        if num_evs > 0
            ev_data_struct_array = table2struct(ev_data);
        else
            ev_data_struct_array = struct([]); % 空结构体数组
        end
    else % isstruct
        if num_evs > 0 && ~all(isfield(ev_data, {'Latitude', 'Longitude'}))
             error('ev_data 结构体数组中必须包含 "Latitude" 和 "Longitude" 字段。');
        end
        num_evs = length(ev_data);
        ev_data_struct_array = ev_data;
    end

    % 为每个EV确定其功能区
    classified_areas = cell(num_evs, 1);
    for i = 1:num_evs
        lat = ev_data_struct_array(i).Latitude;
        lon = ev_data_struct_array(i).Longitude;
        if isnumeric(lat) && isnumeric(lon) && ~isnan(lat) && ~isnan(lon)
            classified_areas{i} = classifyLocationByFunctionalArea(lat, lon, zone_definitions);
        else
            classified_areas{i} = 'InvalidCoordinates'; 
        end
    end
    
    unique_functional_areas = unique(classified_areas);
    ev_spatial_stats = struct();

    % --- (可选) 计算总体统计数据 ---
    all_stats = struct();
    all_stats.total_charging_sessions = num_evs;
    
    if num_evs > 0 % 仅当有数据时计算统计量
        if isfield(ev_data_struct_array, 'ChargingAmount') && isnumeric([ev_data_struct_array.ChargingAmount])
            all_stats.total_charging_amount = sum([ev_data_struct_array.ChargingAmount], 'omitnan');
            all_stats.average_charging_amount = mean([ev_data_struct_array.ChargingAmount], 'omitnan');
        end
        if isfield(ev_data_struct_array, 'InitialSOC') && isnumeric([ev_data_struct_array.InitialSOC])
            all_stats.average_initial_SOC = mean([ev_data_struct_array.InitialSOC], 'omitnan');
        end
        if isfield(ev_data_struct_array, 'ChargingMode')
            modes_data = {ev_data_struct_array.ChargingMode};
            if ~all(cellfun(@ischar, modes_data)) && ~iscategorical([modes_data{:}]) && ~isstring([modes_data{:}])
                 warning('ChargingMode 包含非字符/分类/字符串数据，分布统计可能不准确。');
                 all_stats.charging_mode_distribution = table();
            else
                [mode_counts, mode_names] = groupcounts(categorical(modes_data'));
                all_stats.charging_mode_distribution = table(mode_names, mode_counts, mode_counts/sum(mode_counts)*100, 'VariableNames', {'模式', '次数', '百分比'});
            end
        end
        if isfield(ev_data_struct_array, 'ChargingDuration') && isnumeric([ev_data_struct_array.ChargingDuration])
            all_stats.average_charging_duration = mean([ev_data_struct_array.ChargingDuration], 'omitnan');
        end
        if isfield(ev_data_struct_array, 'StartTime') % 假设StartTime是数值型小时 (0-23.99)
             if isnumeric([ev_data_struct_array.StartTime])
                start_hours_all = [ev_data_struct_array.StartTime];
                all_stats.start_time_hourly_distribution_counts = histcounts(start_hours_all, 0:24);
                all_stats.start_time_hourly_distribution_prob = histcounts(start_hours_all, 0:24, 'Normalization', 'probability');
             elseif isdatetime([ev_data_struct_array.StartTime])
                start_hours_all = hour([ev_data_struct_array.StartTime]);
                all_stats.start_time_hourly_distribution_counts = histcounts(start_hours_all, -0.5:23.5); % Datetime hour is 0-23
                all_stats.start_time_hourly_distribution_prob = histcounts(start_hours_all, -0.5:23.5, 'Normalization', 'probability');
             end
        end
    end
    ev_spatial_stats.Overall = all_stats;
    % --- 结束总体统计 ---

    % 按功能区进行统计
    for k = 1:length(unique_functional_areas)
        area_name_str = unique_functional_areas{k};
        valid_area_field_name = matlab.lang.makeValidName(area_name_str);
        
        indices_in_area = strcmp(classified_areas, area_name_str);
        area_ev_data_struct = ev_data_struct_array(indices_in_area);
        num_in_area = sum(indices_in_area);

        if num_in_area == 0
            stats = struct('total_charging_sessions', 0); % 即使区域为空也创建一个条目
            ev_spatial_stats.(valid_area_field_name) = stats;
            continue; 
        end

        stats = struct();
        stats.total_charging_sessions = num_in_area;

        if isfield(area_ev_data_struct, 'ChargingAmount') && isnumeric([area_ev_data_struct.ChargingAmount])
            stats.total_charging_amount = sum([area_ev_data_struct.ChargingAmount], 'omitnan');
            stats.average_charging_amount = mean([area_ev_data_struct.ChargingAmount], 'omitnan');
        end
        if isfield(area_ev_data_struct, 'InitialSOC') && isnumeric([area_ev_data_struct.InitialSOC])
            stats.average_initial_SOC = mean([area_ev_data_struct.InitialSOC], 'omitnan');
        end
        if isfield(area_ev_data_struct, 'ChargingMode')
            area_modes_data = {area_ev_data_struct.ChargingMode};
             if ~all(cellfun(@ischar, area_modes_data)) && ~iscategorical([area_modes_data{:}]) && ~isstring([area_modes_data{:}])
                 stats.charging_mode_distribution = table();
            else
                [mode_counts, mode_names] = groupcounts(categorical(area_modes_data'));
                stats.charging_mode_distribution = table(mode_names, mode_counts, mode_counts/sum(mode_counts)*100, 'VariableNames', {'模式', '次数', '百分比'});
            end
        end
        if isfield(area_ev_data_struct, 'ChargingDuration') && isnumeric([area_ev_data_struct.ChargingDuration])
            stats.average_charging_duration = mean([area_ev_data_struct.ChargingDuration], 'omitnan');
        end
        if isfield(area_ev_data_struct, 'StartTime')
            if isnumeric([area_ev_data_struct.StartTime])
                area_start_hours = [area_ev_data_struct.StartTime];
                stats.start_time_hourly_distribution_counts = histcounts(area_start_hours, 0:24);
                stats.start_time_hourly_distribution_prob = histcounts(area_start_hours, 0:24, 'Normalization', 'probability');
            elseif isdatetime([area_ev_data_struct.StartTime])
                area_start_hours = hour([area_ev_data_struct.StartTime]);
                stats.start_time_hourly_distribution_counts = histcounts(area_start_hours, -0.5:23.5);
                stats.start_time_hourly_distribution_prob = histcounts(area_start_hours, -0.5:23.5, 'Normalization', 'probability');
            end
        end
        ev_spatial_stats.(valid_area_field_name) = stats;
    end
    disp('已完成电动汽车空间特性分析。');
end