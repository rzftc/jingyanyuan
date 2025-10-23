% analyzeACSpatialCharacteristics.m
function ac_spatial_stats = analyzeACSpatialCharacteristics(ac_data, zone_definitions)
% analyzeACSpatialCharacteristics 分析空调的空间特性，包括设定温度分布。
%
% 输入:
%   ac_data: 一个表 (table) 或结构体数组，包含AC数据。
%            必须包含 'Latitude' 和 'Longitude' 列/字段。
%            必须包含 'TsetMax' 和 'TsetMin' 用于温度分布分析。
%   zone_definitions: 用于 classifyLocationByFunctionalArea 的功能区定义。
%
% 输出:
%   ac_spatial_stats: 一个结构体，其字段名是各功能区的名称 (以及一个 'Overall' 字段)。
%                     每个功能区字段下包含该区域的AC统计信息。

    if ~(istable(ac_data) || isstruct(ac_data))
        error('ac_data 必须是一个表或结构体数组。');
    end

    required_fields = {'Latitude', 'Longitude', 'TsetMax', 'TsetMin'};
    if istable(ac_data)
        if ~all(ismember(required_fields, ac_data.Properties.VariableNames))
            error('ac_data 表中必须包含 "Latitude", "Longitude", "TsetMax", "TsetMin" 列。');
        end
        num_acs = height(ac_data);
        ac_data_table = ac_data; 
    else 
        if num_acs > 0 && ~all(isfield(ac_data, required_fields)) % isstruct
            error('ac_data 结构体数组中必须包含 "Latitude", "Longitude", "TsetMax", "TsetMin" 字段。');
        end
        num_acs = length(ac_data);
        if num_acs > 0
            ac_data_table = struct2table(ac_data);
        else
            % 创建一个包含所需列的空表，以避免后续函数出错
            varTypes = {'double', 'double', 'double', 'double'};
            ac_data_table = table('Size', [0, length(required_fields)], ...
                                  'VariableTypes', varTypes, ...
                                  'VariableNames', required_fields);
        end
    end

    % 为每个AC确定其功能区并添加到表中
    classified_area_column_data = cell(num_acs, 1);
    if num_acs > 0
        for i = 1:num_acs
            lat = ac_data_table.Latitude(i);
            lon = ac_data_table.Longitude(i);
             if isnumeric(lat) && isnumeric(lon) && ~isnan(lat) && ~isnan(lon)
                classified_area_column_data{i} = classifyLocationByFunctionalArea(lat, lon, zone_definitions);
            else
                classified_area_column_data{i} = 'InvalidCoordinates';
            end
        end
        ac_data_table.FunctionalArea = categorical(classified_area_column_data);
    else
        % 如果没有AC数据，创建一个带FunctionalArea列的空表
         ac_data_table.FunctionalArea = categorical(cell(0,1));
    end

    % 1. 聚合设定点计数
    try
        ac_setpoint_counts = aggregateSetpointCountsByArea(ac_data_table);
    catch ME
        error('调用 aggregateSetpointCountsByArea 失败: %s。请确保该函数在路径中且输入正确。', ME.message);
    end

    % 2. 计算设定点概率
    try
        ac_setpoint_probabilities = calculateSetpointProbabilitiesByArea(ac_setpoint_counts);
    catch ME
        error('调用 calculateSetpointProbabilitiesByArea 失败: %s。请确保该函数在路径中且输入正确。', ME.message);
    end
    
    ac_spatial_stats = struct();
    
    % --- (可选) 计算总体统计数据 ---
    ac_spatial_stats.Overall.total_ac_units = num_acs;
    % 可以创建一个不包含FunctionalArea列的临时表来计算总体温度分布
    if num_acs > 0
        temp_ac_data_for_overall = ac_data_table;
        temp_ac_data_for_overall.FunctionalArea = repmat(categorical({'OverallAnalysis'}), num_acs, 1);
        overall_counts = aggregateSetpointCountsByArea(temp_ac_data_for_overall);
        overall_probs = calculateSetpointProbabilitiesByArea(overall_counts);
        if isfield(overall_probs, 'OverallAnalysis')
             ac_spatial_stats.Overall.setpoint_temperature_distribution = overall_probs.OverallAnalysis;
        else
             ac_spatial_stats.Overall.setpoint_temperature_distribution = struct('TsetMax_values', [], 'P_TsetMax', [], 'TsetMin_values', [], 'P_TsetMin', []);
        end
    else
        ac_spatial_stats.Overall.setpoint_temperature_distribution = struct('TsetMax_values', [], 'P_TsetMax', [], 'TsetMin_values', [], 'P_TsetMin', []);
    end


    % 将按区域的温度分布存储到输出结构体
    functional_area_names_from_probs = fieldnames(ac_setpoint_probabilities);
    for i = 1:length(functional_area_names_from_probs)
        area_name = functional_area_names_from_probs{i};
        % 获取该区域的AC数量
        if num_acs > 0
            num_in_area = sum(ac_data_table.FunctionalArea == categorical(cellstr(area_name)));
        else
            num_in_area = 0;
        end
        
        ac_spatial_stats.(area_name).total_ac_units = num_in_area;
        ac_spatial_stats.(area_name).setpoint_temperature_distribution = ac_setpoint_probabilities.(area_name);
        
        % 可在此处添加其他AC空间特性分析
        if num_in_area > 0 && istable(ac_data_table) && ismember('RatedPower', ac_data_table.Properties.VariableNames)
           area_ac_data_subset = ac_data_table(ac_data_table.FunctionalArea == categorical(cellstr(area_name)), :);
           ac_spatial_stats.(area_name).average_rated_power = mean(area_ac_data_subset.RatedPower, 'omitnan');
        end
    end
    
    % 确保所有定义的功能区（以及Unknown, InvalidCoordinates）都在输出结构体中，即使没有数据
    all_possible_area_names = cellfun(@(x) x.name, zone_definitions, 'UniformOutput', false);
    all_possible_area_names = [all_possible_area_names; 'Unknown'; 'InvalidCoordinates'];
    all_possible_area_names = unique(all_possible_area_names); % 去重

    for i = 1:length(all_possible_area_names)
        area_name_check = matlab.lang.makeValidName(all_possible_area_names{i});
        if ~isfield(ac_spatial_stats, area_name_check)
            ac_spatial_stats.(area_name_check).total_ac_units = 0;
            ac_spatial_stats.(area_name_check).setpoint_temperature_distribution = struct('TsetMax_values', [], 'P_TsetMax', [], 'TsetMin_values', [], 'P_TsetMin', []);
        end
    end

    disp('已完成空调空间特性分析。');
end