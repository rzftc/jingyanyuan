% aggregateSetpointCountsByArea.m
function area_counts = aggregateSetpointCountsByArea(ac_data)
% aggregateSetpointCountsByArea 按功能区域聚合空调设定温度的计数。
%
% 输入:
%   ac_data: 一个表 (table)，至少包含以下列:
%            'FunctionalArea' - 分类或字符串类型，指定区域类型 (例如, 'Residential', 'Commercial')。
%                               此列由 classifyLocationByFunctionalArea 函数生成或预先存在。
%            'TsetMax' - 数值类型，最高温度设定点。
%            'TsetMin' - 数值类型，最低温度设定点。
%
% 输出:
%   area_counts: 一个结构体，其每个字段名对应 ac_data 中一个唯一的功能区域名称。
%                每个字段包含一个子结构体，内容如下:
%                - TsetMax_values: 该区域内观察到的唯一 TsetMax 设定点向量。
%                - TsetMax_counts: 对应 TsetMax_values 中每个值的计数向量。
%                - total_TsetMax:  该区域内 TsetMax 观测总数。
%                - TsetMin_values: 该区域内观察到的唯一 TsetMin 设定点向量。
%                - TsetMin_counts: 对应 TsetMin_values 中每个值的计数向量。
%                - total_TsetMin:  该区域内 TsetMin 观测总数。

    if ~istable(ac_data) || ~all(ismember({'FunctionalArea', 'TsetMax', 'TsetMin'}, ac_data.Properties.VariableNames))
        error('输入 ac_data 必须是包含 ''FunctionalArea'', ''TsetMax'', 和 ''TsetMin'' 列的表 (table)。');
    end

    % 确保FunctionalArea是categorical类型，以便unique正确处理
    if ~iscategorical(ac_data.FunctionalArea)
        unique_areas = unique(string(ac_data.FunctionalArea)); % 先转string再unique
    else
        unique_areas = unique(ac_data.FunctionalArea);
    end
    
    area_counts = struct();

    for i = 1:length(unique_areas)
        current_area_val = unique_areas(i);
        current_area_name_str = string(current_area_val); %转换为字符串以用于字段名
        % 确保字段名是有效的MATLAB标识符
        valid_area_field_name = matlab.lang.makeValidName(current_area_name_str);

        % 根据FunctionalArea的值筛选数据
        if iscategorical(ac_data.FunctionalArea)
            area_data = ac_data(ac_data.FunctionalArea == current_area_val, :);
        else % 如果原始是cellstr或char array转换来的string
            area_data = ac_data(string(ac_data.FunctionalArea) == current_area_name_str, :);
        end


        area_stats = struct();

        % 处理 TsetMax
        if ~isempty(area_data.TsetMax) && any(~isnan(area_data.TsetMax))
            valid_TsetMax = area_data.TsetMax(~isnan(area_data.TsetMax)); % 移除NaN值
            if ~isempty(valid_TsetMax)
                [G_max, unique_TsetMax] = findgroups(valid_TsetMax);
                counts_TsetMax = splitapply(@numel, valid_TsetMax, G_max);
                area_stats.TsetMax_values = unique_TsetMax;
                area_stats.TsetMax_counts = counts_TsetMax;
                area_stats.total_TsetMax = sum(counts_TsetMax);
            else
                area_stats.TsetMax_values = [];
                area_stats.TsetMax_counts = [];
                area_stats.total_TsetMax = 0;
            end
        else
            area_stats.TsetMax_values = [];
            area_stats.TsetMax_counts = [];
            area_stats.total_TsetMax = 0;
        end

        % 处理 TsetMin
        if ~isempty(area_data.TsetMin) && any(~isnan(area_data.TsetMin))
             valid_TsetMin = area_data.TsetMin(~isnan(area_data.TsetMin)); % 移除NaN值
            if ~isempty(valid_TsetMin)
                [G_min, unique_TsetMin] = findgroups(valid_TsetMin);
                counts_TsetMin = splitapply(@numel, valid_TsetMin, G_min);
                area_stats.TsetMin_values = unique_TsetMin;
                area_stats.TsetMin_counts = counts_TsetMin;
                area_stats.total_TsetMin = sum(counts_TsetMin);
            else
                area_stats.TsetMin_values = [];
                area_stats.TsetMin_counts = [];
                area_stats.total_TsetMin = 0;
            end
        else
            area_stats.TsetMin_values = [];
            area_stats.TsetMin_counts = [];
            area_stats.total_TsetMin = 0;
        end
        
        area_counts.(valid_area_field_name) = area_stats;
    end
    disp('已按区域聚合空调设定点计数。');
end