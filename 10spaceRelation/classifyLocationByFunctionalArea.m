% classifyLocationByFunctionalArea.m
function functional_area = classifyLocationByFunctionalArea(latitude, longitude, zone_definitions)
% classifyLocationByFunctionalArea 根据经纬度和功能区定义，对地点进行分类。
%
% 输入:
%   latitude:         单个地点的纬度 (数值)。
%   longitude:        单个地点的经度 (数值)。
%   zone_definitions: 一个元胞数组 (cell array)，每个元胞是一个结构体，定义一个功能区。
%                     每个功能区结构体应包含:
%                       - name: 功能区的名称 (字符串, 例如 'Residential', 'Commercial')。
%                       - polygon_x: 定义该区域边界的多边形顶点的 X 坐标 (或经度) 向量。
%                                    (注意：向量的第一个和最后一个点应相同以闭合多边形)
%                       - polygon_y: 定义该区域边界的多边形顶点的 Y 坐标 (或纬度) 向量。
%                                    (注意：向量的第一个和最后一个点应相同以闭合多边形)
%
% 输出:
%   functional_area:  地点所属功能区的名称 (字符串)。
%                     如果没有落在任何定义的功能区内，则返回 'Unknown'。

    functional_area = 'Unknown'; % 默认值

    if ~iscell(zone_definitions)
        error('zone_definitions 必须是一个元胞数组。');
    end

    for i = 1:length(zone_definitions)
        zone = zone_definitions{i};
        if ~isstruct(zone) || ~all(isfield(zone, {'name', 'polygon_x', 'polygon_y'}))
            warning('zone_definitions 中的第 %d 个元素格式不正确，已跳过。', i);
            continue;
        end

        % 确保多边形坐标是数值向量
        if ~isnumeric(zone.polygon_x) || ~isvector(zone.polygon_x) || ...
           ~isnumeric(zone.polygon_y) || ~isvector(zone.polygon_y)
            warning('功能区 "%s" 的多边形坐标格式不正确，已跳过。', zone.name);
            continue;
        end
        
        % 确保多边形闭合 (首尾点相同) - 可选检查，inpolygon通常能处理
        % if zone.polygon_x(1) ~= zone.polygon_x(end) || zone.polygon_y(1) ~= zone.polygon_y(end)
        %     warning('功能区 "%s" 的多边形未闭合，结果可能不准确。', zone.name);
        % end

        % 使用 inpolygon 函数判断点是否在多边形内
        % 注意: inpolygon 需要 Xq, Yq, Xv, Yv 的顺序
        % 这里假设 longitude 对应 X, latitude 对应 Y
        try
            is_inside = inpolygon(longitude, latitude, zone.polygon_x, zone.polygon_y);
        catch ME_inpolygon
            warning('功能区 "%s" 调用inpolygon时出错: %s，已跳过。', zone.name, ME_inpolygon.message);
            continue;
        end

        if is_inside
            functional_area = zone.name;
            return; % 一旦找到所属区域，即可返回
        end
    end
end