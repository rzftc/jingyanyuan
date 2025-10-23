% calculateSetpointProbabilitiesByArea.m
function area_probabilities = calculateSetpointProbabilitiesByArea(area_counts)
% calculateSetpointProbabilitiesByArea 计算空调设定温度的分布。
%
% 输入:
%   area_counts: 从 aggregateSetpointCountsByArea 函数得到的输出结构体。
%                每个字段 (功能区域) 包含 TsetMax/TsetMin 的值、计数和总数。
%
% 输出:
%   area_probabilities: 一个与 area_counts 结构相似的结构体，但包含概率而非计数。
%                       每个字段 (功能区域) 包含:
%                       - TsetMax_values: 唯一的 TsetMax 设定点向量。
%                       - P_TsetMax:      对应每个值的概率 P(TsetMax = X) 向量。
%                       - TsetMin_values: 唯一的 TsetMin 设定点向量。
%                       - P_TsetMin:      对应每个值的概率 P(TsetMin = X) 向量。

    area_probabilities = struct();
    functional_areas = fieldnames(area_counts);

    for i = 1:length(functional_areas)
        area_name = functional_areas{i};
        current_counts = area_counts.(area_name);
        
        area_dist = struct();

        % 计算 P_max
        if current_counts.total_TsetMax > 0
            area_dist.TsetMax_values = current_counts.TsetMax_values;
            area_dist.P_TsetMax = current_counts.TsetMax_counts(:) / current_counts.total_TsetMax; % 确保是列向量
        else
            area_dist.TsetMax_values = [];
            area_dist.P_TsetMax = [];
        end

        % 计算 P_min
        if current_counts.total_TsetMin > 0
            area_dist.TsetMin_values = current_counts.TsetMin_values;
            area_dist.P_TsetMin = current_counts.TsetMin_counts(:) / current_counts.total_TsetMin; % 确保是列向量
        else
            area_dist.TsetMin_values = [];
            area_dist.P_TsetMin = [];
        end
        
        area_probabilities.(area_name) = area_dist;
    end
    disp('已按区域计算空调设定点概率。');
end