function participation_level = calculateParticipation(incentive_price, base_price)
    % calculateParticipation: 计算用户的参与度水平
    %
    % 输入：
    %   incentive_price - 激励电价（标量或向量，非负）
    %   base_price      - 基准电价（标量，必须 > 0）
    %
    % 输出：
    %   participation_level - 计算出的参与度水平，范围 [0, 1]
    
    % 输入验证
    if base_price <= 0
        error('base_price 必须大于 0');
    end
    if any(incentive_price < 0)
        error('incentive_price 必须非负');
    end

    % 计算参与度 (基于式21)
    participation_level = - (1 / (4 * base_price.^2)) * incentive_price.^2 ...
                          + (1 / base_price) * incentive_price;

    % 约束参与度范围 [0,1]
    participation_level = max(min(participation_level, 1), 0);
end
