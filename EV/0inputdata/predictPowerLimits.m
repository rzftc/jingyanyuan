function [P_min, P_max] = predictPowerLimits(EVs, current_time, H)
    % 示例实现：确保P_min <= P_max
    P_max = zeros(H,1);
    P_min = zeros(H,1);
    
    for t = 1:H
        % 计算当前可用EV数量
        available_EVs = sum([EVs.t_dep] > (current_time + t));
        
        % 功率上下限计算
        P_max(t) = available_EVs * max([EVs.P_N]); % 最大可能功率
        P_min(t) = 0.2 * P_max(t); % 最小功率设为最大值的20%
        
        % 确保P_min <= P_max
        if P_min(t) > P_max(t)
            P_min(t) = P_max(t);
        end
    end
end
