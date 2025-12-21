function [SDCI, Rho] = calculate_metrics_post(P_AC, P_EV)
    % 计算 SDCI 和 Spearman Rho 指标
    % 对应 3.4 节的定义
    
    % SDCI (同向互补性指数)
    numerator = sum(min(P_AC, P_EV));
    denominator = sum(max(P_AC, P_EV));
    if denominator == 0
        SDCI = 0;
    else
        SDCI = numerator / denominator;
    end
    
    % Spearman Rho (秩相关系数)
    % 使用 MATLAB 内置 corr
    if length(P_AC) > 2
        Rho = corr(P_AC, P_EV, 'Type', 'Spearman');
    else
        Rho = 0;
    end
end