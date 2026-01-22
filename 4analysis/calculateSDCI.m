function [SDCI] = calculateSDCI(n_AC, n_EV, deltaP_AC, deltaP_EV)
% 计算同方向互补性指数 (SDCI)
% 输入参数：
%   n_AC   - 空调参与数量向量 [T×1]
%   n_EV   - EV参与数量向量 [T×1]
%   deltaP_AC - 单台空调调节能力向量 [T×1] (正值)
%   deltaP_EV - 单台EV调节能力向量 [T×1] (正值)
% 输出：
%   SDCI - 调节互补性指数

    % 计算总调节能力
    AC = n_AC .* deltaP_AC;
    EV = n_EV .* deltaP_EV;
    
    % 计算逐时段min/max
    min_vals = min(AC, EV);
    max_vals = max(AC, EV);
    
    % 处理分母为零的情况
    total_max = sum(max_vals);
    if total_max == 0
        SDCI = 0;
    else
        SDCI = sum(min_vals) / total_max;
    end
end