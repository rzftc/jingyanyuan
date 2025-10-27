function [deltaE_up, deltaE_down] = incentiveTempEV_updown(p, p_min, p_max, p_min_prime, p_max_prime, E_tar_max)
% incentiveTempEV: 根据激励电价计算能量灵活性窗口的上下边界。
%
% 输入:
%   p           - 激励电价，可以是标量或向量
%   p_min       - 原始激励下限 (用于下调灵活性)
%   p_max       - 原始激励上限 (用于下调灵活性)
%   p_min_prime - 调整后激励下限 (用于上调灵活性)
%   p_max_prime - 调整后激励上限 (用于上调灵活性)
%   E_tar_max   - 该车辆愿意提供的最大电量变化（灵活性）
%
% 输出:
%   deltaE_up   - 上调灵活性：允许的能量目标增加量 (p > p_min_prime 时 > 0)
%   deltaE_down - 下调灵活性：允许的能量目标减少量 (p > p_min 时 > 0)

    % 确保p是列向量以便于处理
    p = p(:);
    N = length(p);
    deltaE_up = zeros(N,1);
    deltaE_down = zeros(N,1);

    % 为每个激励价格计算灵活性
    for i = 1:N
        % --- 计算上调灵活性 (deltaE_up) ---
        % 基于公式(22)，代表目标电量可增加的量
        if p(i) >= p_min_prime && p(i) <= p_max_prime
            deltaE_up(i) = (E_tar_max(i) / (p_max_prime - p_min_prime)) * (p(i) - p_min_prime);
        elseif p(i) > p_max_prime
            deltaE_up(i) = E_tar_max(i);
        else
            deltaE_up(i) = 0;
        end
        
        % --- 计算下调灵活性 (deltaE_down) ---
        % 基于公式(23)，代表目标电量可减少的量
        if p(i) >= p_min && p(i) <= p_max
            deltaE_down(i) = (E_tar_max(i) / (p_max - p_min)) * (p(i) - p_min);
        elseif p(i) > p_max
            deltaE_down(i) = E_tar_max(i);
        else
            deltaE_down(i) = 0;
        end
    end
end