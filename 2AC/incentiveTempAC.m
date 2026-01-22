function [T_set_up, T_set_down, deltaT] = incentiveTempAC(p, p_min, p_max, p_min_prime, p_max_prime, T_set_max)
    % incentiveTempAC: 根据激励电价计算温度变化的上下限
    % 并返回 deltaT（温度变化值）
    %
    % 输入:
    %   p - 激励电价，可以是标量或向量
    %   p_min - 最低电价 p_min
    %   p_max - 最高电价 p_max
    %   p_min_prime - 最低电价修正值 p_min'
    %   p_max_prime - 最高电价修正值 p_max'
    %   T_set_max - 最大温度设定变化值
    %
    % 输出:
    %   T_set_up - 温度变化的上限
    %   T_set_down - 温度变化的下限
    %   deltaT - 温度变化值

    % 初始化输出
    T_set_up = zeros(size(p));
    T_set_down = zeros(size(p));
    deltaT = zeros(size(p));

    % 计算 T_set_up（公式 22）
    for i = 1:length(p)
        if p(i) >= 0 && p(i) <= p_max_prime
            T_set_up(i) = (T_set_max / (p_max_prime - p_min_prime)) * (p(i) - p_min_prime);
        elseif p(i) > p_max_prime
            T_set_up(i) = T_set_max;
        else
            T_set_up(i) = 0; % 对于无效电价，设为 0
        end
    end

    % 计算 T_set_down（公式 23）
    for i = 1:length(p)
        if p(i) >= 0 && p(i) <= p_max
            T_set_down(i) = (T_set_max / (p_max - p_min)) * (p(i) - p_min);
        elseif p(i) > p_max
            T_set_down(i) = T_set_max;
        else
            T_set_down(i) = 0; % 对于无效电价，设为 0
        end
    end

    % 计算 deltaT（温度变化值，上下限之间的随机值）
    for i = 1:length(p)
        deltaT(i) = T_set_down(i) + (T_set_up(i) - T_set_down(i)) * rand();
    end
end
