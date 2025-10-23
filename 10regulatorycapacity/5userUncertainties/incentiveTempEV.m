function [E_tar_up, E_tar_down, deltaE] = incentiveTempEV(p, p_min, p_max, p_min_prime, p_max_prime, E_tar_max)
    % 输入参数说明:
    % p - 激励电价向量 [Nx1]
    % E_tar_max - 每个EV的最大电量变化向量 [Nx1]
    
    N = length(p);
    E_tar_up = zeros(N,1);
    E_tar_down = zeros(N,1);
    deltaE = zeros(N,1);
    
    for i = 1:N
        % 计算 E_tar_up
        if p(i) >= p_min_prime && p(i) <= p_max_prime
            E_tar_up(i) = (E_tar_max(i) / (p_max_prime - p_min_prime)) * (p(i) - p_min_prime);
        elseif p(i) > p_max_prime
            E_tar_up(i) = E_tar_max(i);
        else
            E_tar_up(i) = 0;
        end
        
        % 计算 E_tar_down
        if p(i) >= p_min && p(i) <= p_max
            E_tar_down(i) = (E_tar_max(i) / (p_max - p_min)) * (p(i) - p_min);
        elseif p(i) > p_max
            E_tar_down(i) = E_tar_max(i);
        else
            E_tar_down(i) = 0;
        end
        
        % 生成随机变化量
        if E_tar_up(i)>E_tar_down(i)
            deltaE(i) = E_tar_down(i) + (E_tar_up(i) - E_tar_down(i)) * rand();
        elseif E_tar_up(i)==E_tar_down(i)
            deltaE(i) = E_tar_up(i);
        end
    end
end



% function [E_tar_up, E_tar_down, deltaE] = incentiveTempEV(p, p_min, p_max, p_min_prime, p_max_prime, E_tar_max)
%     % incentiveTempEV: 根据激励电价计算目标能量的上下限变化
%     % 并返回 deltaE（目标能量变化值）
%     %
%     % 输入:
%     %   p - 激励电价，可以是标量或向量
%     %   p_min - 最低电价 p_min
%     %   p_max - 最高电价 p_max
%     %   p_min_prime - 最低电价修正值 p_min'
%     %   p_max_prime - 最高电价修正值 p_max'
%     %   E_tar_max - 最大目标能量变化值
%     %
%     % 输出:
%     %   E_tar_up - 目标能量的上限变化
%     %   E_tar_down - 目标能量的下限变化
%     %   deltaE - 实际目标能量变化值（上下限之间的随机值）
% 
%     % 初始化输出
%     E_tar_up = zeros(size(p));
%     E_tar_down = zeros(size(p));
%     deltaE = zeros(size(p));
% 
%     % 计算 E_tar_up（公式 22）
%     for i = 1:length(p)
%         if p(i) >= 0 && p(i) <= p_max_prime
%             E_tar_up(i) = (E_tar_max / (p_max_prime - p_min_prime)) * (p(i) - p_min_prime);
%         elseif p(i) > p_max_prime
%             E_tar_up(i) = E_tar_max;
%         else
%             E_tar_up(i) = 0; % 对于无效电价，设为 0
%         end
%     end
% 
%     % 计算 E_tar_down（公式 23）
%     for i = 1:length(p)
%         if p(i) >= 0 && p(i) <= p_max
%             E_tar_down(i) = (E_tar_max / (p_max - p_min)) * (p(i) - p_min);
%         elseif p(i) > p_max
%             E_tar_down(i) = E_tar_max;
%         else
%             E_tar_down(i) = 0; % 对于无效电价，设为 0
%         end
%     end
% 
%     % 计算 deltaE（目标能量的随机变化值，上下限之间随机取值）
%     for i = 1:length(p)
%         deltaE(i) = E_tar_down(i) + (E_tar_up(i) - E_tar_down(i)) * rand();
%     end
% end
