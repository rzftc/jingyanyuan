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

