function [lambda_star] = aggregateEVs(EVs, P_tar)
    % 分离处于 LockON 和 LockOFF 状态的电动汽车
    lockON_mask  = strcmp({EVs.state}, 'LockON');
    lockOFF_mask = strcmp({EVs.state}, 'LockOFF');
    P_lockON     = sum([EVs(lockON_mask).P_N]);
    
    % 处理正常状态（非 LockON/LockOFF）的电动汽车
    normal_EVs = EVs(~lockON_mask & ~lockOFF_mask);
    % normal_EVs = EVs(~lockOFF_mask);
    % normal_EVs = EVs;
    
    % 生成 S' 特征值列表并按升序排列
    S_prime_list = arrayfun(@(ev) getSPrime(ev), normal_EVs);
    P_N_list     = [normal_EVs.P_N];
    [S_prime_sorted, sort_idx] = sort(S_prime_list, 'descend');  % S' 降序排序
    P_N_sorted   = P_N_list(sort_idx);
    
    % 生成从左到右的累积功率（cumulative_P(i) = sum(P_N_sorted(1:i))）
    cumulative_P    = cumsum(P_N_sorted);
    total_normal_P  = sum(P_N_sorted);
    
    % 计算剩余功率需求
    remaining_P = P_tar - P_lockON;
    % remaining_P = P_tar;
    
    % ================= 核心逻辑判断 =================
    if remaining_P <= 0
        % 情况 1：仅需 LockON 功率即可满足需求
        lambda_star = 1;  % 不激活任何正常状态 EV
    elseif remaining_P >= total_normal_P
        % 情况 2：需要所有正常状态 EV 参与响应
        lambda_star = -1; % 激活全部正常状态 EV
    else
        % 情况 3：通过二分查找确定临界点
        % 寻找第一个满足 cumulative_P(k) >= remaining_P 的索引 k
        k = find(cumulative_P >= remaining_P, 1, 'first');
        
        if isempty(k)
            lambda_star = 1;
        else
            % 确定 lambda_star 为当前 EV 的 S' 特征值
            lambda_star = S_prime_sorted(k);
            
            % 确保索引不越界
            if k > 1
                lambda_low  = S_prime_sorted(k-1);
            else
                lambda_low  = -1;
            end
            lambda_high = S_prime_sorted(k);
            
            % 取区间上限作为 lambda 临界值
            lambda_star = lambda_high;
            if lambda_star>=1
                lambda_star = 1;
            elseif lambda_star<=-1
                lambda_star = -1;
            end
        end
    end
end

% 备用算法实现（保持注释原样）
% function [P_agg, lambda_star, S_agg_next] = aggregateEVs(EVs, P_tar, M1, M2, M3, S_agg_current)
%     % 分离处于 LockON 和 LockOFF 状态的电动汽车
%     lockON_mask  = strcmp({EVs.state}, 'LockON');
%     lockOFF_mask = strcmp({EVs.state}, 'LockOFF');
%     P_lockON     = sum([EVs(lockON_mask).P_N]);
% 
%     % 处理正常状态（非 LockON/LockOFF）的电动汽车
%     normal_EVs = EVs(~lockON_mask & ~lockOFF_mask);
% 
%     % 生成 S' 特征值列表并按升序排列
%     S_prime_list = arrayfun(@(ev) getSPrime(ev), normal_EVs);
%     P_N_list     = [normal_EVs.P_N];
%     [S_prime_sorted, sort_idx] = sort(S_prime_list, 'ascend');
%     P_N_sorted   = P_N_list(sort_idx);
% 
%     % 生成从右向左的累积功率曲线
%     cumulative_P = zeros(1, length(P_N_sorted));
%     if ~isempty(cumulative_P)
%         cumulative_P(end) = P_N_sorted(end);
%         for i = length(P_N_sorted)-1:-1:1
%             cumulative_P(i) = cumulative_P(i+1) + P_N_sorted(i);
%         end
%     end
% 
%     % 计算 lambda 临界值
%     remaining_P = P_tar - P_lockON;
% 
%     if remaining_P <= 0
%         lambda_star = inf;
%         P_agg       = P_lockON;
%     elseif isempty(cumulative_P) || remaining_P >= cumulative_P(1)
%         lambda_star = -inf;
%         P_agg       = P_lockON + sum(P_N_sorted);
%     else
%         % 寻找最后一个满足 cumulative_P >= remaining_P 的索引
%         k = find(cumulative_P >= remaining_P, 1, 'last');
% 
%         if isempty(k)
%             lambda_star = inf;
%             P_agg       = P_lockON;
%         else
%             lambda_star      = S_prime_sorted(k);
%             activated_EVs    = k:length(P_N_sorted);
%             P_agg            = P_lockON + sum(P_N_sorted(activated_EVs));
%         end
%     end
% 
%     % 更新聚合 SOC 状态（带饱和限制）
%     S_agg_next = (P_agg - M2 * S_agg_current - M3) / M1;
%     S_agg_next = max(min(S_agg_next, 1), -1);
% end
