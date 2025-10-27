function [MI_plus, MI_minus, Constraints] = calculateMutualInfo2(n_AC, n_EV, AC_plus, AC_minus, EV_plus, EV_minus, N_AC_max, N_EV_max)
    % 参数设置
    K = 3;                          % 离散化区间数
    T_val = size(AC_plus,1);        % 时段总数
    M = 1e4;                        % 大M常数
    epsilon = 1e-6;                 % 防重叠常数
    n_points = 8;                   % 分段线性化点数
    
    % ================== 合并计算框架 ==================
    % 1. 合并调节能力计算
    Q_AC = {sum(n_AC.*AC_plus,2), sum(n_AC.*abs(AC_minus),2)}; % 上下调节
    Q_EV = {sum(n_EV.*EV_plus,2), sum(n_EV.*abs(EV_minus),2)};
    
    % 2. 统一离散化设置
    Q_max = max([N_AC_max*max(AC_plus(:)), N_EV_max*max(EV_plus(:)),...
                N_AC_max*max(abs(AC_minus(:))), N_EV_max*max(abs(EV_minus(:)))]);
    bin_edges = linspace(0, Q_max*(1+epsilon), K+1);
    
    % 3. 创建合并的二元变量（T×K×K×2）
    z = binvar(T_val, K, K, 2, 'full'); 
    
    % 4. 构建联合分箱约束
    Constraints = [];
    for t = 1:T_val
        for k = 1:2 % 上下调节
            Constraints = [Constraints, sum(sum(z(t,:,:,k))) == 1];
            for i = 1:K
                for j = 1:K
                    Constraints = [Constraints,
                        Q_AC{k}(t) >= bin_edges(i) - M*(1 - z(t,i,j,k)),
                        Q_EV{k}(t) >= bin_edges(j) - M*(1 - z(t,i,j,k)),
                        Q_AC{k}(t) <= bin_edges(i+1) + M*(1 - z(t,i,j,k)) - epsilon,
                        Q_EV{k}(t) <= bin_edges(j+1) + M*(1 - z(t,i,j,k)) - epsilon
                    ];
                end
            end
        end
    end
    
    % 5. 概率计算优化（K×K×2）
    P = sdpvar(K,K,2,'full'); 
    for k = 1:2
        sum_z = sum(z(:,:,:,k),1);
        Constraints = [Constraints, P(:,:,k) == reshape(sum_z, K, K)/T_val];
    end
    
    % 6. 对数近似优化（关键修正）
    log_points = linspace(epsilon, 1, n_points)';
    log_values = log2(log_points);
    
    % 定义所有对数变量
    log_P = sdpvar(K,K,2,'full');       % log(P)
    log_P_AC = sdpvar(K,2);             % log(sum(P,2))
    log_P_EV = sdpvar(K,2);             % log(sum(P,1))
    
    % 预计算选择器变量
    selector = binvar(n_points, K, K, 2, 'full');
    Constraints = [Constraints, sum(selector,1) == 1];
    
    % 向量化约束
    for k = 1:2
        for i = 1:K
            for j = 1:K
                % P(i,j,k)的线性化
                Constraints = [Constraints,
                    P(i,j,k) == selector(:,i,j,k)'*log_points,
                    log_P(i,j,k) == selector(:,i,j,k)'*log_values
                ];
            end
            
            % 边际概率对数（每个k和i/j只计算一次）
            if j == 1 % P_AC约束
                selector_AC = binvar(n_points,1);
                Constraints = [Constraints,
                    sum(selector_AC) == 1,
                    sum(P(i,:,k)) == selector_AC'*log_points,
                    log_P_AC(i,k) == selector_AC'*log_values
                ];
            end
            
            if i == 1 % P_EV约束
                selector_EV = binvar(n_points,1);
                Constraints = [Constraints,
                    sum(selector_EV) == 1,
                    sum(P(:,j,k)) == selector_EV'*log_points,
                    log_P_EV(j,k) == selector_EV'*log_values
                ];
            end
        end
    end
    
    % 7. 互信息计算（修正维度）
    MI_plus = sum(sum(P(:,:,1) .* (log_P(:,:,1) - repmat(log_P_AC(:,1),1,K) - repmat(log_P_EV(:,1)',K,1))));
    MI_minus = sum(sum(P(:,:,2) .* (log_P(:,:,2) - repmat(log_P_AC(:,2),1,K) - repmat(log_P_EV(:,2)',K,1))));
    
    % ================== 全局约束 ==================
    Constraints = [Constraints,
        sum(P(:)) == 2,  % 合并归一化约束
        P >= epsilon
    ];
end
