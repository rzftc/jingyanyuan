function [MI_plus, MI_minus, Constraints] = calculateMutualInfo3(n_AC, n_EV, AC_plus, AC_minus, EV_plus, EV_minus, N_AC_max, N_EV_max)
    % 参数设置
    K = 3;                          % 离散化区间数
    T_val = size(AC_plus,1);        
    M = 1e4;                        % 大M常数
    epsilon = 1e-6;                 
    n_points = 8;                   
    
    % ================== 调节能力计算 ==================
    Q_AC = {sum(n_AC.*AC_plus,2), sum(n_AC.*abs(AC_minus),2)};
    Q_EV = {sum(n_EV.*EV_plus,2), sum(n_EV.*abs(EV_minus),2)};
    
    % ================== 优化分箱策略 ==================
    Q_max = max([N_AC_max*max(AC_plus(:)), N_EV_max*max(EV_plus(:)),...
                N_AC_max*max(abs(AC_minus(:))), N_EV_max*max(abs(EV_minus(:)))]);
    bin_edges = linspace(0, Q_max*(1+epsilon), K+1);
    
    % ================== 变量定义 ==================
    z = binvar(T_val, K, K, 2, 'full'); 
    P = sdpvar(K,K,2,'full');
    
    % ================== 精简约束构建 ==================
    Constraints = [];
    
    % 1. 合并分箱约束（减少重复边界约束）
    for t = 1:T_val
        for k = 1:2
            Constraints = [Constraints, sum(sum(z(t,:,:,k))) == 1];
            for i = 1:K
                for j = 1:K
                    % 合并同类约束
                    Constraints = [Constraints,
                        [Q_AC{k}(t); Q_EV{k}(t)] >= [bin_edges(i); bin_edges(j)] - M*(1 - z(t,i,j,k)),
                        [Q_AC{k}(t); Q_EV{k}(t)] <= [bin_edges(i+1); bin_edges(j+1)] + M*(1 - z(t,i,j,k)) - epsilon
                    ];
                end
            end
        end
    end
    
    % 2. 概率计算优化
    for k = 1:2
        Constraints = [Constraints, P(:,:,k) == reshape(sum(z(:,:,:,k),1), K, K)/T_val];
    end
    
    % 3. 对数近似优化（共享选择器变量）
    log_points = linspace(epsilon, 1, n_points)';
    log_values = log2(log_points);
    
    selector = binvar(n_points, K, K, 2, 'full');
    Constraints = [Constraints, sum(selector,1) == 1];
    
    % 统一构建对数近似约束
    for k = 1:2
        for i = 1:K
            for j = 1:K
                Constraints = [Constraints,
                    P(i,j,k) == selector(:,i,j,k)'*log_points,
                    log2(P(i,j,k)) == selector(:,i,j,k)'*log_values  % 直接使用log2关系
                ];
            end
        end
        % 边际概率计算（移除非必要的重复约束）
        Constraints = [Constraints,
            sum(P(:,:,k),2) >= epsilon,
            sum(P(:,:,k),1) >= epsilon
        ];
    end
    
    % ================== 互信息计算 ==================
    MI_plus = sum(sum(P(:,:,1).*(log2(P(:,:,1)) - repmat(log2(sum(P(:,:,1),2)),1,K) - repmat(log2(sum(P(:,:,1),1)),K,1))));
    MI_minus = sum(sum(P(:,:,2).*(log2(P(:,:,2)) - repmat(log2(sum(P(:,:,2),2)),1,K) - repmat(log2(sum(P(:,:,2),1)),K,1))));
    
    % ================== 必要约束 ==================
    Constraints = [Constraints, P >= epsilon];
end
