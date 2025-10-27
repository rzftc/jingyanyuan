function [MI_plus, MI_minus, Constraints] = calculateMutualInfo(n_AC, n_EV, AC_plus, AC_minus, EV_plus, EV_minus, N_AC_max, N_EV_max)
    % 参数设置
    K = 3;                          % 离散化区间数
    T_val = size(AC_plus,1);        % 时段总数（数值）
    M = 1e4;                        % 大M常数
    epsilon = 1e-6;                 % 防重叠常数
    n_points = 20;                  % 对数函数分段线性化点数
    
    % ================== 上调节互信息计算 ==================
    % 1. 定义调节能力变量
    Q_AC_plus = sdpvar(T_val,1);
    Q_EV_plus = sdpvar(T_val,1);
    Constraints = [
        Q_AC_plus == sum(n_AC.*AC_plus,2), 
        Q_EV_plus == sum(n_EV.*EV_plus,2)
    ];

    % 2. 离散化设置
    Q_max_plus = max(N_AC_max*max(AC_plus(:)), N_EV_max*max(EV_plus(:)));
    bin_edges = linspace(0, Q_max_plus*(1+epsilon), K+1);
    
    % 3. 创建二元变量（修正维度）
    z_plus = binvar(T_val, K, K, 'full');
    
    % 4. 构建分箱约束
    for t = 1:T_val
        Constraints = [Constraints, sum(sum(z_plus(t,:,:))) == 1];
        for i = 1:K
            for j = 1:K
                Constraints = [Constraints,
                    Q_AC_plus(t) >= bin_edges(i) - M*(1 - z_plus(t,i,j)),
                    Q_EV_plus(t) >= bin_edges(j) - M*(1 - z_plus(t,i,j)),
                    Q_AC_plus(t) <= bin_edges(i+1) + M*(1 - z_plus(t,i,j)) - epsilon,
                    Q_EV_plus(t) <= bin_edges(j+1) + M*(1 - z_plus(t,i,j)) - epsilon
                ];
            end
        end
    end
    
    % 5. 概率计算（修正维度匹配）
    P_plus = sdpvar(K,K,'full');
    sum_z_plus = sum(z_plus,1);  % 1×K×K
    
    % 将三维sum_z_plus转换为二维
    sum_z_plus_2d = reshape(sum_z_plus, K, K);
    Constraints = [Constraints, 
        P_plus == sum_z_plus_2d/T_val
    ];
    
    % 边际概率计算
    P_AC_plus = sdpvar(K,1);
    P_EV_plus = sdpvar(1,K);
    Constraints = [Constraints,
        P_AC_plus == sum(P_plus,2),
        P_EV_plus == sum(P_plus,1)
    ];
    
    % 6. 对数函数线性化
    log_points = linspace(epsilon, 1, n_points);
    log_values = log2(log_points);
    
    % 定义对数变量
    log_P_plus = sdpvar(K,K,'full');
    log_P_AC = sdpvar(K,1);
    log_P_EV = sdpvar(1,K);
    
    % 添加分段线性约束
    for i = 1:K
        for j = 1:K
            % 创建选择器变量
            selector = binvar(n_points,1);
            Constraints = [Constraints,
                sum(selector) == 1,
                P_plus(i,j) == selector'*log_points',
                log_P_plus(i,j) == selector'*log_values'
            ];
        end
        
        % 边际概率对数
        selector_AC = binvar(n_points,1);
        selector_EV = binvar(n_points,1);
        Constraints = [Constraints,
            sum(selector_AC) == 1,
            P_AC_plus(i) == selector_AC'*log_points',
            log_P_AC(i) == selector_AC'*log_values',
            
            sum(selector_EV) == 1,
            P_EV_plus(i) == selector_EV'*log_points',
            log_P_EV(i) == selector_EV'*log_values'
        ];
    end
    
    % 7. 互信息计算（修正维度）
    MI_plus = sum(sum(P_plus .* (log_P_plus - repmat(log_P_AC,1,K) - repmat(log_P_EV,K,1))));
    
    % ================== 下调节互信息计算 ==================
    % 1. 定义调节能力变量
    Q_AC_minus = sum(n_AC.*abs(AC_minus),2);
    Q_EV_minus = sum(n_EV.*abs(EV_minus),2);
    
    % 2. 离散化设置
    Q_max_minus = max(N_AC_max*max(abs(AC_minus(:))), N_EV_max*max(abs(EV_minus(:))));
    bin_edges_minus = linspace(0, Q_max_minus*(1+epsilon), K+1);
    
    % 3. 创建二元变量
    z_minus = binvar(T_val, K, K, 'full');
    
    % 4. 构建分箱约束
    for t = 1:T_val
        Constraints = [Constraints, sum(sum(z_minus(t,:,:))) == 1];
        for i = 1:K
            for j = 1:K
                Constraints = [Constraints,
                    Q_AC_minus(t) >= bin_edges_minus(i) - M*(1 - z_minus(t,i,j)),
                    Q_EV_minus(t) >= bin_edges_minus(j) - M*(1 - z_minus(t,i,j)),
                    Q_AC_minus(t) <= bin_edges_minus(i+1) + M*(1 - z_minus(t,i,j)) - epsilon,
                    Q_EV_minus(t) <= bin_edges_minus(j+1) + M*(1 - z_minus(t,i,j)) - epsilon
                ];
            end
        end
    end
    
    % 5. 概率计算
    P_minus = sdpvar(K,K,'full');
    sum_z_minus = sum(z_minus,1);
    sum_z_minus_2d = reshape(sum_z_minus, K, K);
    Constraints = [Constraints, 
        P_minus == sum_z_minus_2d/T_val
    ];
    
    % 边际概率
    P_AC_minus = sdpvar(K,1);
    P_EV_minus = sdpvar(1,K);
    Constraints = [Constraints,
        P_AC_minus == sum(P_minus,2),
        P_EV_minus == sum(P_minus,1)
    ];
    
    % 6. 对数近似
    log_P_minus = sdpvar(K,K,'full');
    log_P_AC_minus = sdpvar(K,1);
    log_P_EV_minus = sdpvar(1,K);
    
    for i = 1:K
        for j = 1:K
            selector = binvar(n_points,1);
            Constraints = [Constraints,
                sum(selector) == 1,
                P_minus(i,j) == selector'*log_points',
                log_P_minus(i,j) == selector'*log_values'
            ];
        end
        
        selector_AC = binvar(n_points,1);
        selector_EV = binvar(n_points,1);
        Constraints = [Constraints,
            sum(selector_AC) == 1,
            P_AC_minus(i) == selector_AC'*log_points',
            log_P_AC_minus(i) == selector_AC'*log_values',
            
            sum(selector_EV) == 1,
            P_EV_minus(i) == selector_EV'*log_points',
            log_P_EV_minus(i) == selector_EV'*log_values'
        ];
    end
    
    % 7. 互信息计算
    MI_minus = sum(sum(P_minus .* (log_P_minus - repmat(log_P_AC_minus,1,K) - repmat(log_P_EV_minus,K,1))));
    
    % ================== 全局约束 ==================
    Constraints = [Constraints,
        % 概率归一化
        sum(P_plus(:)) == 1, sum(P_minus(:)) == 1,
        % 非负约束
        P_plus >= epsilon, P_minus >= epsilon
    ];
end
