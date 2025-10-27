function [n_AC, n_EV] = optimizeComplementarity(AC_plus, AC_minus, EV_plus, EV_minus, ...
                                               mu, lambda1, lambda2, lambda3, ...
                                               N_AC, N_EV, E_EV_req, T_AC_off, ...
                                               DeltaP_grid_plus_max, DeltaP_grid_minus_max, ...
                                               theta1, theta2)
    % 获取正确的矩阵维度
    [T_AC, G] = size(AC_plus);
    [T_EV, H] = size(EV_plus);
    T = max([T_AC, T_EV, length(mu), length(DeltaP_grid_plus_max)]); % 使用最大的时间维度
    
    % 定义优化变量
    n_AC = sdpvar(G, T, 'full');
    n_EV = sdpvar(H, T, 'full');
    
    % 计算互补性指标
    [SDCI_plus, SDCI_minus, ~] = calculateComplementarity(AC_plus, AC_minus, EV_plus, EV_minus);
    
    % 计算相关性指标
    [rho_plus, rho_minus] = calculateSpearmanCorrelation(AC_plus, AC_minus, EV_plus, EV_minus);
    
    % 计算DTW距离
    [DTW_plus, DTW_minus] = calculateDTW(AC_plus, AC_minus, EV_plus, EV_minus);
    
    % 目标函数
    objective = 0;
    for t = 1:T
        % 使用安全索引访问
        current_mu = mu(min(t,length(mu)));
        for g = 1:G
            if current_mu == 1
                p_AC = AC_plus(min(t,T_AC), g);
            else
                p_AC = AC_minus(min(t,T_AC), g);
            end
            objective = objective + current_mu * n_AC(g,t) * p_AC;
        end
        
        for h = 1:H
            if current_mu == 1
                p_EV = EV_plus(min(t,T_EV), h);
            else
                p_EV = EV_minus(min(t,T_EV), h);
            end
            objective = objective + current_mu * n_EV(h,t) * p_EV;
        end
    end
    
    % 添加惩罚项
    objective = objective - lambda1 * (SDCI_plus + SDCI_minus) ...
                          - lambda2 * (rho_plus + rho_minus) ...
                          - lambda3 * (DTW_plus + DTW_minus);
    
    % 约束条件（确保所有约束都是有效的）
    constraints = [];
    
    % 设备数量限制
    for g = 1:G
        constraints = [constraints, 0 <= n_AC(g,:) <= N_AC(g)];
    end
    
    for h = 1:H
        constraints = [constraints, 0 <= n_EV(h,:) <= N_EV(h)];
    end
    
    % 用户需求约束（确保EV_minus不为空）
    if ~isempty(EV_minus) && ~isempty(n_EV)
        constraints = [constraints, ...
            sum(n_EV(:) .* abs(EV_minus(:))) <= E_EV_req];
        constraints = [constraints, ...
            sum(n_EV(:) .* abs(EV_plus(:))) <= E_EV_req];
    else
        constraints = [constraints, ...
            sum(n_EV(:)) * 0.25 <= E_EV_req];
    end
    
    constraints = [constraints, sum(n_AC(:)) * 0.25 <= T_AC_off];
    
    %% 修正后的电网安全约束
    for t = 1:T
        if t <= length(DeltaP_grid_plus_max) && t <= length(DeltaP_grid_minus_max)
            % 约束值有效性检查
            if DeltaP_grid_plus_max(t) < 0 || DeltaP_grid_minus_max(t) < 0
                error('电网约束值不能为负，时间点: %d', t);
            end
            
            % 计算各设备贡献
            if t <= T_AC
                ac_plus_term = sum(n_AC(:,t) .* AC_plus(t,:)');
                ac_minus_term = sum(n_AC(:,t) .* AC_minus(t,:)');
            else
                ac_plus_term = 0;
                ac_minus_term = 0;
            end
            
            if t <= T_EV
                ev_plus_term = sum(n_EV(:,t) .* EV_plus(t,:)');
                ev_minus_term = sum(n_EV(:,t) .* EV_minus(t,:)');
            else
                ev_plus_term = 0;
                ev_minus_term = 0;
            end
            
            % 添加有效约束
            constraints = [constraints, ...
                ac_plus_term + ev_plus_term <= DeltaP_grid_plus_max(t), ...  % 上调约束
                ac_minus_term + ev_minus_term <= DeltaP_grid_minus_max(t)];  % 下调约束（关键修改）
        end
    end
    
    % 互补性动态约束
    if ~isempty(SDCI_plus) && ~isempty(SDCI_minus)
        constraints = [constraints, ...
            SDCI_plus <= theta1, ...
            SDCI_minus <= theta2];
    end
    
    % 求解优化问题
    options = sdpsettings('solver', 'cplex', 'verbose', 1);
    optimize(constraints, -objective, options);
    
    % 返回结果
    n_AC = value(n_AC);
    n_EV = value(n_EV);
end
