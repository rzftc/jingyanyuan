function [S_agg_opt, P_base_agg_opt] = EVbaseP_aggregate(EVs, S_agg_current, H, dt)
    % 预分配参数存储
    M1 = zeros(H,1);
    M2 = zeros(H,1);
    M3 = zeros(H,1);
    P_agg_max = zeros(H,1);
    num=length(EVs);
    %% 计算时变参数（核心修改）
    for t = 1:H
        for i = 1:num
            EVs(i) = updateLockState(EVs(i), t*60);
        end
        % 计算当前时间步参数
        [M1(t), M2(t), M3(t), P_agg_max(t)] = calculateAggregationCoefficients(EVs, dt);
    end
    
    %% 优化模型构建（使用时变参数）
    S_agg = sdpvar(H,1);
    P_base_agg = sdpvar(H,1);
    Objective = sum(S_agg(2:end).^2);
    Constraints = [];
    
    for t = 1:H
        if t == 1
            S_prev = S_agg_current;
        else
            S_prev = S_agg(t-1);
        end
        % 动态方程（带时变参数）
        current_P = M1(t)*S_agg(t) + M2(t)*S_prev + M3(t);
        Constraints = [Constraints, 
            -1 <= S_agg(t) <= 1,
            0 <= current_P <= P_agg_max(t),
            current_P == P_base_agg(t)
        ];
    end

    % 求解器配置
    options = sdpsettings('solver','cplex','verbose',0);
    optimize(Constraints, Objective, options);

    % 结果处理
    if ~isempty(value(S_agg))
        S_agg_opt = value(S_agg);
        P_base_agg_opt = value(P_base_agg);
    else
        % 防御性处理：返回50%功率和线性SOC
        warning('无可行解，返回保守策略');
        S_agg_opt = linspace(S_agg_current, 0, H)';
        P_base_agg_opt = 0.5 * P_agg_max * ones(H,1);
    end
end
