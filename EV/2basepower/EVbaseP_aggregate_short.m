function [S_agg_opt, P_base_agg_opt] = EVbaseP_aggregate_short(EVs, S_agg_current, H, dt_long)
    % 输入参数说明:
    % H - 预测时域（小时）
    % dt_long - 长时间步长（分钟）
    
    %% 参数转换
    steps_per_hour = 60/dt_long;  % 每小时步数（例如dt_long=10分钟时为6步）
    total_steps = H * steps_per_hour;  % 总预测步数（24*6=144步）
    
    %% 预分配参数存储
    M1 = zeros(total_steps,1);
    M2 = zeros(total_steps,1);
    M3 = zeros(total_steps,1);
    P_agg_max = zeros(total_steps,1);
    num = length(EVs);
    
    %% 计算时变参数（按实际步长计算）
    for t = 1:total_steps
        current_time = (t-1)*dt_long;  % 当前时间（分钟）
        
        % 更新所有EV状态
        for i = 1:num
            EVs(i) = updateLockState(EVs(i), current_time);
        end
        
        % 计算当前步参数
        [M1(t), M2(t), M3(t), P_agg_max(t)] = calculateAggregationCoefficients(EVs, dt_long);
    end
    
    %% 优化模型构建
    S_agg = sdpvar(total_steps,1);
    P_base_agg = sdpvar(total_steps,1);
    Objective = sum(S_agg(2:end).^2);
    Constraints = [];
    
    for t = 1:total_steps
        if t == 1
            S_prev = S_agg_current;
        else
            S_prev = S_agg(t-1);
        end
        
        current_P = M1(t)*S_agg(t) + M2(t)*S_prev + M3(t);
        Constraints = [Constraints, 
            -1 <= S_agg(t) <= 1,
            0 <= current_P <= P_agg_max(t),
            current_P == P_base_agg(t)
        ];
    end
    
    %% 求解优化问题
    options = sdpsettings('solver','cplex','verbose',0);
    optimize(Constraints, Objective, options);
    
    %% 结果处理
    if ~isempty(value(S_agg))
        S_agg_opt = value(S_agg);
        P_base_agg_opt = value(P_base_agg);
    else
        warning('无可行解，返回保守策略');
        S_agg_opt = linspace(S_agg_current, 0, total_steps)';
        P_base_agg_opt = 0.5 * P_agg_max;
    end
end
