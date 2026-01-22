function [P_cmd, integral_error_next] = calculatePIPower_single(SOC_current, SOC_target, integral_error_prev, Kp, Ki, dt, P_max_limit, P_min_limit)
    % calculatePIPower_single: 变频空调PI控制 (高增益震荡版)
    
    % 1. 计算误差 (目标 - 当前)
    error_val = SOC_target - SOC_current;
    
    % 2. 预计算积分项
    % 积分累积速度由 dt 决定，为了加快消除误差，这里直接累积
    integral_term_trial = integral_error_prev + error_val * dt;
    
    % 3. 计算原始功率
    P_raw = Kp * error_val + Ki * integral_term_trial;
    
    % 4. 抗积分饱和 (优化版：允许在未饱和区域快速积分)
    if P_raw > P_max_limit
        P_cmd = P_max_limit;
        % 只有当误差方向试图让功率回归安全区时，才允许积分变化
        if error_val < 0 
            integral_error_next = integral_term_trial;
        else
            integral_error_next = integral_error_prev; % 冻结积分，防止发散
        end
    elseif P_raw < P_min_limit
        P_cmd = P_min_limit;
        if error_val > 0
            integral_error_next = integral_term_trial;
        else
            integral_error_next = integral_error_prev;
        end
    else
        P_cmd = P_raw;
        integral_error_next = integral_term_trial;
    end
    
    % 放宽安全钳位：允许积分项在更大范围内工作，以消除该死的静差
    % SOC范围是0-1， +/- 10.0 的积分库已经非常大了
    integral_error_next = max(-10, min(10, integral_error_next));
    
end