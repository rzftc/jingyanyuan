function [S_opt, P_base_opt] = EVbaseP_single_longstep(C_EV, eta, E_tar, E_in, t_dep, t_in, time_step, r, P_ON, S_curr, H)
    % 参数验证
    if nargin < 11
        error('缺少输入参数：必须包含H参数');
    end
    
    % 定义扩展的SOC序列（包含初始状态）
    S = sdpvar(H+1, 1);  % S[1]~S[H+1]
    P_base = sdpvar(H, 1);
    
    % 初始化约束和目标函数
    Constraints = [S(1) == S_curr];  % 初始状态约束
    Objective = sum(S(2:end).^2);     % 最小化SOC波动
    
    % 预计算总充电时段
    total_charge_time = max(0, t_dep - t_in);
    
    % 主循环处理每个时间步
    for t = 1:H
        current_time = (t-1)*time_step;  % 绝对时间计算
        online = (current_time >= t_in) && (current_time <= t_dep);
        
        if online && (total_charge_time > 0)
            % 动态参数计算（修复参数错误）
            
            [m1,m2,m3]=calculateEVABC_single(C_EV, eta, E_tar, E_in, t_dep, t_in, time_step, r);
            
            % 状态方程约束（修正时间索引）
            Constraints = [Constraints, 
                P_base(t) == m1*S(t+1) + m2*S(t) + m3,  % 式(20)
                -1 <= S(t+1) <= 1,                        % SOC边界
                0 <= P_base(t) <= P_ON                    % 功率边界
            ];
        else
            % 非充电时段处理
            Constraints = [Constraints,
                P_base(t) == 0,                           % 无充放电
                S(t+1) == S(t)                            % SOC保持
            ];
        end
    end
    
    % 求解优化问题
   options = sdpsettings('solver','cplex','verbose',0);
    optimize(Constraints, Objective, options);
    
    % 结果处理
    if ~isempty(value(S))
        S_opt = value(S(2:end));
        P_base_opt = value(P_base);
    else
        S_opt = [];
        P_base_opt = [];
    end
end