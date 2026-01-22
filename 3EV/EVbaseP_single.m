function [S_opt, P_base_opt] = EVbaseP_single(m1, m2, m3, P_ON, S_curr)
    % EVBASEP_SINGLE 计算单个EV的基准功率与最优SOC
    % 输入参数:
    %   m1, m2, m3 - 动态方程系数 (由calculateEVABC_single生成)
    %   P_ON       - EV额定功率
    %   S_curr     - 当前SOC (S_k)
    % 输出参数:
    %   S_opt      - 最优下一时刻SOC (S_{k+1})
    %   P_base_opt - 最优基准功率

    % 定义优化变量 (S_{k+1})
    S_next = sdpvar(1,1);  % 下一时刻SOC

    % 计算基线功率表达式
    P_base = m1 * S_next + m2 * S_curr + m3; % 关键修改：P_base为表达式

    % 目标函数: 最小化 S_{k+1}^2 (式20)
    Objective = S_next^2;

    % 约束条件
    Constraints = [...
        % SOC边界约束
        -1 <= S_next <= 1, ...
        % 功率边界约束
        0 <= P_base <= P_ON, ...
        % 防御性约束: 确保方程有解
        m1*S_next + m2*S_curr + m3 >= 0, ... 
        m1*S_next + m2*S_curr + m3 <= P_ON ...
    ];

    % 求解器配置
    options = sdpsettings('solver', 'cplex', 'verbose', 0);
    optimize(Constraints, Objective, options);

    % 结果提取
    if ~isempty(value(S_next))
        S_opt = value(S_next);
        P_base_opt = value(P_base); % 直接输出表达式值
    else
        S_opt = [];
        P_base_opt = [];
    end
end
