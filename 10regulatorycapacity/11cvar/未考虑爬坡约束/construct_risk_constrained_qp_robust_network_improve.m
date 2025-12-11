function [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_robust_network_improve(...
    P_req, S_AC, S_EV, R_AC, R_EV, cost_p, risk_p, net_p)
    % construct_risk_constrained_qp_robust_network (最终修正版)
    % 输入:
    %   R_AC, R_EV: 必须传入物理最大容量(Physical Max)，而非保守的可靠容量，
    %               否则无法体现风险权衡。
    
    [T, N_scenarios] = size(S_AC);
    N_line = size(net_p.PTDF, 1);
    
    % --- 变量索引定义 ---
    % X = [P_AC(T); P_EV(T); P_slack(T); eta(1); z(N)]
    idx_P_AC   = 1:T;
    idx_P_EV   = (T+1):2*T;
    idx_Slack  = (2*T+1):3*T;
    idx_eta    = 3*T + 1;
    idx_z      = (3*T + 2) : (3*T + 1 + N_scenarios);
    
    num_vars = idx_z(end);
    
    %% 1. 目标函数
    H = sparse(num_vars, num_vars);
    
    % 功率成本二次项
    H(sub2ind(size(H), idx_P_AC, idx_P_AC)) = 2 * cost_p.c2_ac + 1e-6;
    H(sub2ind(size(H), idx_P_EV, idx_P_EV)) = 2 * cost_p.c2_ev + 1e-6;
    
    % [重要修复] Slack 变量只有线性惩罚，二次项设为微小正则化
    % 之前这里是 20000，导致 Slack 成本被错误放大
    H(sub2ind(size(H), idx_Slack, idx_Slack)) = 1e-6; 
    
    % 风险变量正则化 (防止 Beta=0 时奇异)
    H(idx_eta, idx_eta) = 1e-6; 
    H(sub2ind(size(H), idx_z, idx_z)) = 1e-6;

    f = zeros(num_vars, 1);
    f(idx_P_AC) = cost_p.c1_ac;
    f(idx_P_EV) = cost_p.c1_ev;
    f(idx_Slack) = cost_p.c_slack; % 线性惩罚
    f(idx_eta)  = risk_p.beta;
    f(idx_z)    = risk_p.beta / (N_scenarios * (1 - risk_p.confidence));
    
    %% 2. 等式约束 (功率平衡)
    Aeq = sparse(T, num_vars);
    for t = 1:T
        Aeq(t, idx_P_AC(t)) = 1;
        Aeq(t, idx_P_EV(t)) = 1;
        Aeq(t, idx_Slack(t)) = 1; 
    end
    beq = P_req;
    
    %% 3. 不等式约束
    rho = risk_p.rho_pen;
    
    % CVaR 约束组 1: Loss - eta - z <= 0
    % Loss = rho_pen * (P_total - Capacity_Scenario)
    A_cvar1 = sparse(N_scenarios, num_vars);
    b_cvar1 = zeros(N_scenarios, 1);
    
    for s = 1:N_scenarios
        A_cvar1(s, idx_P_AC) = rho; 
        A_cvar1(s, idx_P_EV) = rho;
        A_cvar1(s, idx_eta)  = -1;
        A_cvar1(s, idx_z(s)) = -1;
        
        % 这里的 Capacity 对应场景 s 的实际能力
        total_limit_s = sum(S_AC(:,s)) + sum(S_EV(:,s));
        b_cvar1(s) = rho * total_limit_s;
    end

    % CVaR 约束组 2: z >= -eta
    A_cvar2 = sparse(N_scenarios, num_vars);
    b_cvar2 = zeros(N_scenarios, 1);
    for s = 1:N_scenarios
        A_cvar2(s, idx_eta)  = -1;
        A_cvar2(s, idx_z(s)) = -1;
    end
    
    % 网络潮流约束
    Sens_AC = net_p.PTDF * net_p.AcDist; 
    Sens_EV = net_p.PTDF * net_p.EvDist; 
    num_net_con = 2 * N_line * T;
    A_net = sparse(num_net_con, num_vars);
    b_net = zeros(num_net_con, 1);
    
    count = 0;
    for t = 1:T
        if isfield(net_p, 'BaseFlow') && size(net_p.BaseFlow, 2) == T
            P_base_t = net_p.BaseFlow(:, t);
        else
            P_base_t = zeros(N_line, 1);
        end
        for l = 1:N_line
            count = count + 1;
            A_net(count, idx_P_AC(t)) = Sens_AC(l);
            A_net(count, idx_P_EV(t)) = Sens_EV(l);
            b_net(count) = net_p.LineLimit(l) - P_base_t(l);
            
            count = count + 1;
            A_net(count, idx_P_AC(t)) = -Sens_AC(l);
            A_net(count, idx_P_EV(t)) = -Sens_EV(l);
            b_net(count) = net_p.LineLimit(l) + P_base_t(l);
        end
    end
    
    A = [A_cvar1; A_cvar2; A_net];
    b = [b_cvar1; b_cvar2; b_net];
    
    %% 4. 变量边界
    lb = -inf(num_vars, 1);
    ub = inf(num_vars, 1);
    
    lb(idx_P_AC) = 0; ub(idx_P_AC) = R_AC; % R_AC 应为物理最大容量
    lb(idx_P_EV) = 0; ub(idx_P_EV) = R_EV;
    lb(idx_Slack) = 0; ub(idx_Slack) = inf;  % [修正] Slack必须非负，代表切负荷
    
    lb(idx_eta) = -1e10; 
    lb(idx_z) = 0;      
    
    %% 输出索引
    info.idx_P_AC = idx_P_AC;
    info.idx_P_EV = idx_P_EV;
    info.idx_Slack = idx_Slack;
    info.idx_z = idx_z;
    info.idx_eta = idx_eta;
end