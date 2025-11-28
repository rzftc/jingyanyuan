function [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_robust(...
    P_req, S_AC, S_EV, R_AC, R_EV, cost_p, risk_p)
    % 增强版构建函数 (V2 - 修复负成本Bug)
    % 修复：增加 CVaR 非负截断约束，防止因安全裕度过大导致的负风险成本

    [T, N_scenarios] = size(S_AC);
    
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
    H(sub2ind(size(H), idx_P_AC, idx_P_AC)) = 2 * cost_p.c2_ac;
    H(sub2ind(size(H), idx_P_EV, idx_P_EV)) = 2 * cost_p.c2_ev;
    H(sub2ind(size(H), idx_Slack, idx_Slack)) = 2 * 10000; % 松弛变量二次惩罚
    
    f = zeros(num_vars, 1);
    f(idx_P_AC) = cost_p.c1_ac;
    f(idx_P_EV) = cost_p.c1_ev;
    f(idx_Slack) = cost_p.c_slack; 
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
    
    %% 3. 不等式约束 (CVaR)
    rho = risk_p.rho_pen;
    
    % 约束组 1: z_s >= Loss_s - eta  =>  Loss_s - eta - z_s <= 0
    % Loss_s = rho * (P_total - Limit_s)
    % 展开: rho*P_AC + rho*P_EV - eta - z_s <= rho * Limit_s
    A_cvar1 = sparse(N_scenarios, num_vars);
    b_cvar1 = zeros(N_scenarios, 1);
    
    for s = 1:N_scenarios
        A_cvar1(s, idx_P_AC) = rho; 
        A_cvar1(s, idx_P_EV) = rho;
        A_cvar1(s, idx_eta)  = -1;
        A_cvar1(s, idx_z(s)) = -1;
        
        total_limit_s = sum(S_AC(:,s)) + sum(S_EV(:,s));
        b_cvar1(s) = rho * total_limit_s;
    end

    % 约束组 2 (新增 - Fix Bug): z_s >= 0 - eta  => -eta - z_s <= 0
    % 这确保了我们计算的是 max(0, Loss) 的 CVaR，避免因 Loss 为负数导致的负成本
    A_cvar2 = sparse(N_scenarios, num_vars);
    b_cvar2 = zeros(N_scenarios, 1);
    for s = 1:N_scenarios
        A_cvar2(s, idx_eta)  = -1;
        A_cvar2(s, idx_z(s)) = -1;
        % 右端项为 0
    end
    
    % 合并不等式约束
    A = [A_cvar1; A_cvar2];
    b = [b_cvar1; b_cvar2];
    
    %% 4. 变量边界
    lb = -inf(num_vars, 1);
    ub = inf(num_vars, 1);
    
    lb(idx_P_AC) = 0; ub(idx_P_AC) = R_AC;
    lb(idx_P_EV) = 0; ub(idx_P_EV) = R_EV;
    lb(idx_Slack) = 0; ub(idx_Slack) = inf; 
    
    lb(idx_eta) = -inf; % eta 可以是任意值
    lb(idx_z) = 0;      % z 必须非负
    
    %% 输出索引
    info.idx_P_AC = idx_P_AC;
    info.idx_P_EV = idx_P_EV;
    info.idx_Slack = idx_Slack;
end