function [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_robust_network(...
    P_req, S_AC, S_EV, R_AC, R_EV, cost_p, risk_p, net_p)
    % construct_risk_constrained_qp_robust_network
    % 构建包含网络潮流约束(PTDF)的风险约束QP模型
    %
    % 新增输入 net_p 结构体包含:
    %   net_p.PTDF: [N_line * N_bus] 矩阵
    %   net_p.LineLimit: [N_line * 1] 线路容量上限
    %   net_p.BaseFlow: [N_line * T] 线路基础潮流 (可选，默认为0)
    %   net_p.AcDist: [N_bus * 1] AC功率在各节点的分布比例 (和为1)
    %   net_p.EvDist: [N_bus * 1] EV功率在各节点的分布比例 (和为1)

    [T, N_scenarios] = size(S_AC);
    N_line = size(net_p.PTDF, 1);
    
    % --- 变量索引定义 (保持不变) ---
    % X = [P_AC(T); P_EV(T); P_slack(T); eta(1); z(N)]
    idx_P_AC   = 1:T;
    idx_P_EV   = (T+1):2*T;
    idx_Slack  = (2*T+1):3*T;
    idx_eta    = 3*T + 1;
    idx_z      = (3*T + 2) : (3*T + 1 + N_scenarios);
    
    num_vars = idx_z(end);
    
    %% 1. 目标函数 (保持不变)
    H = sparse(num_vars, num_vars);
    H(sub2ind(size(H), idx_P_AC, idx_P_AC)) = 2 * cost_p.c2_ac;
    H(sub2ind(size(H), idx_P_EV, idx_P_EV)) = 2 * cost_p.c2_ev;
    H(sub2ind(size(H), idx_Slack, idx_Slack)) = 2 * 10000; 
    
    f = zeros(num_vars, 1);
    f(idx_P_AC) = cost_p.c1_ac;
    f(idx_P_EV) = cost_p.c1_ev;
    f(idx_Slack) = cost_p.c_slack; 
    f(idx_eta)  = risk_p.beta;
    f(idx_z)    = risk_p.beta / (N_scenarios * (1 - risk_p.confidence));
    
    %% 2. 等式约束 (功率平衡) (保持不变)
    Aeq = sparse(T, num_vars);
    for t = 1:T
        Aeq(t, idx_P_AC(t)) = 1;
        Aeq(t, idx_P_EV(t)) = 1;
        Aeq(t, idx_Slack(t)) = 1; 
    end
    beq = P_req;
    
    %% 3. 不等式约束 A*x <= b
    
    % --- 3.1 CVaR 约束 (原有的) ---
    rho = risk_p.rho_pen;
    
    % 约束组 1: Loss_s - eta - z_s <= 0
    A_cvar1 = sparse(N_scenarios, num_vars);
    b_cvar1 = zeros(N_scenarios, 1);
    for s = 1:N_scenarios
        A_cvar1(s, idx_P_AC) = rho; 
        A_cvar1(s, idx_P_EV) = rho;
        A_cvar1(s, idx_eta)  = -1;
        A_cvar1(s, idx_z(s)) = -1;
        b_cvar1(s) = rho * (sum(S_AC(:,s)) + sum(S_EV(:,s)));
    end

    % 约束组 2: -eta - z_s <= 0
    A_cvar2 = sparse(N_scenarios, num_vars);
    b_cvar2 = zeros(N_scenarios, 1);
    for s = 1:N_scenarios
        A_cvar2(s, idx_eta)  = -1;
        A_cvar2(s, idx_z(s)) = -1;
    end
    
    % --- 3.2 网络潮流约束 (PTDF) [新增] ---
    % 计算聚合灵敏度因子 (N_line x 1)
    %这是为了将聚合功率 P_AC(t) 映射到线路潮流上
    % Sens_AC(l) = sum( PTDF(l, n) * AcDist(n) )
    Sens_AC = net_p.PTDF * net_p.AcDist; 
    Sens_EV = net_p.PTDF * net_p.EvDist; 
    
    % 初始化网络约束矩阵
    % 约束数量 = 2 * N_line * T (每条线每个时刻有上下限)
    num_net_con = 2 * N_line * T;
    A_net = sparse(num_net_con, num_vars);
    b_net = zeros(num_net_con, 1);
    
    count = 0;
    for t = 1:T
        % 获取当前时刻的基础潮流
        if isfield(net_p, 'BaseFlow') && size(net_p.BaseFlow, 2) == T
            P_base_t = net_p.BaseFlow(:, t);
        else
            P_base_t = zeros(N_line, 1);
        end
        
        % 假设 VPP 是提供上调服务 (减少负荷/增加注入)
        % Flow = BaseFlow + PTDF * (P_inj)
        % Flow = BaseFlow + Sens_AC * P_AC(t) + Sens_EV * P_EV(t)
        
        % 约束 1: Flow <= LineLimit
        % Sens_AC * P_AC(t) + Sens_EV * P_EV(t) <= LineLimit - BaseFlow
        for l = 1:N_line
            count = count + 1;
            A_net(count, idx_P_AC(t)) = Sens_AC(l);
            A_net(count, idx_P_EV(t)) = Sens_EV(l);
            b_net(count) = net_p.LineLimit(l) - P_base_t(l);
        end
        
        % 约束 2: Flow >= -LineLimit  =>  -Flow <= LineLimit
        % -Sens_AC * P_AC(t) - Sens_EV * P_EV(t) <= LineLimit + BaseFlow
        for l = 1:N_line
            count = count + 1;
            A_net(count, idx_P_AC(t)) = -Sens_AC(l);
            A_net(count, idx_P_EV(t)) = -Sens_EV(l);
            b_net(count) = net_p.LineLimit(l) + P_base_t(l);
        end
    end
    
    % 合并所有不等式约束
    A = [A_cvar1; A_cvar2; A_net];
    b = [b_cvar1; b_cvar2; b_net];
    
    %% 4. 变量边界 (保持不变)
    lb = -inf(num_vars, 1);
    ub = inf(num_vars, 1);
    
    lb(idx_P_AC) = 0; ub(idx_P_AC) = R_AC;
    lb(idx_P_EV) = 0; ub(idx_P_EV) = R_EV;
    lb(idx_Slack) = 0; ub(idx_Slack) = inf; 
    
    lb(idx_eta) = -inf; 
    lb(idx_z) = 0;      
    
    %% 输出索引
    info.idx_P_AC = idx_P_AC;
    info.idx_P_EV = idx_P_EV;
    info.idx_Slack = idx_Slack;
end