function [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast(...
    P_req, S_AC, S_EV, R_AC, R_EV, R_Gen, R_Shed, cost_p, risk_p, net_p)
    % construct_risk_constrained_qp_fast (物理约束增强版)
    
    [T, N_scenarios] = size(S_AC);
    N_line = size(net_p.PTDF, 1);
    
    % --- 1. 变量索引 ---
    idx_P_AC   = 1:T;
    idx_P_EV   = (T+1):2*T;
    idx_P_Gen  = (2*T+1):3*T;   % 火电调节
    idx_P_Shed = (3*T+1):4*T;   % 切负荷
    idx_eta    = 4*T + 1;
    idx_z      = (4*T + 2) : (4*T + 1 + N_scenarios);
    
    num_vars = idx_z(end);
    
    %% 2. 目标函数
    H = sparse(num_vars, num_vars);
    H(sub2ind(size(H), idx_P_AC, idx_P_AC))   = 2 * cost_p.c2_ac + 1e-6;
    H(sub2ind(size(H), idx_P_EV, idx_P_EV))   = 2 * cost_p.c2_ev + 1e-6;
    H(sub2ind(size(H), idx_P_Gen, idx_P_Gen)) = 2 * cost_p.c2_gen + 1e-6;
    H(sub2ind(size(H), idx_P_Shed, idx_P_Shed)) = 1e-6; 
    H(idx_eta, idx_eta) = 1e-6; 
    H(sub2ind(size(H), idx_z, idx_z)) = 1e-6;

    f = zeros(num_vars, 1);
    f(idx_P_AC)   = cost_p.c1_ac;
    f(idx_P_EV)   = cost_p.c1_ev;
    f(idx_P_Gen)  = cost_p.c1_gen;
    f(idx_P_Shed) = cost_p.c1_shed;
    f(idx_eta)    = risk_p.beta;
    f(idx_z)      = risk_p.beta / (N_scenarios * (1 - risk_p.confidence));
    
    %% 3. 等式约束
    Aeq = sparse(T, num_vars);
    for t = 1:T
        Aeq(t, idx_P_AC(t))   = 1;
        Aeq(t, idx_P_EV(t))   = 1;
        Aeq(t, idx_P_Gen(t))  = 1;
        Aeq(t, idx_P_Shed(t)) = 1;
    end
    beq = P_req;
    
    %% 4. 不等式约束
    rho = risk_p.rho_pen;
    
    % CVaR 约束 (AC/EV)
    rows_pac = reshape(repmat((1:N_scenarios)', 1, T)', [], 1);
    cols_pac = repmat(idx_P_AC(:), N_scenarios, 1);
    vals_pac = repmat(rho, length(rows_pac), 1);
    rows_pev = rows_pac; 
    cols_pev = repmat(idx_P_EV(:), N_scenarios, 1);
    vals_pev = repmat(rho, length(rows_pev), 1);
    rows_eta = (1:N_scenarios)'; cols_eta = repmat(idx_eta, N_scenarios, 1); vals_eta = repmat(-1, N_scenarios, 1);
    rows_z = (1:N_scenarios)'; cols_z = idx_z(:); vals_z = repmat(-1, N_scenarios, 1);
    
    A_cvar1 = sparse([rows_pac; rows_pev; rows_eta; rows_z], ...
                     [cols_pac; cols_pev; cols_eta; cols_z], ...
                     [vals_pac; vals_pev; vals_eta; vals_z], ...
                     N_scenarios, num_vars);
    b_cvar1 = rho * (sum(S_AC, 1)' + sum(S_EV, 1)');

    i_cvar2 = [rows_eta; rows_z]; j_cvar2 = [cols_eta; cols_z]; v_cvar2 = [vals_eta; vals_z];
    A_cvar2 = sparse(i_cvar2, j_cvar2, v_cvar2, N_scenarios, num_vars);
    b_cvar2 = zeros(N_scenarios, 1);
    
    % 网络约束
    Sens_AC = net_p.PTDF * net_p.AcDist; Sens_EV = net_p.PTDF * net_p.EvDist;
    Sens_Gen = net_p.PTDF * net_p.GenDist; Sens_Shed = net_p.PTDF * net_p.ShedDist;
    
    num_net_con = 2 * N_line * T;
    A_net = sparse(num_net_con, num_vars);
    b_net = zeros(num_net_con, 1);
    
    count = 0;
    for t = 1:T
        P_base_t = zeros(N_line, 1);
        if isfield(net_p, 'BaseFlow'), P_base_t = net_p.BaseFlow(:, t); end
        
        for l = 1:N_line
            count = count + 1;
            A_net(count, idx_P_AC(t))   = Sens_AC(l); A_net(count, idx_P_EV(t))   = Sens_EV(l);
            A_net(count, idx_P_Gen(t))  = Sens_Gen(l); A_net(count, idx_P_Shed(t)) = Sens_Shed(l);
            b_net(count) = net_p.LineLimit(l) - P_base_t(l);
            
            count = count + 1;
            A_net(count, idx_P_AC(t))   = -Sens_AC(l); A_net(count, idx_P_EV(t))   = -Sens_EV(l);
            A_net(count, idx_P_Gen(t))  = -Sens_Gen(l); A_net(count, idx_P_Shed(t)) = -Sens_Shed(l);
            b_net(count) = net_p.LineLimit(l) + P_base_t(l);
        end
    end
    
    A = [A_cvar1; A_cvar2; A_net];
    b = [b_cvar1; b_cvar2; b_net];
    
    %% 5. 变量边界 (物理约束核心)
    lb = -inf(num_vars, 1);
    ub = inf(num_vars, 1);
    
    lb(idx_P_AC) = 0; ub(idx_P_AC) = R_AC;
    lb(idx_P_EV) = 0; ub(idx_P_EV) = R_EV;
    
    % [重要] 火电和切负荷的物理边界 (此处接收向量并正确赋值)
    lb(idx_P_Gen) = 0;  ub(idx_P_Gen)  = R_Gen; 
    lb(idx_P_Shed) = 0; ub(idx_P_Shed) = R_Shed;
    
    lb(idx_eta) = -1e10; 
    lb(idx_z) = 0;      
    
    %% 输出索引
    info.idx_P_AC   = idx_P_AC;
    info.idx_P_EV   = idx_P_EV;
    info.idx_P_Gen  = idx_P_Gen;
    info.idx_P_Shed = idx_P_Shed;
    info.idx_Slack  = idx_P_Shed;
    info.idx_z      = idx_z;
    info.idx_eta    = idx_eta;
end