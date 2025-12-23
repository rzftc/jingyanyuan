function [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast_ramp_tly(...
    P_req, S_AC, S_EV_Scenarios, R_AC, R_EV, R_Gen, R_Shed, cost_p, risk_p, net_p)
% construct_risk_constrained_qp_fast_ramp_tly 
% 最终修正版：
% 1. 采用【增量模型 (Deviation Model)】消除基线积分漂移，恢复上调能力。
% 2. 保留软约束 (Soft Constraints) 以处理调节过程中的累积偏差。

    [T, N_scenarios] = size(S_AC);
    N_line = size(net_p.PTDF, 1);
    
    % 提取参数
    ac_p = net_p.AC_Params;
    ev_p = net_p.EV_Params; 
    dir_sig = net_p.Direction_Signal; 

    %% 1. 变量索引
    idx_P_AC   = 1:T;
    idx_P_EV   = (T+1):2*T;
    idx_P_Gen  = (2*T+1):3*T;
    idx_P_Shed = (3*T+1):4*T;
    idx_eta    = 4*T + 1;
    idx_z      = (4*T + 2) : (4*T + 1 + N_scenarios);
    
    % 状态变量
    idx_S_AC   = (idx_z(end) + 1) : (idx_z(end) + T);
    idx_S_EV   = (idx_S_AC(end) + 1) : (idx_S_AC(end) + T);
    
    % 松弛变量 (用于 EV 状态软约束)
    idx_Slk_EV = (idx_S_EV(end) + 1) : (idx_S_EV(end) + T);
    
    % 变量总数
    num_vars = idx_Slk_EV(end);

    %% 2. 目标函数
    H = sparse(num_vars, num_vars);
    H(sub2ind(size(H), idx_P_AC,   idx_P_AC))   = 2 * cost_p.c2_ac   + 1e-6;
    H(sub2ind(size(H), idx_P_EV,   idx_P_EV))   = 2 * cost_p.c2_ev   + 1e-6;
    H(sub2ind(size(H), idx_P_Gen,  idx_P_Gen))  = 2 * cost_p.c2_gen  + 1e-6;
    H(sub2ind(size(H), idx_P_Shed, idx_P_Shed)) = 1e-6;
    H(idx_eta, idx_eta)                         = 1e-6;
    H(sub2ind(size(H), idx_z, idx_z))           = 1e-6;
    
    % 状态正则
    H(sub2ind(size(H), idx_S_AC, idx_S_AC))     = 1e-9;
    H(sub2ind(size(H), idx_S_EV, idx_S_EV))     = 1e-9;
    
    % 松弛变量二次正则
    H(sub2ind(size(H), idx_Slk_EV, idx_Slk_EV)) = 1e-6;

    f = zeros(num_vars, 1);
    f(idx_P_AC)   = cost_p.c1_ac;
    f(idx_P_EV)   = cost_p.c1_ev;
    f(idx_P_Gen)  = cost_p.c1_gen;
    f(idx_P_Shed) = cost_p.c1_shed;
    f(idx_eta)    = risk_p.beta;
    f(idx_z)      = risk_p.beta / (N_scenarios * (1 - risk_p.confidence));
    
    % 对松弛变量施加惩罚 (保证尽可能不越界)
    Penalty_Weight = 1e6; 
    f(idx_Slk_EV) = Penalty_Weight;

    %% 3. 等式约束
    
    % --- 3.1 功率平衡 ---
    num_bal_con = T;
    
    % --- 3.2 状态演化 (AC + EV) ---
    num_state_ac = T; 
    num_state_ev = T; 
    total_eq = num_bal_con + num_state_ac + num_state_ev;
    
    Aeq = sparse(total_eq, num_vars);
    beq = zeros(total_eq, 1);
    
    S_AC_0 = 0.5; 
    S_EV_0 = 0;   
    
    % 3.1 功率平衡
    for t = 1:T
        Aeq(t, idx_P_AC(t))   = 1;
        Aeq(t, idx_P_EV(t))   = 1;
        Aeq(t, idx_P_Gen(t))  = 1;
        Aeq(t, idx_P_Shed(t)) = 1;
        beq(t) = P_req(t);
    end
    
    % 3.2 AC 状态约束
    row_offset = T;
    ac_A = ac_p.A; ac_B = ac_p.B; ac_C = ac_p.C;
    for t = 1:T
        row_ac = row_offset + t;
        Aeq(row_ac, idx_S_AC(t)) = 1;
        if t > 1
            Aeq(row_ac, idx_S_AC(t-1)) = -ac_A;
            rhs_ac = ac_C;
        else
            rhs_ac = ac_A * S_AC_0 + ac_C;
        end
        coef_p_ac = ac_B * dir_sig(t);
        Aeq(row_ac, idx_P_AC(t)) = -coef_p_ac;
        beq(row_ac) = rhs_ac;
    end
    
    % 3.3 [修改] EV 状态约束 (改为增量模型)
    % 物理含义：S_EV 代表相对于基线的【状态偏差量 Delta_S】
    % 公式：Delta_S(t) = A(t)*Delta_S(t-1) + B(t)*P_dispatch(t)
    % 注意：C(t) 和 P_base(t) 被数学抵消，从而消除了基线漂移。
    
    row_offset = 2 * T;
    ev_A = ev_p.A; ev_B = ev_p.B; 
    % [注意] 不需要 ev_C 和 P_base_ev 了
    
    for t = 1:T
        row_ev = row_offset + t;
        at = ev_A(t); bt = ev_B(t); 
        
        Aeq(row_ev, idx_S_EV(t)) = 1;
        
        % 右侧 constant term (rhs) 设为 0
        % 这确保了当 P_dispatch=0 时，S_EV (偏差) 保持为 0，不会漂移。
        if t > 1
            Aeq(row_ev, idx_S_EV(t-1)) = -at;
            rhs_ev = 0; 
        else
            rhs_ev = 0; % 初始偏差假设为 0
        end
        
        coef_p_ev = bt * dir_sig(t);
        Aeq(row_ev, idx_P_EV(t)) = -coef_p_ev;
        beq(row_ev) = rhs_ev;
    end

    %% 4. 不等式约束
    rho = risk_p.rho_pen;
    
    % --- 4.1 CVaR 约束 ---
    num_con_cvar = T * N_scenarios;
    row_idx = (1:num_con_cvar)';
    t_idx = repmat((1:T)', N_scenarios, 1); 
    s_idx = ceil(row_idx / T); 
    
    i_pac = row_idx; j_pac = idx_P_AC(t_idx)'; v_pac = repmat(rho, num_con_cvar, 1);
    i_pev = row_idx; j_pev = idx_P_EV(t_idx)'; v_pev = repmat(rho, num_con_cvar, 1);
    i_eta = row_idx; j_eta = repmat(idx_eta, num_con_cvar, 1); v_eta = repmat(-1, num_con_cvar, 1);
    i_z = row_idx; j_z = idx_z(s_idx)'; v_z = repmat(-1, num_con_cvar, 1);
    
    A_cvar1 = sparse([i_pac; i_pev; i_eta; i_z], [j_pac; j_pev; j_eta; j_z], [v_pac; v_pev; v_eta; v_z], num_con_cvar, num_vars);
    
    R_AC_bound = zeros(num_con_cvar, 1);
    R_EV_bound = zeros(num_con_cvar, 1);
    for s = 1:N_scenarios
        start_idx = (s-1)*T + 1;
        end_idx = s*T;
        R_AC_bound(start_idx:end_idx) = S_AC(:, s);
        R_EV_bound(start_idx:end_idx) = S_EV_Scenarios(:, s);
    end
    b_cvar1 = rho * (R_AC_bound + R_EV_bound);

    rows_eta = (1:N_scenarios)'; cols_eta = repmat(idx_eta, N_scenarios, 1); vals_eta = repmat(-1, N_scenarios, 1);
    rows_z = (1:N_scenarios)'; cols_z = idx_z(:); vals_z = repmat(-1, N_scenarios, 1);
    A_cvar2 = sparse([rows_eta; rows_z], [cols_eta; cols_z], [vals_eta; vals_z], N_scenarios, num_vars);
    b_cvar2 = zeros(N_scenarios, 1);

    % --- 4.2 PTDF 约束 ---
    Sens_AC = net_p.PTDF * net_p.AcDist(:);
    Sens_EV = net_p.PTDF * net_p.EvDist(:);
    Sens_Gen = zeros(N_line, 1); if isfield(net_p, 'GenDist'), Sens_Gen = net_p.PTDF * net_p.GenDist(:); end
    Sens_Shed = zeros(N_line, 1); if isfield(net_p, 'ShedDist'), Sens_Shed = net_p.PTDF * net_p.ShedDist(:); end

    num_net_con = 2 * N_line * T;
    A_net = sparse(num_net_con, num_vars);
    b_net = zeros(num_net_con, 1);
    count = 0;
    for t = 1:T
        if isfield(net_p, 'BaseFlow') && size(net_p.BaseFlow, 2) >= t, P_base_t = net_p.BaseFlow(:, t); else, P_base_t = zeros(N_line, 1); end
        for l = 1:N_line
            count = count + 1;
            A_net(count, idx_P_AC(t)) = Sens_AC(l); A_net(count, idx_P_EV(t)) = Sens_EV(l); A_net(count, idx_P_Gen(t)) = Sens_Gen(l); A_net(count, idx_P_Shed(t)) = Sens_Shed(l);
            b_net(count) = net_p.LineLimit(l) - P_base_t(l);
            count = count + 1;
            A_net(count, idx_P_AC(t)) = -Sens_AC(l); A_net(count, idx_P_EV(t)) = -Sens_EV(l); A_net(count, idx_P_Gen(t)) = -Sens_Gen(l); A_net(count, idx_P_Shed(t)) = -Sens_Shed(l);
            b_net(count) = net_p.LineLimit(l) + P_base_t(l);
        end
    end

    % --- 4.3 爬坡约束 ---
    A_ramp = sparse(0, num_vars); b_ramp = zeros(0, 1);
    if T > 1
        if isfield(risk_p, 'ramp_ac') && ~isempty(risk_p.ramp_ac), ramp_ac = risk_p.ramp_ac(:); else, ramp_ac = 999 * ones(T-1, 1); end
        if isfield(risk_p, 'ramp_ev') && ~isempty(risk_p.ramp_ev), ramp_ev = risk_p.ramp_ev(:); else, ramp_ev = 999 * ones(T-1, 1); end
        if isfield(risk_p, 'ramp_gen') && ~isempty(risk_p.ramp_gen), ramp_gen = risk_p.ramp_gen(:); else, ramp_gen = 999 * ones(T-1, 1); end
        
        max_ramp_con = 2 * 3 * (T - 1);
        A_ramp = sparse(max_ramp_con, num_vars);
        b_ramp = zeros(max_ramp_con, 1);
        k = 0;
        for t = 2:T
             tt = t-1;
             k=k+1; A_ramp(k, idx_P_AC(t))=1; A_ramp(k, idx_P_AC(t-1))=-1; b_ramp(k)=ramp_ac(min(tt,end));
             k=k+1; A_ramp(k, idx_P_AC(t))=-1; A_ramp(k, idx_P_AC(t-1))=1; b_ramp(k)=ramp_ac(min(tt,end));
             k=k+1; A_ramp(k, idx_P_EV(t))=1; A_ramp(k, idx_P_EV(t-1))=-1; b_ramp(k)=ramp_ev(min(tt,end));
             k=k+1; A_ramp(k, idx_P_EV(t))=-1; A_ramp(k, idx_P_EV(t-1))=1; b_ramp(k)=ramp_ev(min(tt,end));
             k=k+1; A_ramp(k, idx_P_Gen(t))=1; A_ramp(k, idx_P_Gen(t-1))=-1; b_ramp(k)=ramp_gen(min(tt,end));
             k=k+1; A_ramp(k, idx_P_Gen(t))=-1; A_ramp(k, idx_P_Gen(t-1))=1; b_ramp(k)=ramp_gen(min(tt,end));
        end
        A_ramp = A_ramp(1:k, :); b_ramp = b_ramp(1:k);
    end
    
    % --- 4.4 EV 状态软约束 ---
    % 约束对象：偏差量 Delta_S
    % 范围：[-1.2, 1.2] (归一化偏差)
    
    Limit_Val = 1.2; 
    
    num_soft_con = 2 * T;
    A_soft = sparse(num_soft_con, num_vars);
    b_soft = zeros(num_soft_con, 1);
    
    for t = 1:T
        % 上限: Delta_S - Slk <= 1.2
        row_u = t;
        A_soft(row_u, idx_S_EV(t))   = 1;
        A_soft(row_u, idx_Slk_EV(t)) = -1;
        b_soft(row_u) = Limit_Val;
        
        % 下限: -Delta_S - Slk <= 1.2
        row_l = T + t;
        A_soft(row_l, idx_S_EV(t))   = -1;
        A_soft(row_l, idx_Slk_EV(t)) = -1;
        b_soft(row_l) = Limit_Val;
    end

    A = [A_cvar1; A_cvar2; A_net; A_ramp; A_soft];
    b = [b_cvar1; b_cvar2; b_net; b_ramp; b_soft];

    %% 5. 变量边界
    lb = -inf(num_vars, 1); ub = inf(num_vars, 1);
    
    lb(idx_P_AC) = 0; ub(idx_P_AC) = R_AC(:);
    lb(idx_P_EV) = 0; ub(idx_P_EV) = R_EV(:);
    lb(idx_P_Gen) = 0; ub(idx_P_Gen) = R_Gen(:);
    lb(idx_P_Shed) = 0; ub(idx_P_Shed) = R_Shed(:);
    lb(idx_eta) = -1e10; lb(idx_z) = 0;

    lb(idx_S_AC) = 0; ub(idx_S_AC) = 1;
    
    % EV 状态 (偏差) 无硬边界
    lb(idx_S_EV) = -inf; ub(idx_S_EV) = inf;
    
    % 松弛变量非负
    lb(idx_Slk_EV) = 0; ub(idx_Slk_EV) = inf;

    %% 6. 输出信息
    info.idx_P_AC = idx_P_AC; info.idx_P_EV = idx_P_EV; info.idx_P_Gen = idx_P_Gen; info.idx_P_Shed = idx_P_Shed;
    info.idx_eta = idx_eta; info.idx_z = idx_z;
    info.idx_S_AC = idx_S_AC; 
    info.idx_S_EV = idx_S_EV; 
end