function [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast_ramp_tly_origin(...
    P_req, S_AC, S_EV, R_AC, R_EV, R_Gen, R_Shed, cost_p, risk_p, net_p)
% construct_risk_constrained_qp_fast_ramp_tly 
% 修改：仅保留 AC 灰盒约束，移除 EV 状态约束。
% [新增] 增加 EV 累积能量约束
% [修正] 爬坡约束逻辑：仅在同向调节时施加约束，允许方向切换时的快速跳变 (方案A)

    [T, N_scenarios] = size(S_AC);
    N_line = size(net_p.PTDF, 1);
    
    % 提取参数
    ac_p = net_p.AC_Params;
    dir_sig = net_p.Direction_Signal; 
    
    % [新增] 提取能量边界参数
    if isfield(net_p, 'Reliable_EV_E_Up') && isfield(net_p, 'dt')
        E_Up_Bound = net_p.Reliable_EV_E_Up;
        E_Down_Bound = net_p.Reliable_EV_E_Down;
        dt = net_p.dt;
        enable_energy_constraint = true;
    else
        enable_energy_constraint = false;
        warning('未检测到能量边界参数，EV能量约束未启用。');
    end

    %% 1. 变量索引
    idx_P_AC   = 1:T;
    idx_P_EV   = (T+1):2*T;
    idx_P_Gen  = (2*T+1):3*T;
    idx_P_Shed = (3*T+1):4*T;
    idx_eta    = 4*T + 1;
    idx_z      = (4*T + 2) : (4*T + 1 + N_scenarios);
    
    % [修改] 仅保留 AC 状态变量
    idx_S_AC   = (idx_z(end) + 1) : (idx_z(end) + T);
    
    % 变量总数
    num_vars = idx_S_AC(end);
    
    %% 2. 目标函数
    H = sparse(num_vars, num_vars);
    H(sub2ind(size(H), idx_P_AC,   idx_P_AC))   = 2 * cost_p.c2_ac   + 1e-6;
    H(sub2ind(size(H), idx_P_EV,   idx_P_EV))   = 2 * cost_p.c2_ev   + 1e-6;
    H(sub2ind(size(H), idx_P_Gen,  idx_P_Gen))  = 2 * cost_p.c2_gen  + 1e-6;
    H(sub2ind(size(H), idx_P_Shed, idx_P_Shed)) = 1e-6;
    H(idx_eta, idx_eta)                         = 1e-6;
    H(sub2ind(size(H), idx_z, idx_z))           = 1e-6;
    
    % AC状态变量微小正则
    H(sub2ind(size(H), idx_S_AC, idx_S_AC))     = 1e-9;

    f = zeros(num_vars, 1);
    f(idx_P_AC)   = cost_p.c1_ac;
    f(idx_P_EV)   = cost_p.c1_ev;
    f(idx_P_Gen)  = cost_p.c1_gen;
    f(idx_P_Shed) = cost_p.c1_shed;
    f(idx_eta)    = risk_p.beta;
    f(idx_z)      = risk_p.beta / (N_scenarios * (1 - risk_p.confidence));

    %% 3. 等式约束
    
    % --- 3.1 功率平衡 ---
    num_bal_con = T;
    
    % --- 3.2 状态演化 (仅 AC) ---
    num_state_con = T; 
    total_eq = num_bal_con + num_state_con;
    
    Aeq = sparse(total_eq, num_vars);
    beq = zeros(total_eq, 1);
    
    % 初始状态
    S_AC_0 = 0.5;
    
    % 填充功率平衡 (1 : T)
    for t = 1:T
        Aeq(t, idx_P_AC(t))   = 1;
        Aeq(t, idx_P_EV(t))   = 1;
        Aeq(t, idx_P_Gen(t))  = 1;
        Aeq(t, idx_P_Shed(t)) = 1;
        beq(t) = P_req(t);
    end
    
    % 填充 AC 状态约束 (T+1 : 2T)
    row_offset = T;
    ac_A = ac_p.A; ac_B = ac_p.B; ac_C = ac_p.C;
    
    for t = 1:T
        % === AC State Constraint ===
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

    %% 4. 不等式约束
    rho = risk_p.rho_pen;
    
    % CVaR 约束
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
        R_EV_bound(start_idx:end_idx) = S_EV(:, s); 
    end
    b_cvar1 = rho * (R_AC_bound + R_EV_bound);

    rows_eta = (1:N_scenarios)'; cols_eta = repmat(idx_eta, N_scenarios, 1); vals_eta = repmat(-1, N_scenarios, 1);
    rows_z = (1:N_scenarios)'; cols_z = idx_z(:); vals_z = repmat(-1, N_scenarios, 1);
    A_cvar2 = sparse([rows_eta; rows_z], [cols_eta; cols_z], [vals_eta; vals_z], N_scenarios, num_vars);
    b_cvar2 = zeros(N_scenarios, 1);

    % PTDF 约束
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

    % --- [关键修改] 爬坡约束 (仅同向约束) ---
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
             % [修改点] 仅当 t 与 t-1 方向一致时，施加爬坡约束
             % 若方向不一致(dir_sig 跳变)，则视为过零点，不施加差分约束
             if dir_sig(t) == dir_sig(t-1)
                 tt = t-1;
                 % AC
                 k=k+1; A_ramp(k, idx_P_AC(t))=1; A_ramp(k, idx_P_AC(t-1))=-1; b_ramp(k)=ramp_ac(min(tt,end));
                 k=k+1; A_ramp(k, idx_P_AC(t))=-1; A_ramp(k, idx_P_AC(t-1))=1; b_ramp(k)=ramp_ac(min(tt,end));
                 % EV
                 k=k+1; A_ramp(k, idx_P_EV(t))=1; A_ramp(k, idx_P_EV(t-1))=-1; b_ramp(k)=ramp_ev(min(tt,end));
                 k=k+1; A_ramp(k, idx_P_EV(t))=-1; A_ramp(k, idx_P_EV(t-1))=1; b_ramp(k)=ramp_ev(min(tt,end));
                 % Gen
                 k=k+1; A_ramp(k, idx_P_Gen(t))=1; A_ramp(k, idx_P_Gen(t-1))=-1; b_ramp(k)=ramp_gen(min(tt,end));
                 k=k+1; A_ramp(k, idx_P_Gen(t))=-1; A_ramp(k, idx_P_Gen(t-1))=1; b_ramp(k)=ramp_gen(min(tt,end));
             end
        end
        % 截断未使用的预分配行
        A_ramp = A_ramp(1:k, :); b_ramp = b_ramp(1:k);
    end

    % [新增] 累积能量约束 (Cumulative Energy Constraint)
    % 约束逻辑：-E_Down(k) <= sum(dir(t)*P(t)*dt) <= E_Up(k)
    A_energy = sparse(0, num_vars); b_energy = zeros(0, 1);
    
    if enable_energy_constraint
        % 我们有 2*T 个约束：T 个上限，T 个下限
        num_energy_con = 2 * T;
        A_energy = sparse(num_energy_con, num_vars);
        b_energy = zeros(num_energy_con, 1);
        
        for k = 1:T
            % 构建累积系数向量：对于 t=1...k，系数为 dir_sig(t) * dt
            % 注意：dir_sig > 0 表示上调(充电)，能量增加；dir_sig < 0 表示下调(放电)，能量减少
            
            % 上限约束行 (Row k)
            % sum_{t=1}^k [dir(t)*dt] * P_EV(t) <= E_Up_Bound(k)
            for t = 1:k
                A_energy(k, idx_P_EV(t)) = dir_sig(t) * dt;
            end
            b_energy(k) = E_Up_Bound(k);
            
            % 下限约束行 (Row T+k)
            % sum_{t=1}^k [dir(t)*dt] * P_EV(t) >= -E_Down_Bound(k)
            % => sum_{t=1}^k [-dir(t)*dt] * P_EV(t) <= E_Down_Bound(k)
            for t = 1:k
                A_energy(T+k, idx_P_EV(t)) = -dir_sig(t) * dt;
            end
            b_energy(T+k) = E_Down_Bound(k);
        end
    end

    A = [A_cvar1; A_cvar2; A_net; A_ramp; A_energy];
    b = [b_cvar1; b_cvar2; b_net; b_ramp; b_energy];

    %% 5. 变量边界
    lb = -inf(num_vars, 1); ub = inf(num_vars, 1);
    
    lb(idx_P_AC) = 0; ub(idx_P_AC) = R_AC(:);
    lb(idx_P_EV) = 0; ub(idx_P_EV) = R_EV(:);
    lb(idx_P_Gen) = 0; ub(idx_P_Gen) = R_Gen(:);
    lb(idx_P_Shed) = 0; ub(idx_P_Shed) = R_Shed(:);
    lb(idx_eta) = -1e10; lb(idx_z) = 0;

    lb(idx_S_AC) = 0; ub(idx_S_AC) = 1;

    %% 6. 输出信息
    info.idx_P_AC = idx_P_AC; info.idx_P_EV = idx_P_EV; info.idx_P_Gen = idx_P_Gen; info.idx_P_Shed = idx_P_Shed;
    info.idx_eta = idx_eta; info.idx_z = idx_z;
    info.idx_S_AC = idx_S_AC; 
end