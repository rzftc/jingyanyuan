function [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast_ramp_tly(...
    P_req, S_AC, S_EV, R_AC, R_EV, R_Gen, R_Shed, cost_p, risk_p, net_p)
% construct_risk_constrained_qp_fast_ramp_tly 
% 修改说明：
% 1. 移除了 AC 和 EV 的爬坡约束。
% 2. 仅对火电机组 (Gen) 施加物理出力爬坡约束。
% 3. 设定了合理的默认爬坡阈值 (15 MW/15min)。

    [T, N_scenarios] = size(S_AC);
    N_line = size(net_p.PTDF, 1);
    
    % 提取参数
    ac_p = net_p.AC_Params;
    dir_sig = net_p.Direction_Signal; 
    
    % 提取能量边界参数
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
    
    % AC 状态变量
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

    % --- [核心修改] 爬坡约束 (仅限制火电机组物理爬坡) ---
    A_ramp = sparse(0, num_vars); b_ramp = zeros(0, 1);
    
    if T > 1
        % 1. 设定合理的火电物理爬坡阈值
        % 依据：IEEE 30 节点系统，假设最大机组约 300+ MW。
        % 取 15 MW/15min (即 1 MW/min) 作为合理的物理爬坡限制。
        % 这代表了机组出力的平滑变化要求。
        default_gen_ramp = 15; 
        
        if isfield(risk_p, 'ramp_gen') && ~isempty(risk_p.ramp_gen)
            ramp_gen_limit = risk_p.ramp_gen(:); 
        else
            ramp_gen_limit = default_gen_ramp * ones(T-1, 1);
        end
        
        % 2. 获取基线功率 (用于计算物理总出力)
        % 注意：请确保主程序 net_params 中包含 P_Gen_Base
        if isfield(net_p, 'P_Gen_Base')
            P_base = net_p.P_Gen_Base;
        else
            % 缺失处理：默认为 0，但会打印警告
            warning('[Ramp] net_p.P_Gen_Base 缺失，火电将按照调节量爬坡计算（不准确）！');
            P_base = zeros(T, 1);
        end

        max_ramp_con = 2 * (T - 1); 
        A_ramp = sparse(max_ramp_con, num_vars);
        b_ramp = zeros(max_ramp_con, 1);
        
        k = 0;
        for t = 2:T
             % 逻辑说明：
             % P_total(t) = P_base(t) + Sign(t) * P_QP(t)
             % dir_sig(t) == 1  (系统缺负荷) => Gen 下调 => P_QP 为减量 => Sign = -1
             % dir_sig(t) == -1 (系统缺电)   => Gen 上调 => P_QP 为增量 => Sign = +1
             
             s_t = -1; 
             if dir_sig(t) == -1, s_t = 1; end
             
             s_prev = -1;
             if dir_sig(t-1) == -1, s_prev = 1; end
             
             % 计算基线变化量
             dBase = P_base(t) - P_base(t-1);
             limit_val = ramp_gen_limit(min(t-1, end));
             
             % 约束 1: P_total(t) - P_total(t-1) <= Limit
             % => (Base_t + s_t*P_t) - (Base_p + s_p*P_p) <= Limit
             % => s_t*P_t - s_p*P_p <= Limit - dBase
             k = k + 1;
             A_ramp(k, idx_P_Gen(t))   = s_t;
             A_ramp(k, idx_P_Gen(t-1)) = -s_prev;
             b_ramp(k) = limit_val - dBase;
             
             % 约束 2: P_total(t) - P_total(t-1) >= -Limit
             % => -(P_total(t) - P_total(t-1)) <= Limit
             % => -(s_t*P_t - s_p*P_p + dBase) <= Limit
             % => -s_t*P_t + s_p*P_p <= Limit + dBase
             k = k + 1;
             A_ramp(k, idx_P_Gen(t))   = -s_t;
             A_ramp(k, idx_P_Gen(t-1)) = s_prev;
             b_ramp(k) = limit_val + dBase;
        end
    end

    % [新增] 累积能量约束 (Cumulative Energy Constraint)
    A_energy = sparse(0, num_vars); b_energy = zeros(0, 1);
    
    if enable_energy_constraint
        num_energy_con = 2 * T;
        A_energy = sparse(num_energy_con, num_vars);
        b_energy = zeros(num_energy_con, 1);
        
        for k = 1:T
            % 上限约束行 (Row k)
            for t = 1:k
                A_energy(k, idx_P_EV(t)) = dir_sig(t) * dt;
            end
            b_energy(k) = E_Up_Bound(k);
            
            % 下限约束行 (Row T+k)
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