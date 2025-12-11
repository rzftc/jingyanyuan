function [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast_ramp_cvar(...
    P_req, S_AC, S_EV, R_AC, R_EV, R_Gen, R_Shed, cost_p, risk_p, net_p)
% construct_risk_constrained_qp_fast_ramp (物理约束增强版 + 爬坡约束 + 逐时刻CVaR)
%
% 修改说明：
%   已将 CVaR 约束修改为针对每个时间步 t 独立施加约束，防止 EV 利用全天总能量约束在个别时刻大幅越界。
%
% 变量说明（按时间 t = 1..T）：
%   P_AC(t)   - AC 聚合调节功率 (MW，上调为正)
%   P_EV(t)   - EV 聚合调节功率 (MW，上调为正)
%   P_Gen(t)  - 火电调节功率 (MW，上调为正)
%   P_Shed(t) - 切负荷功率 (MW，上调为正)
%   eta(t)    - CVaR 辅助变量 (每个时刻 t 一个)
%   z(s, t)   - CVaR 场景变量 (每个时刻 t、每个场景 s 一个)

    [T, N_scenarios] = size(S_AC);
    N_line = size(net_p.PTDF, 1);

    %% 1. 变量索引
    idx_P_AC   = 1:T;
    idx_P_EV   = (T+1):2*T;
    idx_P_Gen  = (2*T+1):3*T;   
    idx_P_Shed = (3*T+1):4*T;   
    
    % 修改：eta 和 z 扩展为与时间相关
    idx_eta    = 4*T + (1:T);                      % T 个 eta
    idx_z      = 5*T + (1 : N_scenarios * T);      % N*T 个 z (排列顺序：t=1的所有s, t=2的所有s...)

    num_vars = idx_z(end);

    %% 2. 目标函数
    H = sparse(num_vars, num_vars);
    % 二次成本 + 轻微正则项
    H(sub2ind(size(H), idx_P_AC,   idx_P_AC))   = 2 * cost_p.c2_ac   + 1e-6;
    H(sub2ind(size(H), idx_P_EV,   idx_P_EV))   = 2 * cost_p.c2_ev   + 1e-6;
    H(sub2ind(size(H), idx_P_Gen,  idx_P_Gen))  = 2 * cost_p.c2_gen  + 1e-6;
    H(sub2ind(size(H), idx_P_Shed, idx_P_Shed)) = 1e-6;   
    H(sub2ind(size(H), idx_eta,    idx_eta))    = 1e-6;
    H(sub2ind(size(H), idx_z,      idx_z))      = 1e-6;

    f = zeros(num_vars, 1);
    f(idx_P_AC)   = cost_p.c1_ac;
    f(idx_P_EV)   = cost_p.c1_ev;
    f(idx_P_Gen)  = cost_p.c1_gen;
    f(idx_P_Shed) = cost_p.c1_shed;
    
    % CVaR 目标：Sum_t [ beta * eta_t + beta/(N(1-a)) * Sum_s z_{s,t} ]
    f(idx_eta)    = risk_p.beta;  
    f(idx_z)      = risk_p.beta / (N_scenarios * (1 - risk_p.confidence));

    %% 3. 等式约束：每个时间步功率平衡
    Aeq = sparse(T, num_vars);
    for t = 1:T
        Aeq(t, idx_P_AC(t))   = 1;
        Aeq(t, idx_P_EV(t))   = 1;
        Aeq(t, idx_P_Gen(t))  = 1;
        Aeq(t, idx_P_Shed(t)) = 1;
    end
    beq = P_req(:);

    %% 4. 不等式约束
    rho = risk_p.rho_pen;

    % ---------- 4.1 CVaR 约束 (逐时刻独立施加) ----------
    % 约束 1: rho * (P_AC(t) + P_EV(t)) - eta(t) - z(s,t) <= rho * Capacity(t,s)
    % 约束 2: -eta(t) - z(s,t) <= 0  <=> z(s,t) >= -eta(t)
    
    % 辅助向量构造：t_vec = [1,1...1, 2,2...2, ... T,T...T]' (每个重复 N_scenarios 次)
    t_vec = kron((1:T)', ones(N_scenarios, 1));
    
    % 对应的变量索引向量 (长度 N*T)
    rows_range = (1 : N_scenarios * T)';
    
    idx_P_AC_vec = idx_P_AC(t_vec)';
    idx_P_EV_vec = idx_P_EV(t_vec)';
    idx_eta_vec  = idx_eta(t_vec)';
    idx_z_vec    = idx_z(:);  % z 已经是按 t=1, t=2... 顺序排列
    
    % --- 构建 A_cvar1 ---
    % rho * P_AC + rho * P_EV - eta - z <= RHS
    cols_cvar1 = [idx_P_AC_vec; idx_P_EV_vec; idx_eta_vec; idx_z_vec];
    rows_cvar1 = [rows_range;   rows_range;   rows_range;  rows_range];
    vals_cvar1 = [repmat(rho, N_scenarios*T, 1); 
                  repmat(rho, N_scenarios*T, 1); 
                  repmat(-1,  N_scenarios*T, 1); 
                  repmat(-1,  N_scenarios*T, 1)];
              
    A_cvar1 = sparse(rows_cvar1, cols_cvar1, vals_cvar1, N_scenarios*T, num_vars);

    % --- 计算 RHS (Total Capacity Limit per step) ---
    % S_AC: T x N -> 变换为 N*T x 1 (按 t=1...排列)
    % S_AC' (:)? No. S_AC is T rows, N cols. 
    % We need [S_AC(1,:), S_AC(2,:)...] which is S_AC'(:).
    S_AC_vec = reshape(S_AC', [], 1); 
    S_EV_vec = reshape(S_EV', [], 1);
    
    base_cap_vec = S_AC_vec + S_EV_vec;

    tight_factor = 1;
    if isfield(risk_p, 'tight_factor') && ~isempty(risk_p.tight_factor)
        tight_factor = risk_p.tight_factor;
    end

    margin = 0;
    if isfield(risk_p, 'margin') && ~isempty(risk_p.margin)
        margin = risk_p.margin;
    end

    total_limit_vec = tight_factor * base_cap_vec - margin;
    total_limit_vec = max(total_limit_vec, 0);

    b_cvar1 = rho * total_limit_vec;

    % --- 构建 A_cvar2 ---
    % -eta - z <= 0
    cols_cvar2 = [idx_eta_vec; idx_z_vec];
    rows_cvar2 = [rows_range;  rows_range];
    vals_cvar2 = [repmat(-1, N_scenarios*T, 1); 
                  repmat(-1, N_scenarios*T, 1)];
              
    A_cvar2 = sparse(rows_cvar2, cols_cvar2, vals_cvar2, N_scenarios*T, num_vars);
    b_cvar2 = zeros(N_scenarios*T, 1);

    % ---------- 4.2 网络潮流约束 (基于 PTDF) ----------
    Sens_AC = net_p.PTDF * net_p.AcDist(:);
    Sens_EV = net_p.PTDF * net_p.EvDist(:);

    if isfield(net_p, 'GenDist') && ~isempty(net_p.GenDist)
        Sens_Gen = net_p.PTDF * net_p.GenDist(:);
    else
        Sens_Gen = zeros(N_line, 1);
    end

    if isfield(net_p, 'ShedDist') && ~isempty(net_p.ShedDist)
        Sens_Shed = net_p.PTDF * net_p.ShedDist(:);
    else
        Sens_Shed = zeros(N_line, 1);
    end

    num_net_con = 2 * N_line * T;
    A_net = sparse(num_net_con, num_vars);
    b_net = zeros(num_net_con, 1);

    count = 0;
    for t = 1:T
        if isfield(net_p, 'BaseFlow') && size(net_p.BaseFlow, 2) >= t
            P_base_t = net_p.BaseFlow(:, t);
        else
            P_base_t = zeros(N_line, 1);
        end

        for l = 1:N_line
            count = count + 1;
            A_net(count, idx_P_AC(t))   =  Sens_AC(l);
            A_net(count, idx_P_EV(t))   =  Sens_EV(l);
            A_net(count, idx_P_Gen(t))  =  Sens_Gen(l);
            A_net(count, idx_P_Shed(t)) =  Sens_Shed(l);
            b_net(count) = net_p.LineLimit(l) - P_base_t(l);

            count = count + 1;
            A_net(count, idx_P_AC(t))   = -Sens_AC(l);
            A_net(count, idx_P_EV(t))   = -Sens_EV(l);
            A_net(count, idx_P_Gen(t))  = -Sens_Gen(l);
            A_net(count, idx_P_Shed(t)) = -Sens_Shed(l);
            b_net(count) = net_p.LineLimit(l) + P_base_t(l);
        end
    end

    % ---------- 4.3 爬坡约束 (AC / EV / Gen) ----------
    if T > 1
        if isfield(risk_p, 'ramp_ac') && ~isempty(risk_p.ramp_ac)
            ramp_ac = risk_p.ramp_ac(:);
            if isscalar(ramp_ac), ramp_ac = ramp_ac * ones(T-1, 1); end
            if length(ramp_ac) < T-1, ramp_ac = [ramp_ac; repmat(ramp_ac(end), T-1 - length(ramp_ac), 1)]; end
            ramp_ac = ramp_ac(1:T-1);
        else
            CAP_ac = mean(S_AC, 2);
            dCAP_ac = diff(CAP_ac);
            base_ramp_ac = robust_abs_percentile(dCAP_ac, 95);
            if base_ramp_ac <= 0
                R_AC_vec = R_AC(:);
                base_ramp_ac = 0.5 * mean(R_AC_vec(~isnan(R_AC_vec)));
                if ~isfinite(base_ramp_ac) || base_ramp_ac < 0, base_ramp_ac = 0; end
            end
            ramp_ac = base_ramp_ac * ones(T-1, 1);
        end

        if isfield(risk_p, 'ramp_ev') && ~isempty(risk_p.ramp_ev)
            ramp_ev = risk_p.ramp_ev(:);
            if isscalar(ramp_ev), ramp_ev = ramp_ev * ones(T-1, 1); end
            if length(ramp_ev) < T-1, ramp_ev = [ramp_ev; repmat(ramp_ev(end), T-1 - length(ramp_ev), 1)]; end
            ramp_ev = ramp_ev(1:T-1);
        else
            CAP_ev = mean(S_EV, 2);
            dCAP_ev = diff(CAP_ev);
            base_ramp_ev = robust_abs_percentile(dCAP_ev, 95);
            if base_ramp_ev <= 0
                R_EV_vec = R_EV(:);
                base_ramp_ev = 0.5 * mean(R_EV_vec(~isnan(R_EV_vec)));
                if ~isfinite(base_ramp_ev) || base_ramp_ev < 0, base_ramp_ev = 0; end
            end
            ramp_ev = base_ramp_ev * ones(T-1, 1);
        end

        if isfield(risk_p, 'ramp_gen') && ~isempty(risk_p.ramp_gen)
            ramp_gen = risk_p.ramp_gen(:);
            if isscalar(ramp_gen), ramp_gen = ramp_gen * ones(T-1, 1); end
            if length(ramp_gen) < T-1, ramp_gen = [ramp_gen; repmat(ramp_gen(end), T-1 - length(ramp_gen), 1)]; end
            ramp_gen = ramp_gen(1:T-1);
        else
            R_Gen_vec = R_Gen(:);
            dR = diff(R_Gen_vec);
            base_ramp_gen = robust_abs_percentile(dR, 95);
            if base_ramp_gen <= 0
                avg_Rg = mean(R_Gen_vec(~isnan(R_Gen_vec)));
                if isfinite(avg_Rg) && avg_Rg > 0, base_ramp_gen = 0.1 * avg_Rg; else, base_ramp_gen = 0; end
            end
            ramp_gen = base_ramp_gen * ones(T-1, 1);
        end

        max_ramp_con = 2 * 3 * (T - 1);
        A_ramp = sparse(max_ramp_con, num_vars);
        b_ramp = zeros(max_ramp_con, 1);
        k = 0;

        for t = 2:T
            tt = t - 1;
            if ramp_ac(tt) > 0
                k = k + 1; A_ramp(k, idx_P_AC(t)) = 1; A_ramp(k, idx_P_AC(t-1)) = -1; b_ramp(k) = ramp_ac(tt);
                k = k + 1; A_ramp(k, idx_P_AC(t)) = -1; A_ramp(k, idx_P_AC(t-1)) = 1; b_ramp(k) = ramp_ac(tt);
            end
            if ramp_ev(tt) > 0
                k = k + 1; A_ramp(k, idx_P_EV(t)) = 1; A_ramp(k, idx_P_EV(t-1)) = -1; b_ramp(k) = ramp_ev(tt);
                k = k + 1; A_ramp(k, idx_P_EV(t)) = -1; A_ramp(k, idx_P_EV(t-1)) = 1; b_ramp(k) = ramp_ev(tt);
            end
            if ramp_gen(tt) > 0
                k = k + 1; A_ramp(k, idx_P_Gen(t)) = 1; A_ramp(k, idx_P_Gen(t-1)) = -1; b_ramp(k) = ramp_gen(tt);
                k = k + 1; A_ramp(k, idx_P_Gen(t)) = -1; A_ramp(k, idx_P_Gen(t-1)) = 1; b_ramp(k) = ramp_gen(tt);
            end
        end

        if k > 0, A_ramp = A_ramp(1:k, :); b_ramp = b_ramp(1:k);
        else, A_ramp = sparse(0, num_vars); b_ramp = zeros(0, 1); end
    else
        A_ramp = sparse(0, num_vars); b_ramp = zeros(0, 1);
    end

    % ---------- 4.x 汇总所有不等式约束 ----------
    A = [A_cvar1; A_cvar2; A_net; A_ramp];
    b = [b_cvar1; b_cvar2; b_net; b_ramp];

    %% 5. 变量边界
    lb = -inf(num_vars, 1);
    ub =  inf(num_vars, 1);

    R_AC   = R_AC(:);
    R_EV   = R_EV(:);
    R_Gen  = R_Gen(:);
    R_Shed = R_Shed(:);

    lb(idx_P_AC)   = 0;       ub(idx_P_AC)   = R_AC;
    lb(idx_P_EV)   = 0;       ub(idx_P_EV)   = R_EV;
    lb(idx_P_Gen)  = 0;       ub(idx_P_Gen)  = R_Gen;
    lb(idx_P_Shed) = 0;       ub(idx_P_Shed) = R_Shed;

    lb(idx_eta) = -1e10;      
    lb(idx_z)   = 0;          

    %% 6. 输出索引信息
    info.idx_P_AC   = idx_P_AC;
    info.idx_P_EV   = idx_P_EV;
    info.idx_P_Gen  = idx_P_Gen;
    info.idx_P_Shed = idx_P_Shed;
    info.idx_eta    = idx_eta;
    info.idx_z      = idx_z;
end

function val = robust_abs_percentile(x, p)
    x = x(:);
    x = x(isfinite(x));
    if isempty(x), val = 0; return; end
    ax = abs(x);
    ax = sort(ax);
    n = numel(ax);
    idx = max(1, min(n, round(p / 100 * n)));
    val = ax(idx);
end