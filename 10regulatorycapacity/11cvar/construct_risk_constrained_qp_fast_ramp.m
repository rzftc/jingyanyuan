function [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast_ramp(...
    P_req, S_AC, S_EV, R_AC, R_EV, R_Gen, R_Shed, cost_p, risk_p, net_p)
% construct_risk_constrained_qp_fast (物理约束增强版 + 爬坡约束)
%
% 变量说明（按时间 t = 1..T）：
%   P_AC(t)   - AC 聚合调节功率 (MW，上调为正)
%   P_EV(t)   - EV 聚合调节功率 (MW，上调为正)
%   P_Gen(t)  - 火电调节功率 (MW，上调为正)
%   P_Shed(t) - 切负荷功率 (MW，上调为正，代表未满足的需求，带高惩罚)
%   eta       - CVaR 辅助变量
%   z(s)      - CVaR 场景变量，s = 1..N_scenarios
%
% 新增：爬坡约束（ramp rate）+ CVaR 容量收紧
%   对 AC / EV / Gen 三类资源，加入：
%       -R_ramp(t-1) <= P(t) - P(t-1) <= R_ramp(t-1),  t = 2..T
%   CVaR 约束中对场景容量引入收缩系数 / 裕度：
%       total_limit_s = tight_factor * (sum S_AC + sum S_EV) - margin
%   通过 risk_p.tight_factor / risk_p.margin 控制“多保守”。

    [T, N_scenarios] = size(S_AC);
    N_line = size(net_p.PTDF, 1);

    %% 1. 变量索引
    idx_P_AC   = 1:T;
    idx_P_EV   = (T+1):2*T;
    idx_P_Gen  = (2*T+1):3*T;   % 火电调节
    idx_P_Shed = (3*T+1):4*T;   % 切负荷
    idx_eta    = 4*T + 1;
    idx_z      = (4*T + 2) : (4*T + 1 + N_scenarios);

    num_vars = idx_z(end);

    %% 2. 目标函数
    H = sparse(num_vars, num_vars);
    % 二次成本 + 轻微正则项，保证 H 正定
    H(sub2ind(size(H), idx_P_AC,   idx_P_AC))   = 2 * cost_p.c2_ac   + 1e-6;
    H(sub2ind(size(H), idx_P_EV,   idx_P_EV))   = 2 * cost_p.c2_ev   + 1e-6;
    H(sub2ind(size(H), idx_P_Gen,  idx_P_Gen))  = 2 * cost_p.c2_gen  + 1e-6;
    H(sub2ind(size(H), idx_P_Shed, idx_P_Shed)) = 1e-6;   % 切负荷只加很小正则
    H(idx_eta, idx_eta)                         = 1e-6;
    H(sub2ind(size(H), idx_z, idx_z))           = 1e-6;

    f = zeros(num_vars, 1);
    f(idx_P_AC)   = cost_p.c1_ac;
    f(idx_P_EV)   = cost_p.c1_ev;
    f(idx_P_Gen)  = cost_p.c1_gen;
    f(idx_P_Shed) = cost_p.c1_shed;
    f(idx_eta)    = risk_p.beta;
    f(idx_z)      = risk_p.beta / (N_scenarios * (1 - risk_p.confidence));

    %% 3. 等式约束：每个时间步功率平衡
    %   P_AC(t) + P_EV(t) + P_Gen(t) + P_Shed(t) = P_req(t)
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

    % ---------- 4.1 CVaR 约束 (AC/EV 容量违约) ----------
    % 约束 1: rho * (sum_t P_AC(t) + sum_t P_EV(t)) - eta - z_s <= rho * sum_t Capacity_s(t)
    % 使用稀疏构造，一次性铺满所有时间步
    rows_pac = reshape(repmat((1:N_scenarios)', 1, T)', [], 1);     % (T*N_scenarios x 1)
    cols_pac = repmat(idx_P_AC(:), N_scenarios, 1);
    vals_pac = repmat(rho, length(rows_pac), 1);

    rows_pev = rows_pac;
    cols_pev = repmat(idx_P_EV(:), N_scenarios, 1);
    vals_pev = repmat(rho, length(rows_pev), 1);

    rows_eta = (1:N_scenarios)'; 
    cols_eta = repmat(idx_eta, N_scenarios, 1);
    vals_eta = repmat(-1, N_scenarios, 1);

    rows_z = (1:N_scenarios)'; 
    cols_z = idx_z(:);
    vals_z = repmat(-1, N_scenarios, 1);

    A_cvar1 = sparse([rows_pac; rows_pev; rows_eta; rows_z], ...
                     [cols_pac; cols_pev; cols_eta; cols_z], ...
                     [vals_pac; vals_pev; vals_eta; vals_z], ...
                     N_scenarios, num_vars);

    % ☆☆ 这里是关键修改：对场景能力进行“收紧” ☆☆
    % 基础容量：每个场景 s 的 AC+EV 总能力（T 步之和）
    base_cap = sum(S_AC, 1)' + sum(S_EV, 1)';    % N_scenarios x 1

    % tight_factor: <1 表示收紧；默认 1（不收紧）
    tight_factor = 1;
    if isfield(risk_p, 'tight_factor') && ~isempty(risk_p.tight_factor)
        tight_factor = risk_p.tight_factor;
    end

    % margin: 固定减去的裕度（能量单位，和 base_cap 同单位）
    margin = 0;
    if isfield(risk_p, 'margin') && ~isempty(risk_p.margin)
        margin = risk_p.margin;
    end

    % 总容量上限（被收紧后的）
    total_limit = tight_factor * base_cap - margin;
    % 不允许为负
    total_limit = max(total_limit, 0);

    b_cvar1 = rho * total_limit;

    % 约束 2: z_s >= -eta  <=>  -eta - z_s <= 0
    i_cvar2 = [rows_eta; rows_z];
    j_cvar2 = [cols_eta; cols_z];
    v_cvar2 = [vals_eta; vals_z];
    A_cvar2 = sparse(i_cvar2, j_cvar2, v_cvar2, N_scenarios, num_vars);
    b_cvar2 = zeros(N_scenarios, 1);

    % ---------- 4.2 网络潮流约束 (基于 PTDF) ----------
    % 线性化敏感度：每类资源在各线路上的单位潮流增量
    Sens_AC = net_p.PTDF * net_p.AcDist(:);   % (N_line x 1)
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

    num_net_con = 2 * N_line * T;     % 每条线每时刻上下两个约束
    A_net = sparse(num_net_con, num_vars);
    b_net = zeros(num_net_con, 1);

    count = 0;
    for t = 1:T
        % 背景潮流，包含基础负荷 + AC/EV 基线
        if isfield(net_p, 'BaseFlow') && size(net_p.BaseFlow, 2) >= t
            P_base_t = net_p.BaseFlow(:, t);
        else
            P_base_t = zeros(N_line, 1);
        end

        for l = 1:N_line
            % 上限：P_base_t(l) + DeltaP_l(t) <= LineLimit(l)
            count = count + 1;
            A_net(count, idx_P_AC(t))   =  Sens_AC(l);
            A_net(count, idx_P_EV(t))   =  Sens_EV(l);
            A_net(count, idx_P_Gen(t))  =  Sens_Gen(l);
            A_net(count, idx_P_Shed(t)) =  Sens_Shed(l);
            b_net(count) = net_p.LineLimit(l) - P_base_t(l);

            % 下限：-(P_base_t(l) + DeltaP_l(t)) <= LineLimit(l)
            count = count + 1;
            A_net(count, idx_P_AC(t))   = -Sens_AC(l);
            A_net(count, idx_P_EV(t))   = -Sens_EV(l);
            A_net(count, idx_P_Gen(t))  = -Sens_Gen(l);
            A_net(count, idx_P_Shed(t)) = -Sens_Shed(l);
            b_net(count) = net_p.LineLimit(l) + P_base_t(l);
        end
    end

    % ---------- 4.3 爬坡约束 (AC / EV / Gen) ----------
    % 约束形式：
    %   P_x(t) - P_x(t-1) <= R_x_ramp(t-1)
    %  -(P_x(t) - P_x(t-1)) <= R_x_ramp(t-1)
    % 其中 x ∈ {AC, EV, Gen}
    if T > 1
        % --- 4.3.1 自动构造/读取各资源的爬坡上限向量 (长度 T-1) ---

        % ===== AC: 优先使用 risk_p.ramp_ac，否则用 S_AC 推断 =====
        if isfield(risk_p, 'ramp_ac') && ~isempty(risk_p.ramp_ac)
            ramp_ac = risk_p.ramp_ac(:);
            if isscalar(ramp_ac)
                ramp_ac = ramp_ac * ones(T-1, 1);
            else
                if length(ramp_ac) < T-1
                    ramp_ac = [ramp_ac; repmat(ramp_ac(end), T-1 - length(ramp_ac), 1)];
                end
                ramp_ac = ramp_ac(1:T-1);
            end
        else
            % 用场景容量 S_AC (T x N_scen) 自动估计
            % 1) 用场景均值代表“典型容量”曲线
            CAP_ac = mean(S_AC, 2);        % T x 1
            % 2) 相邻时间步的自然变化
            dCAP_ac = diff(CAP_ac);        % (T-1) x 1
            base_ramp_ac = robust_abs_percentile(dCAP_ac, 95);
            % 若全零，则回退到用 R_AC 的典型值
            if base_ramp_ac <= 0
                R_AC_vec = R_AC(:);
                base_ramp_ac = 0.5 * mean(R_AC_vec(~isnan(R_AC_vec)));
                if ~isfinite(base_ramp_ac) || base_ramp_ac < 0
                    base_ramp_ac = 0;  % 极端退化情况
                end
            end
            ramp_ac = base_ramp_ac * ones(T-1, 1);
        end

        % ===== EV: 优先使用 risk_p.ramp_ev，否则用 S_EV 推断 =====
        if isfield(risk_p, 'ramp_ev') && ~isempty(risk_p.ramp_ev)
            ramp_ev = risk_p.ramp_ev(:);
            if isscalar(ramp_ev)
                ramp_ev = ramp_ev * ones(T-1, 1);
            else
                if length(ramp_ev) < T-1
                    ramp_ev = [ramp_ev; repmat(ramp_ev(end), T-1 - length(ramp_ev), 1)];
                end
                ramp_ev = ramp_ev(1:T-1);
            end
        else
            CAP_ev = mean(S_EV, 2);        % T x 1
            dCAP_ev = diff(CAP_ev);        % (T-1) x 1
            base_ramp_ev = robust_abs_percentile(dCAP_ev, 95);
            if base_ramp_ev <= 0
                R_EV_vec = R_EV(:);
                base_ramp_ev = 0.5 * mean(R_EV_vec(~isnan(R_EV_vec)));
                if ~isfinite(base_ramp_ev) || base_ramp_ev < 0
                    base_ramp_ev = 0;
                end
            end
            ramp_ev = base_ramp_ev * ones(T-1, 1);
        end

        % ===== 火电 Gen: 优先使用 risk_p.ramp_gen，否则用 R_Gen 推断 =====
        if isfield(risk_p, 'ramp_gen') && ~isempty(risk_p.ramp_gen)
            ramp_gen = risk_p.ramp_gen(:);
            if isscalar(ramp_gen)
                ramp_gen = ramp_gen * ones(T-1, 1);
            else
                if length(ramp_gen) < T-1
                    ramp_gen = [ramp_gen; repmat(ramp_gen(end), T-1 - length(ramp_gen), 1)];
                end
                ramp_gen = ramp_gen(1:T-1);
            end
        else
            R_Gen_vec = R_Gen(:);          % T x 1
            dR = diff(R_Gen_vec);          % (T-1) x 1
            base_ramp_gen = robust_abs_percentile(dR, 95);
            % 若备用容量基本不变或过小，则退回到“平均备用的 10%”
            if base_ramp_gen <= 0
                avg_Rg = mean(R_Gen_vec(~isnan(R_Gen_vec)));
                if isfinite(avg_Rg) && avg_Rg > 0
                    base_ramp_gen = 0.1 * avg_Rg;
                else
                    base_ramp_gen = 0;
                end
            end
            ramp_gen = base_ramp_gen * ones(T-1, 1);
        end

        % --- 4.3.2 根据 ramp_* 构造线性约束矩阵 ---
        % 每个时间步 (从 2 开始) 每类资源 2 个约束
        max_ramp_con = 2 * 3 * (T - 1);
        A_ramp = sparse(max_ramp_con, num_vars);
        b_ramp = zeros(max_ramp_con, 1);
        k = 0;

        for t = 2:T
            tt = t - 1;  % ramp 索引

            % --- AC 爬坡 ---
            if ramp_ac(tt) > 0
                % 上爬坡:  P_AC(t) - P_AC(t-1) <= ramp_ac(tt)
                k = k + 1;
                A_ramp(k, idx_P_AC(t))   =  1;
                A_ramp(k, idx_P_AC(t-1)) = -1;
                b_ramp(k) = ramp_ac(tt);

                % 下爬坡:  P_AC(t-1) - P_AC(t) <= ramp_ac(tt)
                k = k + 1;
                A_ramp(k, idx_P_AC(t))   = -1;
                A_ramp(k, idx_P_AC(t-1)) =  1;
                b_ramp(k) = ramp_ac(tt);
            end

            % --- EV 爬坡 ---
            if ramp_ev(tt) > 0
                k = k + 1;
                A_ramp(k, idx_P_EV(t))   =  1;
                A_ramp(k, idx_P_EV(t-1)) = -1;
                b_ramp(k) = ramp_ev(tt);

                k = k + 1;
                A_ramp(k, idx_P_EV(t))   = -1;
                A_ramp(k, idx_P_EV(t-1)) =  1;
                b_ramp(k) = ramp_ev(tt);
            end

            % --- 火电 Gen 爬坡 ---
            if ramp_gen(tt) > 0
                k = k + 1;
                A_ramp(k, idx_P_Gen(t))   =  1;
                A_ramp(k, idx_P_Gen(t-1)) = -1;
                b_ramp(k) = ramp_gen(tt);

                k = k + 1;
                A_ramp(k, idx_P_Gen(t))   = -1;
                A_ramp(k, idx_P_Gen(t-1)) =  1;
                b_ramp(k) = ramp_gen(tt);
            end
        end

        % 截取实际使用行数
        if k > 0
            A_ramp = A_ramp(1:k, :);
            b_ramp = b_ramp(1:k);
        else
            A_ramp = sparse(0, num_vars);
            b_ramp = zeros(0, 1);
        end
    else
        % 只有一个时间步则不需要爬坡约束
        A_ramp = sparse(0, num_vars);
        b_ramp = zeros(0, 1);
    end

    % ---------- 4.x 汇总所有不等式约束 ----------
    % 原有: A = [A_cvar1; A_cvar2; A_net];
    % 现在: 再追加 A_ramp
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

    lb(idx_eta) = -1e10;      % eta 可以为负
    lb(idx_z)   = 0;          % z_s >= 0

    %% 6. 输出索引信息
    info.idx_P_AC   = idx_P_AC;
    info.idx_P_EV   = idx_P_EV;
    info.idx_P_Gen  = idx_P_Gen;
    info.idx_P_Shed = idx_P_Shed;
    info.idx_eta    = idx_eta;
    info.idx_z      = idx_z;
end

%% --- 本文件内部的小工具函数：稳健地计算 |x| 的 p 分位数 ---
function val = robust_abs_percentile(x, p)
    % 去掉 NaN / Inf
    x = x(:);
    x = x(isfinite(x));
    if isempty(x)
        val = 0;
        return;
    end
    ax = abs(x);
    ax = sort(ax);
    n = numel(ax);
    idx = max(1, min(n, round(p / 100 * n)));
    val = ax(idx);
end
