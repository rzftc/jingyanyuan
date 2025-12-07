function [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast(...
    P_req, S_AC, S_EV, R_AC, R_EV, cost_p, risk_p, net_p)
    % construct_risk_constrained_qp_fast (最终修正版 - 矢量化构建)
    % 解决大规模场景下构建 A 矩阵过慢的问题
    
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
    
    %% 3. 不等式约束 (矢量化加速版)
    rho = risk_p.rho_pen;
    
    % --- 3.1 CVaR 约束组 1: Loss - eta - z <= 0 (矢量化构建) ---
    % 原始逻辑：
    % A_cvar1(s, idx_P_AC) = rho; 
    % A_cvar1(s, idx_P_EV) = rho;
    % A_cvar1(s, idx_eta)  = -1;
    % A_cvar1(s, idx_z(s)) = -1;
    
    % 1. 构造 P_AC 部分的索引和值
    % 行索引: 1..N 重复 T 次 -> [1,1.., 2,2.., N,N..]
    rows_pac = reshape(repmat((1:N_scenarios)', 1, T)', [], 1);
    % 列索引: 1..T 重复 N 次 -> [1,2..T, 1,2..T, ...]
    cols_pac = repmat(idx_P_AC(:), N_scenarios, 1);
    vals_pac = repmat(rho, length(rows_pac), 1);
    
    % 2. 构造 P_EV 部分的索引和值 (逻辑同上)
    rows_pev = rows_pac; 
    cols_pev = repmat(idx_P_EV(:), N_scenarios, 1);
    vals_pev = repmat(rho, length(rows_pev), 1);
    
    % 3. 构造 eta 部分 (每行 s，在 idx_eta 列上值为 -1)
    rows_eta = (1:N_scenarios)';
    cols_eta = repmat(idx_eta, N_scenarios, 1);
    vals_eta = repmat(-1, N_scenarios, 1);
    
    % 4. 构造 z 部分 (对角线: 每行 s，在 idx_z(s) 列上值为 -1)
    rows_z = (1:N_scenarios)';
    cols_z = idx_z(:); 
    vals_z = repmat(-1, N_scenarios, 1);
    
    % 合并三元组构建稀疏矩阵 A_cvar1 (一次性构建，极大提升速度)
    A_cvar1 = sparse([rows_pac; rows_pev; rows_eta; rows_z], ...
                     [cols_pac; cols_pev; cols_eta; cols_z], ...
                     [vals_pac; vals_pev; vals_eta; vals_z], ...
                     N_scenarios, num_vars);
                     
    % 计算右端项 b_cvar1 (矢量化计算)
    % total_limit_s = sum(S_AC(:,s)) + sum(S_EV(:,s));
    b_cvar1 = rho * (sum(S_AC, 1)' + sum(S_EV, 1)');

    % --- 3.2 CVaR 约束组 2: -eta - z <= 0 (z >= -eta) ---
    % 原始逻辑：
    % A_cvar2(s, idx_eta)  = -1;
    % A_cvar2(s, idx_z(s)) = -1;
    
    % 复用上面的 eta 和 z 的索引与值
    i_cvar2 = [rows_eta; rows_z];
    j_cvar2 = [cols_eta; cols_z];
    v_cvar2 = [vals_eta; vals_z];
    
    A_cvar2 = sparse(i_cvar2, j_cvar2, v_cvar2, N_scenarios, num_vars);
    b_cvar2 = zeros(N_scenarios, 1);
    
    % --- 3.3 网络潮流约束 (保持原逻辑，因不随场景数 N 增加) ---
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
    
    % 合并所有约束
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