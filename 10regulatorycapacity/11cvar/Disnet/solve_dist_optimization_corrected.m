function [P_AC_opt, P_EV_opt, P_slack_opt, U_opt, fval] = solve_dist_optimization_corrected(...
    P_req, net, cost_p, risk_p, S_AC, S_EV, R_AC, R_EV, lambda_SDCI, lambda_Rho, direction_vec)
% solve_dist_optimization_corrected
% 修正版求解器：支持混合方向 (上调/下调)
% 
% 输入:
%   direction_vec: T x 1 向量, 1表示上调(Up, 增负荷), -1表示下调(Down, 减负荷)
%   S_AC, S_EV: 场景数据 (需与当前时刻的方向对应, 若t时刻为Up, 则传入Up的场景)

    [T, N_scenarios] = size(S_AC);
    N_bus = net.N_bus;
    N_branch = length(net.branch_r);
    
    % 变量定义: [P_line(N_br); Q_line(N_br); U_sq(N_bus); P_AC(1); P_EV(1); P_slack(1)]
    n_vars_t = 2*N_branch + N_bus + 3;
    
    idx_start_P_line = 1;
    idx_start_Q_line = N_branch + 1;
    idx_start_U_sq   = 2*N_branch + 1;
    idx_P_AC_agg     = 2*N_branch + N_bus + 1;
    idx_P_EV_agg     = 2*N_branch + N_bus + 2;
    idx_P_slack      = 2*N_branch + N_bus + 3;
    
    num_phys_vars = T * n_vars_t;
    
    % CVaR 辅助变量
    idx_eta = num_phys_vars + 1;
    idx_z   = idx_eta + 1 : idx_eta + N_scenarios;
    num_vars_total = idx_z(end);
    
    %% 1. 目标函数 (Cost + Penalty + CVaR)
    H = sparse(num_vars_total, num_vars_total);
    f = zeros(num_vars_total, 1);
    
    for t = 1:T
        offset = (t-1) * n_vars_t;
        id_ac = offset + idx_P_AC_agg;
        id_ev = offset + idx_P_EV_agg;
        id_sl = offset + idx_P_slack;
        
        % 运行成本
        H(id_ac, id_ac) = 2 * cost_p.c2_ac;
        H(id_ev, id_ev) = 2 * cost_p.c2_ev;
        H(id_sl, id_sl) = 1e-4; % 松弛变量正则化
        
        f(id_ac) = cost_p.c1_ac;
        f(id_ev) = cost_p.c1_ev;
        f(id_sl) = cost_p.c1_shed;
        
        % 协同惩罚 (SDCI/Rho 代理)
        if lambda_SDCI > 0 || lambda_Rho > 0
            pen = lambda_SDCI + lambda_Rho;
            H(id_ac, id_ev) = pen;
            H(id_ev, id_ac) = pen;
        end
    end
    
    % CVaR 成本
    f(idx_eta) = risk_p.beta;
    f(idx_z)   = risk_p.beta / (N_scenarios * (1 - risk_p.confidence));
    
    %% 2. 等式约束 (LinDistFlow + 聚合平衡)
    % 预估非零元
    Aeq = sparse([], [], [], 0, num_vars_total, T*N_bus*15);
    beq = [];
    row = 0;
    
    for t = 1:T
        offset = (t-1) * n_vars_t;
        dir_t = direction_vec(t); % 当前时刻调节方向
        
        % --- 2.1 聚合层 PCC 追踪约束 ---
        % P_AC + P_EV + P_slack = P_req(t)
        row = row + 1;
        Aeq(row, offset + idx_P_AC_agg) = 1;
        Aeq(row, offset + idx_P_EV_agg) = 1;
        Aeq(row, offset + idx_P_slack)  = 1;
        beq(row, 1) = P_req(t);
        
        % --- 2.2 节点功率平衡 (LinDistFlow) ---
        % 决定 P_AC/EV 在节点平衡方程中的符号
        % Up (dir=1):   增加负荷 => P_load_new = P_base + P_adj
        %               Sum(P_in) - Sum(P_out) - P_adj = P_base
        %               系数为 -1
        % Down (dir=-1): 减少负荷 => P_load_new = P_base - P_adj
        %               Sum(P_in) - Sum(P_out) + P_adj = P_base
        %               系数为 +1
        
        if dir_t == 1
            res_coeff = -1;
        else
            res_coeff = 1;
        end
        
        for n = 1:N_bus
            % P 平衡
            row = row + 1;
            % 流入
            in_br = find(net.branch_t == n);
            for k = 1:length(in_br), Aeq(row, offset + idx_start_P_line + in_br(k)-1) = 1; end
            % 流出
            out_br = find(net.branch_f == n);
            for k = 1:length(out_br), Aeq(row, offset + idx_start_P_line + out_br(k)-1) = -1; end
            
            % 资源注入 (分配到节点)
            if net.Dist_AC(n) > 0
                Aeq(row, offset + idx_P_AC_agg) = res_coeff * net.Dist_AC(n);
            end
            if net.Dist_EV(n) > 0
                Aeq(row, offset + idx_P_EV_agg) = res_coeff * net.Dist_EV(n);
            end
            
            beq(row, 1) = net.P_load_base(n);
            
            % Q 平衡 (无资源支持)
            row = row + 1;
            for k = 1:length(in_br), Aeq(row, offset + idx_start_Q_line + in_br(k)-1) = 1; end
            for k = 1:length(out_br), Aeq(row, offset + idx_start_Q_line + out_br(k)-1) = -1; end
            beq(row, 1) = net.Q_load_base(n);
        end
        
        % --- 2.3 电压降方程 ---
        % U_j = U_i - 2(rP + xQ)
        for br = 1:N_branch
            f_node = net.branch_f(br);
            t_node = net.branch_t(br);
            r = net.branch_r(br);
            x = net.branch_x(br);
            
            row = row + 1;
            Aeq(row, offset + idx_start_U_sq + t_node - 1) = 1;
            Aeq(row, offset + idx_start_U_sq + f_node - 1) = -1;
            Aeq(row, offset + idx_start_P_line + br - 1)   = 2 * r;
            Aeq(row, offset + idx_start_Q_line + br - 1)   = 2 * x;
            beq(row, 1) = 0;
        end
        
        % --- 2.4 平衡节点电压 ---
        row = row + 1;
        Aeq(row, offset + idx_start_U_sq + 0) = 1; % Index 1
        beq(row, 1) = 1.0^2;
    end
    
    %% 3. 不等式约束 (CVaR)
    num_cvar = T * N_scenarios;
    A_cvar = sparse(num_cvar, num_vars_total);
    b_cvar = zeros(num_cvar, 1);
    rho = risk_p.rho_pen;
    row = 0;
    
    for t = 1:T
        offset = (t-1) * n_vars_t;
        for s = 1:N_scenarios
            row = row + 1;
            % rho*(P_AC + P_EV) - eta - z_s <= rho * Capacity(t,s)
            A_cvar(row, offset + idx_P_AC_agg) = rho;
            A_cvar(row, offset + idx_P_EV_agg) = rho;
            A_cvar(row, idx_eta) = -1;
            A_cvar(row, idx_z(s)) = -1;
            
            % 容量限制根据当前 S_AC/EV (已处理好对应方向)
            b_cvar(row) = rho * (S_AC(t,s) + S_EV(t,s));
        end
    end
    
    %% 4. 变量边界
    lb = -inf(num_vars_total, 1);
    ub =  inf(num_vars_total, 1);
    
    for t = 1:T
        offset = (t-1) * n_vars_t;
        
        % 电压安全范围 (0.95 ~ 1.05 p.u.)
       % 电压安全范围 (放宽至 0.90 ~ 1.10 p.u. 以适应 IEEE 33 基准状态)
        idx_u = offset + idx_start_U_sq : offset + idx_start_U_sq + N_bus - 1;
        lb(idx_u) = 0.90^2; % 修改为 0.90
        ub(idx_u) = 1.10^2; % 修改为 1.10 (可选)
        
        % 资源出力非负且不超过可靠域上限
        lb(offset + idx_P_AC_agg) = 0; ub(offset + idx_P_AC_agg) = R_AC(t);
        lb(offset + idx_P_EV_agg) = 0; ub(offset + idx_P_EV_agg) = R_EV(t);
        lb(offset + idx_P_slack)  = 0;
    end
    
    % CVaR z_s >= 0
    lb(idx_z) = 0;
    
    %% 5. 求解 QP
    options = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'interior-point-convex');
    [x_opt, fval, exitflag] = quadprog(H, f, A_cvar, b_cvar, Aeq, beq, lb, ub, [], options);
    
    if exitflag ~= 1
        warning('Optimization issue at some step. Exitflag: %d', exitflag);
    end
    
    %% 6. 结果提取
    P_AC_opt = zeros(T, 1);
    P_EV_opt = zeros(T, 1);
    P_slack_opt = zeros(T, 1);
    U_opt = zeros(T, N_bus);
    
    for t = 1:T
        offset = (t-1) * n_vars_t;
        P_AC_opt(t) = x_opt(offset + idx_P_AC_agg);
        P_EV_opt(t) = x_opt(offset + idx_P_EV_agg);
        P_slack_opt(t) = x_opt(offset + idx_P_slack);
        U_opt(t,:) = x_opt(offset + idx_start_U_sq : offset + idx_start_U_sq + N_bus - 1)';
    end
end