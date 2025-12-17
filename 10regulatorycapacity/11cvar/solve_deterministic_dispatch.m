%% --- 辅助函数：确定性调度求解器 ---
function st = solve_deterministic_dispatch(P_grid_demand, direction_signal, ...
    Reliable_AC_Up, Reliable_EV_Up, Reliable_AC_Down, Reliable_EV_Down, ...
    R_Gen_Max, R_Shed_Max, cost_params, dt, options, T_steps)

    % 变量顺序: [P_AC(1:T), P_EV(1:T), P_Gen(1:T), P_Shed(1:T)]
    n_vars = 4 * T_steps;
    
    H = sparse(n_vars, n_vars);
    f = zeros(n_vars, 1);
    
    % 构造目标函数
    idx_ac = 1:T_steps;
    idx_ev = (T_steps+1):2*T_steps;
    idx_gen = (2*T_steps+1):3*T_steps;
    idx_shed = (3*T_steps+1):4*T_steps;
    
    for t = 1:T_steps
        % AC
        i = idx_ac(t);
        H(i,i) = 2 * cost_params.c2_ac * dt;
        f(i) = cost_params.c1_ac * dt;
        % EV
        i = idx_ev(t);
        H(i,i) = 2 * cost_params.c2_ev * dt;
        f(i) = cost_params.c1_ev * dt;
        % Gen
        i = idx_gen(t);
        H(i,i) = 2 * cost_params.c2_gen * dt;
        f(i) = cost_params.c1_gen * dt;
        % Shed
        i = idx_shed(t);
        H(i,i) = 2 * cost_params.c2_shed * dt; % usually 0
        f(i) = cost_params.c1_shed * dt;
    end
    
    % 构造等式约束 (功率平衡)
    % P_AC + P_EV + P_Gen + P_Shed = P_Demand
    Aeq = sparse(T_steps, n_vars);
    beq = P_grid_demand;
    
    for t = 1:T_steps
        Aeq(t, idx_ac(t)) = 1;
        Aeq(t, idx_ev(t)) = 1;
        Aeq(t, idx_gen(t)) = 1;
        Aeq(t, idx_shed(t)) = 1;
    end
    
    % 构造上下限约束 (确定性边界)
    lb = zeros(n_vars, 1);
    ub = zeros(n_vars, 1);
    
    for t = 1:T_steps
        if direction_signal(t) == 1 % Up
            ub(idx_ac(t)) = Reliable_AC_Up(t);
            ub(idx_ev(t)) = Reliable_EV_Up(t);
        else % Down
            ub(idx_ac(t)) = abs(Reliable_AC_Down(t));
            ub(idx_ev(t)) = abs(Reliable_EV_Down(t));
        end
        ub(idx_gen(t)) = R_Gen_Max(t);
        ub(idx_shed(t)) = R_Shed_Max(t);
    end
    
    [x, ~, exitflag] = quadprog(H, f, [], [], Aeq, beq, lb, ub, [], options);
    
    if exitflag > 0
        st.P_AC = x(idx_ac);
        st.P_EV = x(idx_ev);
        st.P_Gen = x(idx_gen);
        st.P_Shed = x(idx_shed);
    else
        warning('方法 1 (确定性) 求解失败');
        st.P_AC = zeros(T_steps,1); st.P_EV = zeros(T_steps,1);
        st.P_Gen = zeros(T_steps,1); st.P_Shed = zeros(T_steps,1);
    end
end