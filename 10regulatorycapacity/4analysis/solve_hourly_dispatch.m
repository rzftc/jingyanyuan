function [u_ac_optimal, u_ev_optimal, optimal_cost, exitflag] = solve_hourly_dispatch(num_ac, num_ev, p_ac, p_ev, c_ac, c_ev, P_req_t)
    % SOLVE_HOURLY_DISPATCH 解决单时段虚拟电厂调度优化问题 (最小化成本版)
    % Inputs:
    %   num_ac: 空调数量
    %   num_ev: 电动汽车数量
    %   p_ac: 空调调节能力向量 (num_ac x 1), p_ac(i) 是第i台空调的调节能力 p_ACi(t)
    %   p_ev: 电动汽车调节能力向量 (num_ev x 1), p_ev(j) 是第j台EV的调节能力 p_EVj(t)
    %   c_ac: 空调单位调节能力成本向量 (num_ac x 1), c_ac(i) 是 c_ACi
    %   c_ev: 电动汽车单位调节能力成本向量 (num_ev x 1), c_ev(j) 是 c_EVj
    %   P_req_t: 电网调节需求 (标量), P_req(t)
    % Outputs:
    %   u_ac_optimal: 优化的空调参与决策 (num_ac x 1, 二进制)
    %   u_ev_optimal: 优化的电动汽车参与决策 (num_ev x 1, 二进制)
    %   optimal_cost: 优化后的最小总成本
    %   exitflag: intlinprog的退出标志

    n_vars = num_ac + num_ev;

    if n_vars == 0
        u_ac_optimal = zeros(num_ac, 1);
        u_ev_optimal = zeros(num_ev, 1);
        optimal_cost = 0;
        exitflag = -200; % 无设备
        return;
    end

    % 目标函数系数 f (用于 min f'*x)
    % 目标是: min sum(u_i * p_i * c_i)
    % 所以 f_i = p_i * c_i
    f_obj_ac_coeffs = [];
    if num_ac > 0 && ~isempty(p_ac) && ~isempty(c_ac)
        f_obj_ac_coeffs = p_ac .* c_ac; % 每台空调参与调节的总成本
    end
    
    f_obj_ev_coeffs = [];
    if num_ev > 0 && ~isempty(p_ev) && ~isempty(c_ev)
        f_obj_ev_coeffs = p_ev .* c_ev; % 每台EV参与调节的总成本
    end
    
    f_obj = [f_obj_ac_coeffs; f_obj_ev_coeffs];

    if isempty(f_obj) && n_vars > 0 % Should not happen if n_vars > 0 and p/c are provided
        warning('solve_hourly_dispatch_mincost: 目标函数系数为空，但有设备。');
        u_ac_optimal = NaN(num_ac, 1);
        u_ev_optimal = NaN(num_ev, 1);
        optimal_cost = NaN;
        exitflag = -201; 
        return;
    end

    % 不等式约束 A*x <= b
    % 原约束: sum(u_ac_i * p_ac_i) + sum(u_ev_j * p_ev_j) >= P_req_t
    % 转换后: -sum(u_ac_i * p_ac_i) - sum(u_ev_j * p_ev_j) <= -P_req_t
    A_ineq_ac_part = [];
    if num_ac > 0 && ~isempty(p_ac)
        A_ineq_ac_part = -p_ac';
    end
    A_ineq_ev_part = [];
    if num_ev > 0 && ~isempty(p_ev)
        A_ineq_ev_part = -p_ev';
    end
    A_ineq = [A_ineq_ac_part, A_ineq_ev_part];
    b_ineq = -P_req_t;

    A_eq = [];
    b_eq = [];

    lb_vars = zeros(n_vars, 1);
    ub_vars = ones(n_vars, 1);
    intcon_vars = 1:n_vars;

    options = optimoptions('intlinprog', 'Display', 'off');
    
    x_sol = NaN(n_vars, 1); 
    fval_sol = NaN;
    exitflag = -99; 
    % output_sol = struct(); % 如果需要可以取消注释


        [x_sol, fval_sol, exitflag] = intlinprog(f_obj, intcon_vars, A_ineq, b_ineq, A_eq, b_eq, lb_vars, ub_vars, [], options);
   
    
    u_ac_optimal = NaN(num_ac, 1);
    u_ev_optimal = NaN(num_ev, 1);
    optimal_cost = NaN;

    if exitflag > 0 
        if num_ac > 0
            u_ac_optimal = x_sol(1:num_ac);
        else
            u_ac_optimal = zeros(0,1);
        end
        if num_ev > 0
            u_ev_optimal = x_sol(num_ac+1:end);
        else
            u_ev_optimal = zeros(0,1);
        end
        optimal_cost = fval_sol; % 直接是最小成本值
        % fprintf('优化成功。最小总成本: %f\n', optimal_cost);
    else
        % warning('intlinprog 未找到解或提前终止, exitflag: %d', exitflag);
    end
end