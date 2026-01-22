function [U_ac_optimal, U_ev_optimal, U_pv_optimal, optimal_total_cost, exitflag, output, population, scores] = ...
    solve_total_time_dispatch_ga_pv(num_ac, num_ev, num_pv, ...
                                 P_ac_individual_t, P_ev_individual_t, P_pv_individual_t, ...
                                 C_ac_individual, C_ev_individual, C_pv_individual, ...
                                 P_req_t, T_horizon, ...
                                 SDCI_origin, rho_origin, ga_options_in)
  
    % --- 1. 定义决策变量的总数 ---
    num_vars_ac_calc = num_ac * T_horizon;
    num_vars_ev_calc = num_ev * T_horizon;
    num_vars_pv_calc = num_pv * T_horizon;
    N_vars_to_ga = num_vars_ac_calc + num_vars_ev_calc + num_vars_pv_calc;

    if N_vars_to_ga == 0 || T_horizon == 0
        U_ac_optimal = zeros(num_ac, T_horizon);
        U_ev_optimal = zeros(num_ev, T_horizon);
        U_pv_optimal = zeros(num_pv, T_horizon);
        optimal_total_cost = 0;
        exitflag = -200; % 表示无优化发生或无法优化
        output = struct('message', 'N_vars_to_ga is 0 or T_horizon is 0. No optimization performed.');
        population = []; scores = [];
        return;
    end

    % --- 2. 定义适应度函数句柄 ---
    % 注意：调用的是支持PV的 objective_function_ga_pv
    fitness_function_handle = @(x_ga) objective_function_ga_pv(x_ga, num_ac, num_ev, num_pv, ...
                                P_ac_individual_t, P_ev_individual_t, P_pv_individual_t, ...
                                C_ac_individual, C_ev_individual, C_pv_individual, ...
                                P_req_t, T_horizon, ...
                                SDCI_origin, rho_origin);

    % --- 3. 设置遗传算法参数 ---
    if nargin < 14 || isempty(ga_options_in)
        options = optimoptions('ga', ...
            'PopulationType', 'bitstring', ...
            'PopulationSize', min(max(50, 2*N_vars_to_ga), 200), ...
            'MaxGenerations', min(max(100, 10*N_vars_to_ga), 500), ...
            'EliteCount', ceil(0.05 * min(max(50, 2*N_vars_to_ga), 200)), ...
            'CrossoverFraction', 0.8, ...
            'MutationFcn', {@mutationuniform, 0.01}, ...
            'Display', 'final', ... 
            'UseParallel', false);
    else
        options = ga_options_in;
        if ~isfield(options, 'PopulationType') || ~strcmp(options.PopulationType, 'bitstring')
             options.PopulationType = 'bitstring';
        end
    end

    % --- 4. 运行遗传算法 ---
    x_sol_ga = []; fval_ga = Inf; exitflag = -1; output = struct(); population = []; scores = []; 
    
    [x_sol_ga, fval_ga, exitflag, output, population, scores] = ga(fitness_function_handle, ...
                                                               N_vars_to_ga, ...
                                                               [], [], [], [], [], [], ... 
                                                               [], ... 
                                                               options);

    % --- 5. 解析并重塑优化结果 ---
    U_ac_optimal = zeros(num_ac, T_horizon);
    U_ev_optimal = zeros(num_ev, T_horizon);
    U_pv_optimal = zeros(num_pv, T_horizon);
    optimal_total_cost = fval_ga;

    if ~isempty(x_sol_ga) && (length(x_sol_ga) == N_vars_to_ga) 
        x_sol_reshaped = x_sol_ga(:);
        
        % 提取 AC
        if num_vars_ac_calc > 0 && num_ac > 0
            ac_vars_flat = x_sol_reshaped(1:num_vars_ac_calc);
            U_ac_optimal = reshape(logical(ac_vars_flat), T_horizon, num_ac)';
        end

        % 提取 EV
        if num_vars_ev_calc > 0 && num_ev > 0
            start_idx_ev = num_vars_ac_calc + 1;
            end_idx_ev = num_vars_ac_calc + num_vars_ev_calc;
            ev_vars_flat = x_sol_reshaped(start_idx_ev : end_idx_ev);
            U_ev_optimal = reshape(logical(ev_vars_flat), T_horizon, num_ev)';
        end

        % 提取 PV
        if num_vars_pv_calc > 0 && num_pv > 0
            start_idx_pv = num_vars_ac_calc + num_vars_ev_calc + 1;
            pv_vars_flat = x_sol_reshaped(start_idx_pv : N_vars_to_ga);
            U_pv_optimal = reshape(logical(pv_vars_flat), T_horizon, num_pv)';
        end
    end
end