% solve_total_time_dispatch_ga_ptdf.m
function [U_ac_optimal, U_ev_optimal, optimal_total_cost, exitflag, output, population, scores] = ...
    solve_total_time_dispatch_ga_ptdf(num_ac, num_ev, P_ac_individual_t, P_ev_individual_t, ...
                                 C_ac_individual, C_ev_individual, P_req_t, T_horizon, ...
                                 SDCI_origin, rho_origin, ga_options_in, ... % 同原函数
                                 Location_AC, Location_EV, PTDF_matrix, ... % 新增
                                 P_Line_Base, P_Line_Max, N_bus, N_line)    % 新增

    % --- 与 solve_total_time_dispatch_ga.m 相同的变量数计算和初始检查 ---
    num_vars_ac_calc = num_ac * T_horizon;
    num_vars_ev_calc = num_ev * T_horizon;
    N_vars_to_ga = num_vars_ac_calc + num_vars_ev_calc;

    if N_vars_to_ga == 0 || T_horizon == 0
        % (返回空结果的逻辑同原函数)
        U_ac_optimal = zeros(num_ac, T_horizon); U_ev_optimal = zeros(num_ev, T_horizon);
        optimal_total_cost = 0; exitflag = -200; output = struct('message', 'No variables or time horizon.'); population = []; scores = [];
        return;
    end

    % --- 定义适应度函数句柄 (调用新版目标函数) ---
    fitness_function_handle = @(x_ga) objective_function_ga_ptdf(x_ga, num_ac, num_ev, ...
                                P_ac_individual_t, P_ev_individual_t, ...
                                C_ac_individual, C_ev_individual, P_req_t, T_horizon, ...
                                SDCI_origin, rho_origin, ...
                                Location_AC, Location_EV, PTDF_matrix, ... % 传递新参数
                                P_Line_Base, P_Line_Max, N_bus, N_line);    % 传递新参数

    % --- 设置遗传算法参数 (同原函数) ---
    if nargin < 11 || isempty(ga_options_in)
        options = optimoptions('ga', 'PopulationType', 'bitstring', 'Display', 'final', 'UseParallel', false);
        % 可以根据需要添加更多默认选项，参考 eco_test_ga.m
    else
        options = ga_options_in;
        if ~isfield(options, 'PopulationType') || ~strcmp(options.PopulationType, 'bitstring')
             options.PopulationType = 'bitstring';
        end
    end

    % --- 运行遗传算法 (同原函数，无 nonlcon 参数) ---
    x_sol_ga = []; fval_ga = Inf; exitflag = -1; output = struct(); population = []; scores = [];
    [x_sol_ga, fval_ga, exitflag, output, population, scores] = ga(fitness_function_handle, ...
                                                               N_vars_to_ga, ...
                                                               [], [], [], [], [], [], ... % A,b,Aeq,beq,lb,ub
                                                               [], ... % nonlcon is []
                                                               options);

    % --- 解析并重塑优化结果 (同原函数) ---
    U_ac_optimal = zeros(num_ac, T_horizon);
    U_ev_optimal = zeros(num_ev, T_horizon);
    optimal_total_cost = fval_ga;

    if ~isempty(x_sol_ga) && (length(x_sol_ga) == N_vars_to_ga)
        x_sol_reshaped = x_sol_ga(:);
        if num_vars_ac_calc > 0 && num_ac > 0
            U_ac_optimal = reshape(logical(x_sol_reshaped(1:num_vars_ac_calc)), T_horizon, num_ac)';
        end
        if num_vars_ev_calc > 0 && num_ev > 0
            U_ev_optimal = reshape(logical(x_sol_reshaped(num_vars_ac_calc+1 : N_vars_to_ga)), T_horizon, num_ev)';
        end
    end
end