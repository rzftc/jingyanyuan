% optimize_hourly_participation_GA_constrained.m (已按 Task 2 修改)
function [n_ac_opt, n_ev_opt] = optimize_hourly_participation_GA_constrained_ptdf( ...
    num_ac_total, num_ev_total, P_ac_hourly, P_ev_hourly, P_req_hourly, steps_per_hour, ga_options, ...
    Location_AC, Location_EV, PTDF_matrix, P_Line_Base_hourly, P_Line_Max, N_bus, N_line)
    % *** 【Task 2 修正】: 增加了新的输入参数 (Location_AC, ..., N_line) ***

    nvars = 2; % [n_ac, n_ev]
    lb = [0, 0];
    ub = [num_ac_total, num_ev_total];
    intcon = 1:2;

    % *** 【Task 2 修正】: 传递网络参数给适应度函数 ***
    fitness_fcn = @(x) hourly_ga_fitness_constrained(x, ...
        P_ac_hourly, P_ev_hourly, P_req_hourly, steps_per_hour, ...
        Location_AC, Location_EV, PTDF_matrix, P_Line_Base_hourly, P_Line_Max, N_bus, N_line);

    % --- 设置 GA 选项 ---
    ga_opts_default = optimoptions('ga', ...
        'PopulationSize', 50, ...
        'MaxGenerations', 80, ...
        'Display', 'off', ...
        'EliteCount', 3, ...
        'UseParallel', true, ... 
        'ConstraintTolerance', 1e-3, ...
        'FunctionTolerance', 1e-4);

    if nargin > 6 && ~isempty(ga_options) && iscell(ga_options)
        ga_opts_final = optimoptions(ga_opts_default, ga_options{:}); 
    else
        ga_opts_final = ga_opts_default;
    end
    
    % --- 运行 GA ---
    try
        [x_opt, fval_opt, exitflag] = ga(fitness_fcn, nvars, [], [], [], [], lb, ub, [], intcon, ga_opts_final);
        if exitflag <= 0
             warning('小时级 GA 未找到最优解或提前终止 (exitflag=%d)。可能无法满足峰值需求。使用最大可用数量作为备选。', exitflag);
             x_opt = ub;
        end
    catch ME_ga
        warning('小时级 GA 运行失败: %s。使用最大可用数量作为备选。');
        x_opt = ub;
    end

    n_ac_opt = round(x_opt(1));
    n_ev_opt = round(x_opt(2));
end