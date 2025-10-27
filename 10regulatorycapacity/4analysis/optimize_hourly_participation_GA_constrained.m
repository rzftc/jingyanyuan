%% 局部函数: 上层GA，优化每小时的参与设备数 (带功率约束)
function [n_ac_opt, n_ev_opt] = optimize_hourly_participation_GA_constrained(num_ac_total, num_ev_total, P_ac_hourly, P_ev_hourly, P_req_hourly, steps_per_hour, ga_options)

    nvars = 2; % [n_ac, n_ev]
    lb = [0, 0];
    ub = [num_ac_total, num_ev_total];
    intcon = 1:2;

    % 目标函数句柄 (包含功率约束检查)
    fitness_fcn = @(x) hourly_ga_fitness_constrained(x, P_ac_hourly, P_ev_hourly, P_req_hourly, steps_per_hour);

    % --- 设置 GA 选项 ---
    % 允许并行计算 (如果主循环是串行)
    ga_opts_default = optimoptions('ga', ...
        'PopulationSize', 50, ... % 可以适当增大
        'MaxGenerations', 80, ... % 可以适当增大
        'Display', 'off', ...
        'EliteCount', 3, ...
        'UseParallel', true, ... % 允许 GA 内部并行
        'ConstraintTolerance', 1e-3, ... % 放宽一点约束容忍度
        'FunctionTolerance', 1e-4);

    % 合并传入的选项和默认选项
    if nargin > 6 && ~isempty(ga_options)
        ga_opts_final = optimoptions(ga_opts_default, ga_options); % User options override defaults
    else
        ga_opts_final = ga_opts_default;
    end
    
    % --- 运行 GA ---
    try
        [x_opt, fval_opt, exitflag] = ga(fitness_fcn, nvars, [], [], [], [], lb, ub, [], intcon, ga_opts_final);
        if exitflag <= 0
             warning('小时级 GA 未找到最优解或提前终止 (exitflag=%d)。可能无法满足峰值需求。使用最大可用数量作为备选。', exitflag);
             % 备选方案：如果GA失败，返回最大数量以保证功率
             x_opt = ub;
        end
    catch ME_ga
        warning('小时级 GA 运行失败: %s。使用最大可用数量作为备选。');
        x_opt = ub; % 备选方案
    end


    n_ac_opt = round(x_opt(1));
    n_ev_opt = round(x_opt(2));
end
