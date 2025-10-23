% Filename: dispatch_vpp_emergency_horizon.m
function results = dispatch_vpp_emergency_horizon(P_emergency_need_series, ac_params, ev_params)
    % 对VPP在给定时间视界内的紧急控制指令进行优化分解。
    % 该函数假设AC和EV聚合体在每个时段的最小/最大调节功率已知。
    %
    % 输入:
    %   P_emergency_need_series: (1xN_T double) 电网紧急控制总需求时序 (kW)。
    %                            正值代表VPP需要向电网供电或等效地减少负荷。
    %   ac_params: 结构体，包含AC聚合体的参数:
    %       .cost_c1: (scalar double) 线性成本系数 (元/kW)
    %       .cost_c2: (scalar double) 二次成本系数 (元/kW^2)
    %       .P_min_series: (1xN_T double) AC在各时段的最小调节功率 (kW)
    %       .P_max_series: (1xN_T double) AC在各时段的最大调节功率 (kW)
    %   ev_params: 结构体，包含EV聚合体的参数:
    %       .cost_c1: (scalar double) 线性成本系数 (元/kW)
    %       .cost_c2: (scalar double) 二次成本系数 (元/kW^2)
    %       .P_min_series: (1xN_T double) EV在各时段的最小调节功率 (kW)
    %       .P_max_series: (1xN_T double) EV在各时段的最大调节功率 (kW)
    %
    % 输出:
    %   results: 结构体，存储每个时段的优化结果:
    %       .P_AC_dispatch_series: (1xN_T double) 优化后分配给AC的功率 (kW)
    %       .P_EV_dispatch_series: (1xN_T double) 优化后分配给EV的功率 (kW)
    %       .lambda_series: (1xN_T double) 各时段的系统边际成本 (元/kW)
    %       .cost_per_step_series: (1xN_T double) 各时段的最小总调节成本 (元)
    %       .total_event_cost: (scalar double) 整个紧急事件的总调节成本 (元)
    %       .status_series: (1xN_T double) 各时段的求解状态标志 (1=成功)

    N_T = length(P_emergency_need_series); % 总时段数

    % 输入参数维度校验 (简化)
    if ~(length(ac_params.P_min_series) == N_T && length(ac_params.P_max_series) == N_T && ...
           length(ev_params.P_min_series) == N_T && length(ev_params.P_max_series) == N_T)
        error('所有输入的时序参数 (P_emergency_need_series, P_min_series, P_max_series for AC/EV) 长度必须等于 N_T.');
    end

    % 初始化结果存储
    results.P_AC_dispatch_series = zeros(1, N_T);
    results.P_EV_dispatch_series = zeros(1, N_T);
    results.lambda_series = NaN(1, N_T);
    results.cost_per_step_series = zeros(1, N_T);
    results.status_series = zeros(1, N_T);

    fprintf('开始逐时段优化分解VPP紧急控制指令 (共 %d 个时段)...\n', N_T);

    for i = 1:N_T
        P_need_current_step = P_emergency_need_series(i);

        % 当前时段AC的参数
        c_ac1 = ac_params.cost_c1;
        c_ac2 = ac_params.cost_c2;
        P_ac_min_current = ac_params.P_min_series(i);
        P_ac_max_current = ac_params.P_max_series(i);

        % 当前时段EV的参数
        c_ev1 = ev_params.cost_c1;
        c_ev2 = ev_params.cost_c2;
        P_ev_min_current = ev_params.P_min_series(i);
        P_ev_max_current = ev_params.P_max_series(i);
        
        % --- 若当前时段需求 <= 0 (假设VPP主要提供正向功率支撑) ---
        if P_need_current_step <= 0 
            results.P_AC_dispatch_series(i) = 0;
            results.P_EV_dispatch_series(i) = 0;
            results.cost_per_step_series(i) = 0; 
            results.lambda_series(i) = min(c_ac1, c_ev1); % 此时边际成本可认为是两者中较小者
            results.status_series(i) = 1; 
            if P_need_current_step < 0
                 fprintf('信息: 时段 %d, 电网需求为负 (%.2f kW)。VPP按0出力调度。\n', i, P_need_current_step);
            end
            continue; 
        end
        
        % --- 构建二次规划问题 ---
        % 目标函数: min (c_ac1*P_ac + c_ac2*P_ac^2) + (c_ev1*P_ev + c_ev2*P_ev^2)
        % quadprog 形式: min 0.5*x'*H_qp*x + f_qp'*x
        % x = [P_ac; P_ev]
        H_qp = diag([2*c_ac2, 2*c_ev2]); 
        f_qp = [c_ac1; c_ev1];          

        % 等式约束: P_ac + P_ev = P_need_current_step  (Aeq*x = beq)
        Aeq = [1, 1];
        beq = P_need_current_step;

        % 边界约束: P_k_min <= P_k <= P_k_max (lb <= x <= ub)
        lb = [P_ac_min_current; P_ev_min_current];
        ub = [P_ac_max_current; P_ev_max_current];
        
        % 检查边界有效性
        if any(lb > ub)
            fprintf('警告: 时段 %d 的调节能力上下限设置无效 (lb > ub)。AC:[%.2f,%.2f], EV:[%.2f,%.2f]. 跳过该时段优化。\n', ...
                    i, lb(1), ub(1), lb(2), ub(2));
            results.status_series(i) = -1; % 自定义错误码：边界无效
            results.cost_per_step_series(i) = Inf;            
            continue;
        end
        
        options_qp = optimoptions('quadprog', 'Display', 'off', 'Algorithm','interior-point-convex');
        
        [x_optimal_step, cost_val_qp, exitflag_qp, ~, lambda_struct_qp] = ...
            quadprog(H_qp, f_qp, [], [], Aeq, beq, lb, ub, [], options_qp);

        if exitflag_qp == 1 % 优化成功
            results.P_AC_dispatch_series(i) = x_optimal_step(1);
            results.P_EV_dispatch_series(i) = x_optimal_step(2);
            results.cost_per_step_series(i) = cost_val_qp; 
            
            % 获取拉格朗日乘子 (边际成本)
            if isfield(lambda_struct_qp, 'eqlin') && ~isempty(lambda_struct_qp.eqlin)
                 results.lambda_series(i) = -lambda_struct_qp.eqlin(1); % quadprog返回的乘子符号可能与定义相反
            else
                 % 备用方法计算lambda (如果求解器不直接返回或仅有一个约束时)
                 % 在KKT条件下，对于未达到边界的资源，其MC应等于lambda
                 mc_ac_at_opt = c_ac1 + 2*c_ac2*x_optimal_step(1);
                 mc_ev_at_opt = c_ev1 + 2*c_ev2*x_optimal_step(2);
                 if abs(x_optimal_step(1) - lb(1)) > 1e-5 && abs(x_optimal_step(1) - ub(1)) > 1e-5 % AC not at bound
                     results.lambda_series(i) = mc_ac_at_opt;
                 elseif abs(x_optimal_step(2) - lb(2)) > 1e-5 && abs(x_optimal_step(2) - ub(2)) > 1e-5 % EV not at bound
                     results.lambda_series(i) = mc_ev_at_opt;
                 else % Both might be at bounds, or one at bound and demand met
                     results.lambda_series(i) = max(mc_ac_at_opt, mc_ev_at_opt); % Or a more sophisticated logic
                 end
            end
            results.status_series(i) = 1; 
            
        elseif exitflag_qp == 0
             fprintf('信息: 时段 %d, quadprog 达到最大迭代次数但可能未收敛。结果可能非最优 (exitflag = %d)。\n', i, exitflag_qp);
             if ~isempty(x_optimal_step) % 即使未完全收敛，也存储结果
                results.P_AC_dispatch_series(i) = x_optimal_step(1);
                results.P_EV_dispatch_series(i) = x_optimal_step(2);
                results.cost_per_step_series(i) = cost_val_qp;
             else
                results.cost_per_step_series(i) = Inf;
             end
             results.lambda_series(i) = NaN; 
             results.status_series(i) = 0;
        else % 优化失败或无可行解 (exitflag_qp < 0)
            fprintf('错误: 时段 %d 优化求解失败或无可行解 (exitflag = %d)。需求=%.2f, AC_cap=[%.2f,%.2f], EV_cap=[%.2f,%.2f]\n', ...
                i, exitflag_qp, P_need_current_step, lb(1), ub(1), lb(2), ub(2));
            results.P_AC_dispatch_series(i) = NaN; % 或根据备用策略设置
            results.P_EV_dispatch_series(i) = NaN;
            results.lambda_series(i) = NaN;
            results.cost_per_step_series(i) = Inf;
            results.status_series(i) = exitflag_qp; 
            
            % 简单备用策略：如果总可用容量不足，则按最大能力分配
            if sum(ub) < P_need_current_step && P_need_current_step > 0
                 results.P_AC_dispatch_series(i) = ub(1);
                 results.P_EV_dispatch_series(i) = ub(2);
                 actual_dispatch_fail_case = sum(ub);
                 fprintf('  备用策略：总容量不足，按最大能力调度 %.2f kW (需求 %.2f kW)。\n', actual_dispatch_fail_case, P_need_current_step);
                 results.cost_per_step_series(i) = (c_ac1*ub(1) + c_ac2*ub(1)^2) + (c_ev1*ub(2) + c_ev2*ub(2)^2);
            end
        end
    end
    results.total_event_cost = sum(results.cost_per_step_series(results.status_series >= 0 & ~isinf(results.cost_per_step_series)));
    fprintf('VPP紧急控制指令优化分解完成。总事件有效成本: %.2f 元\n', results.total_event_cost);
end