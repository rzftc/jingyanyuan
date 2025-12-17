function run_scenario_F_comparison(strategies, P_grid_demand, ...
    Scenarios_AC_Up, Scenarios_EV_Up, Scenarios_AC_Down, Scenarios_EV_Down, ...
    Reliable_AC_Up, Reliable_EV_Up, Reliable_AC_Down, Reliable_EV_Down, ...
    Physical_AC_Up, Physical_EV_Up, Physical_AC_Down, Physical_EV_Down, ...
    R_Gen_Max, R_Shed_Max, cost_params, net_params, direction_signal, dt, options)

    fprintf('\n>>> 场景 F: 不同调度方法的性能对比分析 (Scenario F) <<<\n');
    
    T_steps = length(P_grid_demand);
    [~, N_scenarios] = size(Scenarios_AC_Up);
    N_bus = size(net_params.PTDF, 2);
    N_line = size(net_params.PTDF, 1);
    
    %% 1. 获取三种方法的调度方案
    
    % --- 方法 1: 确定性边界调度 (Deterministic) ---
    fprintf('  - 计算方法 1: 确定性边界调度...\n');
    P_Method1 = solve_deterministic_dispatch(P_grid_demand, direction_signal, ...
        Reliable_AC_Up, Reliable_EV_Up, Reliable_AC_Down, Reliable_EV_Down, ...
        R_Gen_Max, R_Shed_Max, cost_params, dt, options, T_steps);

    % --- 方法 2: 不考虑资源协同 (Uncoordinated, Beta=10, Lambda=0) ---
    fprintf('  - 计算方法 2: 不考虑协同的风险调度 (Beta=10, No SDCI/Rho)...\n');
    risk_p2.beta = 10;
    risk_p2.confidence = 0.95;
    risk_p2.rho_pen = 5000;
    risk_p2.tight_factor = 0.9;
    
    % 构造风险模型但不加协同惩罚
    net_params_safe = net_params;
    net_params_safe.ShedDist = zeros(N_bus, 1);
    
    % 注意：这里需要根据方向选择正确的 Physical Limit，在此简化传入 Up 即可，
    % construct 函数内部逻辑通常处理 capacity，但为保险起见，这里复用 construct_risk_constrained_qp_fast_ramp_tly
    % 需确保传入的 Physical 参数与 Scenario B 一致
    [H2, f2, A2, b2, Aeq2, beq2, lb2, ub2, info2] = construct_risk_constrained_qp_fast_ramp_tly(...
            P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
            Physical_AC_Up, Physical_EV_Up, ...
            R_Gen_Max, R_Shed_Max, ...
            cost_params, risk_p2, net_params_safe);
            
    % 处理方向信号 (Flip Signs)
    start_row_net = 2 * N_scenarios; 
    for t = 1:T_steps
        if direction_signal(t) == 1
            rows_t = start_row_net + (t-1)*2*N_line + (1 : 2*N_line);
            A2(rows_t, info2.idx_P_AC(t)) = -A2(rows_t, info2.idx_P_AC(t));
            A2(rows_t, info2.idx_P_EV(t)) = -A2(rows_t, info2.idx_P_EV(t));
            A2(rows_t, info2.idx_P_Gen(t)) = -A2(rows_t, info2.idx_P_Gen(t));
        end
    end
    
    % 缩放 dt
    idx_pow2 = [info2.idx_P_AC, info2.idx_P_EV, info2.idx_P_Gen, info2.idx_P_Shed];
    for idx = idx_pow2, H2(idx, idx) = H2(idx, idx) * dt; end
    f2(idx_pow2) = f2(idx_pow2) * dt;
    
    % 求解 (仅一次，不进行 SDCI 迭代)
    [x_opt2, ~, flag2] = quadprog(H2, f2, A2, b2, Aeq2, beq2, lb2, ub2, [], options);
    if flag2 > 0
        P_Method2.P_AC = x_opt2(info2.idx_P_AC);
        P_Method2.P_EV = x_opt2(info2.idx_P_EV);
        P_Method2.P_Gen = x_opt2(info2.idx_P_Gen);
        P_Method2.P_Shed = x_opt2(info2.idx_P_Shed);
    else
        warning('方法 2 求解失败，使用零值代替。');
        P_Method2.P_AC = zeros(T_steps,1); P_Method2.P_EV = zeros(T_steps,1);
        P_Method2.P_Gen = zeros(T_steps,1); P_Method2.P_Shed = zeros(T_steps,1);
    end

    % --- 方法 3: 提议方法 (from Scenario B, Beta=10) ---
    fprintf('  - 获取方法 3: 提议方法 (Risk + Coordination)...\n');
    % 假设 strategies{3} 对应 Beta=10 (根据 run_scenario_B_tly 的 beta_values=[0,1,10])
    if length(strategies) >= 3 && ~isempty(strategies{3})
        P_Method3 = strategies{3};
    else
        error('未找到方法 3 的数据，请确保先运行了场景 B 且 beta_values 包含 10。');
    end

    %% 2. 性能评估 (基于 1000 个真实场景的回测)
    methods_list = {P_Method1, P_Method2, P_Method3};
    method_names = {'确定性边界', '无协同调度', '提议方法'};
    
    results = struct('BaseCost', zeros(1,3), 'PenaltyCost', zeros(1,3), 'TotalCost', zeros(1,3));
    
    % 设定违约惩罚价格 (通常比调度成本高，用于评估安全性)
    Penalty_Price = 2000; % 元/MW，用于评估时的虚拟惩罚
    
    for m = 1:3
        st = methods_list{m};
        
        % 1. 基础调度成本 (Operational Cost)
        cost_base = sum((cost_params.c1_ac * st.P_AC + cost_params.c2_ac * st.P_AC.^2)*dt + ...
                        (cost_params.c1_ev * st.P_EV + cost_params.c2_ev * st.P_EV.^2)*dt + ...
                        (cost_params.c1_gen * st.P_Gen + cost_params.c2_gen * st.P_Gen.^2)*dt + ...
                        (cost_params.c1_shed * st.P_Shed)*dt);
        
        % 2. 违约风险评估 (Validation against Scenarios)
        total_violation = 0;
        
        for s = 1:N_scenarios
            % 构建该场景的真实能力
            Cap_AC_s = zeros(T_steps, 1);
            Cap_EV_s = zeros(T_steps, 1);
            for t = 1:T_steps
                if direction_signal(t) == 1
                    Cap_AC_s(t) = Scenarios_AC_Up(t, s);
                    Cap_EV_s(t) = Scenarios_EV_Up(t, s);
                else
                    Cap_AC_s(t) = abs(Scenarios_AC_Down(t, s));
                    Cap_EV_s(t) = abs(Scenarios_EV_Down(t, s));
                end
            end
            
            % 计算违约量
            viol_ac = max(0, st.P_AC - Cap_AC_s);
            viol_ev = max(0, st.P_EV - Cap_EV_s);
            total_violation = total_violation + sum(viol_ac + viol_ev);
        end
        
        avg_violation = total_violation / N_scenarios;
        cost_penalty = avg_violation * Penalty_Price * dt; % 期望违约惩罚
        
        results.BaseCost(m) = cost_base;
        results.PenaltyCost(m) = cost_penalty;
        results.TotalCost(m) = cost_base + cost_penalty;
        
        fprintf('    [%s] 基础成本: %.2f, 期望违约惩罚: %.2f, 综合性能: %.2f\n', ...
            method_names{m}, cost_base, cost_penalty, results.TotalCost(m));
    end
    
    %% 3. 绘图对比
    fig_comp = figure('Name', '三种调度方法性能对比', 'Color', 'w', 'Position', [400, 300, 700, 500]);
    
    y_data = [results.BaseCost; results.PenaltyCost]'; % 堆叠数据
    b = bar(y_data, 'stacked');
    
    b(1).FaceColor = [0.2 0.6 0.8]; % 基础成本颜色
    b(2).FaceColor = [0.8 0.4 0.4]; % 违约风险颜色
    
    set(gca, 'XTickLabel', method_names, 'FontSize', 12, 'FontName', 'Microsoft YaHei');
    ylabel('综合期望成本 (元)', 'FontSize', 12, 'FontName', 'Microsoft YaHei');
    legend({'调度运行成本', '潜在违约风险惩罚'}, 'Location', 'northwest');
    grid on;
    
    % 在柱顶标注总成本
    xtips = 1:3;
    ytips = results.TotalCost;
    labels = string(round(ytips, 0));
    text(xtips, ytips, labels, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 11, 'FontWeight', 'bold');
    
    print(fig_comp, '三种调度方法性能对比.png', '-dpng', '-r300');
    fprintf('  > 绘图已保存: 三种调度方法性能对比.png\n');
end