function run_scenario_G_comparison(beta_val, Max_Iter, N_scenarios, N_bus, N_line, dt, ...
    P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
    Reliable_AC_Up, Reliable_EV_Up, Reliable_AC_Down, Reliable_EV_Down, ... 
    R_Gen_Max, R_Shed_Max, cost_params, net_params, direction_signal, options)

    fprintf('\n==========================================================\n');
    fprintf('>>> 场景 G: 确定性 vs 随机优化 效益对比分析 <<<\n');
    fprintf('==========================================================\n');
    T_steps = length(P_grid_demand);
    
    %% 1. 运行确定性优化 (Deterministic Dispatch)
    % 修改策略：引入保守系数 (Conservative Factor)
    % 原因：确定性优化过于乐观，倾向于用满所有廉价VPP资源。
    %      通过乘以 0.6~0.7 的系数，模拟调度员预留裕度，强制模型动用一部分火电，从而提高成本。
    conservative_factor = 1;  % <--- [关键修改] 这里可以调节 (0.6 表示只敢用 60% 的可靠容量)
    
    fprintf('  正在运行: 确定性优化 (Conservative Factor = %.2f) ...\n', conservative_factor);
    
    % 注意：传入的 Reliable_AC/EV 边界被乘上了 conservative_factor
    st_det = solve_deterministic_dispatch(P_grid_demand, direction_signal, ...
        Reliable_AC_Up * conservative_factor, ...
        Reliable_EV_Up * conservative_factor, ...
        Reliable_AC_Down * conservative_factor, ...
        Reliable_EV_Down * conservative_factor, ... 
        R_Gen_Max, R_Shed_Max, cost_params, dt, options, T_steps);
    
    % --- 事后评估风险 ---
    % 注意：评估风险时，还是用原始场景去检验，看看这个保守策略到底安不安全
    [cvar_det, ~] = evaluate_risk_post_hoc(st_det.P_AC, st_det.P_EV, ...
        Scenarios_AC_Up, Scenarios_EV_Up, beta_val, 0.95);
    
    % 计算运行成本
    cost_det = calculate_operating_cost(st_det, cost_params, dt);
    
    
    %% 2. 运行随机优化 (Stochastic Dispatch)
    % 随机优化不需要 derating，因为它内部通过 CVaR 自动决定用多少
    fprintf('  正在运行: 随机优化 (Stochastic, Beta=%d) ...\n', beta_val);
    
    risk_p.beta = beta_val;
    risk_p.confidence = 0.95;
    risk_p.rho_pen = 5000; 
    
    [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast_ramp_tly(...
        P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
        Reliable_AC_Up, Reliable_EV_Up, ... 
        R_Gen_Max, R_Shed_Max, cost_params, risk_p, net_params);
    
    % 处理方向信号
    start_row_net = 2 * N_scenarios;
    for t = 1:T_steps
        if direction_signal(t) == 1
            rows_t = start_row_net + (t-1)*2*N_line + (1 : 2*N_line);
            A(rows_t, info.idx_P_AC(t)) = -A(rows_t, info.idx_P_AC(t));
            A(rows_t, info.idx_P_EV(t)) = -A(rows_t, info.idx_P_EV(t));
            A(rows_t, info.idx_P_Gen(t)) = -A(rows_t, info.idx_P_Gen(t));
        end
    end
    
    % 求解
    [x_opt, ~, exitflag] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
    
    if exitflag > 0
        st_stoch.P_AC = x_opt(info.idx_P_AC);
        st_stoch.P_EV = x_opt(info.idx_P_EV);
        st_stoch.P_Gen = x_opt(info.idx_P_Gen);
        st_stoch.P_Shed = x_opt(info.idx_P_Shed);
        
        eta_val = x_opt(info.idx_eta);
        z_val   = x_opt(info.idx_z);
        cvar_stoch = eta_val + (1 / (N_scenarios * (1 - risk_p.confidence))) * sum(z_val);
        
        cost_stoch = calculate_operating_cost(st_stoch, cost_params, dt);
    else
        warning('随机优化求解失败');
        cost_stoch = NaN; cvar_stoch = NaN;
    end
    
    %% 3. 结果对比与输出
    fprintf('\n----------------------------------------------------------\n');
    fprintf('指标对比 \t\t| 确定性优化 \t| 随机优化\n');
    fprintf('----------------------------------------------------------\n');
    fprintf('运行成本 (元) \t| %.2f \t| %.2f \t(+%.2f%%)\n', ...
        cost_det, cost_stoch, (cost_stoch - cost_det)/cost_det*100);
    fprintf('CVaR 风险 (MW) \t| %.4f \t| %.4f \t(%.2f%%)\n', ...
        cvar_det, cvar_stoch, (cvar_det - cvar_stoch)/cvar_det*100);
    fprintf('----------------------------------------------------------\n');
    
    % --- 绘图逻辑 ---
    fig = figure('Name', '确定性与随机优化对比', 'Color', 'w', 'Position', [100, 100, 600, 400]);
    
    yyaxis left
    b1 = bar([1, 2], [cost_det, cost_stoch], 0.4, 'FaceColor', [0.2 0.6 0.8]);
    ylabel('运行成本 (元)', 'FontSize', 11);
    ylim([min(cost_det, cost_stoch)*0.8, max(cost_det, cost_stoch)*1.1]);
    
    text(1, cost_det, sprintf('%.0f', cost_det), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'FontWeight', 'bold');
    text(2, cost_stoch, sprintf('%.0f', cost_stoch), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'FontWeight', 'bold');
    
    yyaxis right
    b2 = bar([1.4, 2.4], [cvar_det, cvar_stoch], 0.4, 'FaceColor', [0.8 0.3 0.3]);
    ylabel('CVaR 风险 (MW)', 'FontSize', 11);
    ylim([0, max(cvar_det, cvar_stoch)*1.5]); 
    
    text(1.4, cvar_det, sprintf('%.4f', cvar_det), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'r');
    text(2.4, cvar_stoch, sprintf('%.4f', cvar_stoch), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'r');
    
    set(gca, 'XTick', [1.2, 2.2], 'XTickLabel', {'确定性优化', '随机优化'}, 'FontSize', 11);
    legend({'运行成本', 'CVaR 风险'}, 'Location', 'north', 'Orientation', 'horizontal');
    grid on;
    
    print(fig, '场景G_确定性与随机优化对比.png', '-dpng', '-r600');
    fprintf('  > 已保存对比图: 场景G_确定性与随机优化对比.png\n');
end

%% --- 辅助函数：计算运行成本 (修复向量乘法) ---
function c = calculate_operating_cost(st, p, dt)
    % 使用 .* 确保分时电价向量与功率向量正确相乘
    if length(p.c1_ac) > 1
        term_ac = (p.c1_ac .* st.P_AC + p.c2_ac .* st.P_AC.^2);
        term_ev = (p.c1_ev .* st.P_EV + p.c2_ev .* st.P_EV.^2);
        term_gen = (p.c1_gen .* st.P_Gen + p.c2_gen .* st.P_Gen.^2);
    else
        term_ac = (p.c1_ac * st.P_AC + p.c2_ac * st.P_AC.^2);
        term_ev = (p.c1_ev * st.P_EV + p.c2_ev * st.P_EV.^2);
        term_gen = (p.c1_gen * st.P_Gen + p.c2_gen * st.P_Gen.^2);
    end
    term_shed = (p.c1_shed * st.P_Shed);
    c = sum((term_ac + term_ev + term_gen + term_shed) * dt);
end

%% --- 辅助函数：事后评估风险 ---
function [cvar, risk_val] = evaluate_risk_post_hoc(P_AC, P_EV, S_AC, S_EV, beta, conf)
    [~, N_scenarios] = size(S_AC);
    total_sched = P_AC + P_EV;
    total_cap_s = S_AC + S_EV;
    
    violations = zeros(N_scenarios, 1);
    for s = 1:N_scenarios
        violations(s) = max(0, max(total_sched - total_cap_s(:, s))); 
    end
    
    sorted_v = sort(violations);
    idx_var = ceil(conf * N_scenarios);
    if idx_var > N_scenarios, idx_var = N_scenarios; end
    tail_losses = sorted_v(idx_var:end);
    
    if isempty(tail_losses)
        cvar = 0;
    else
        cvar = mean(tail_losses);
    end
    risk_val = cvar;
end