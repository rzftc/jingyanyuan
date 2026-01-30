function run_scenario_G_comparison1(beta_val, Max_Iter, N_scenarios, N_bus, N_line, dt, ...
    P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
    Reliable_AC_Up, Reliable_EV_Up, Reliable_AC_Down, Reliable_EV_Down, ... 
    Reliable_AC_Base, Reliable_EV_Base, ... 
    R_Gen_Max, R_Shed_Max, cost_params, net_params, direction_signal, options)

    fprintf('\n==========================================================\n');
    fprintf('>>> 场景 G: 确定性 vs 随机优化 效益对比分析 (含分项成本) <<<\n');
    fprintf('==========================================================\n');
    T_steps = length(P_grid_demand);
    
    % [参数设置] 定义统一的违约惩罚价格 (元/MW)，用于统一量纲
    rho_price = 100000; 
    
    %% 1. 运行确定性优化 (Deterministic Dispatch)
    conservative_factor = 1; 
    
    fprintf('  正在运行: 确定性优化 (Conservative Factor = %.2f) ...\n', conservative_factor);
    
    st_det = solve_deterministic_dispatch(P_grid_demand, direction_signal, ...
        Reliable_AC_Up * conservative_factor, ...
        Reliable_EV_Up * conservative_factor, ...
        Reliable_AC_Down * conservative_factor, ...
        Reliable_EV_Down * conservative_factor, ... 
        R_Gen_Max, R_Shed_Max, cost_params, dt, options, T_steps);
    
    % --- 事后评估风险 ---
    % 计算出的 cvar_det_mw 单位是 MW
    [cvar_det_mw, ~] = evaluate_risk_post_hoc(st_det.P_AC, st_det.P_EV, ...
        Scenarios_AC_Up, Scenarios_EV_Up, beta_val, 0.95);
    
    % 将物理风险 (MW) 转换为 风险成本 (元)
    cvar_det_cost = cvar_det_mw * rho_price;
    
    % [新增] 计算分项运行成本
    [cost_det, det_breakdown] = calculate_operating_cost_breakdown(st_det, cost_params, dt);
    
    
    %% 2. 运行随机优化 (Stochastic Dispatch)
    fprintf('  正在运行: 随机优化 (Stochastic, Beta=%d) ...\n', beta_val);
    
    risk_p.beta = beta_val;
    risk_p.confidence = 0.95;
    risk_p.rho_pen = rho_price; % 使用统一定义的惩罚价格
    
    [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast_ramp_tly(...
        P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
        Reliable_AC_Up, Reliable_EV_Up, ... 
        R_Gen_Max, R_Shed_Max, cost_params, risk_p, net_params);
    
    % 处理方向信号
    start_row_net = (T_steps + 1) * N_scenarios;
    for t = 1:T_steps
        if direction_signal(t) == 1
            rows_t = start_row_net + (t-1)*2*N_line + (1 : 2*N_line);
            A(rows_t, info.idx_P_AC(t)) = -A(rows_t, info.idx_P_AC(t));
            A(rows_t, info.idx_P_EV(t)) = -A(rows_t, info.idx_P_EV(t));
            A(rows_t, info.idx_P_Gen(t)) = -A(rows_t, info.idx_P_Gen(t));
        end
    end
    
    % 求解 (使用 Cplex 类对象)
    % -----------------------------------------------------------
    cplex = Cplex('ScenarioG');
    cplex.Model.sense = 'minimize';
    
    cplex.Model.Q = H;
    cplex.Model.obj = f;
    cplex.Model.lb = lb;
    cplex.Model.ub = ub;
    
    % 处理约束
    if isempty(A)
        A_combined = Aeq;
        lhs = beq;
        rhs = beq;
    elseif isempty(Aeq)
        A_combined = A;
        lhs = -inf(size(b));
        rhs = b;
    else
        A_combined = [A; Aeq];
        lhs = [-inf(size(b)); beq];
        rhs = [b; beq];
    end
    cplex.Model.A = A_combined;
    cplex.Model.lhs = lhs;
    cplex.Model.rhs = rhs;
    
    cplex.DisplayFunc = []; 
    
    cplex.solve();
    
    if isfield(cplex.Solution, 'x') && ~isempty(cplex.Solution.x)
        x_opt = cplex.Solution.x;
        status = cplex.Solution.status;
        if status == 1 || status == 101 || status == 102
            exitflag = 1;
        else
            exitflag = -2;
        end
    else
        x_opt = [];
        exitflag = -2;
    end
    % -----------------------------------------------------------
    
    if exitflag > 0
        st_stoch.P_AC = x_opt(info.idx_P_AC);
        st_stoch.P_EV = x_opt(info.idx_P_EV);
        st_stoch.P_Gen = x_opt(info.idx_P_Gen);
        st_stoch.P_Shed = x_opt(info.idx_P_Shed);
        
        eta_val = x_opt(info.idx_eta);
        z_val   = x_opt(info.idx_z);
        % cvar_stoch 的单位已经是 元 (因为 z 和 eta 在约束中已经包含了 rho_price)
        cvar_stoch = eta_val + (1 / (N_scenarios * (1 - risk_p.confidence))) * sum(z_val);
        
        % [新增] 计算分项运行成本
        [cost_stoch, stoch_breakdown] = calculate_operating_cost_breakdown(st_stoch, cost_params, dt);
    else
        warning('随机优化求解失败');
        cost_stoch = NaN; cvar_stoch = NaN;
        stoch_breakdown.base = NaN; stoch_breakdown.res = NaN; stoch_breakdown.shed = NaN;
    end
    
    %% 3. 结果对比与输出
    % 计算总惩罚与风险成本 (切负荷 + CVaR风险)
    total_penalty_det = det_breakdown.shed + cvar_det_cost;
    total_penalty_stoch = stoch_breakdown.shed + cvar_stoch;
    
    fprintf('\n----------------------------------------------------------\n');
    fprintf('指标对比 \t\t| 确定性优化 \t| 随机优化\n');
    fprintf('----------------------------------------------------------\n');
    fprintf('总运行成本 (元) \t| %.2f \t| %.2f \t(+%.2f%%)\n', ...
        cost_det, cost_stoch, (cost_stoch - cost_det)/cost_det*100);
    fprintf('CVaR 风险成本 (元) \t| %.2f \t| %.2f \t(%.2f%%)\n', ...
        cvar_det_cost, cvar_stoch, (cvar_det_cost - cvar_stoch)/cvar_det_cost*100);
    fprintf('----------------------------------------------------------\n');
    fprintf('>>> 分项成本结构对比 <<<\n');
    fprintf('  基础调度成本 (AC+EV) \t| %.2f \t| %.2f\n', det_breakdown.base, stoch_breakdown.base);
    fprintf('  备用/火电成本 (Gen)  \t| %.2f \t| %.2f\n', det_breakdown.res, stoch_breakdown.res);
    fprintf('  惩罚与风险成本       \t| %.2f \t| %.2f\n', total_penalty_det, total_penalty_stoch);
    fprintf('----------------------------------------------------------\n');
    
    %% 4. 绘图
    
    % --- 图 1：成本与风险对比柱状图 ---
    fig = figure('Name', '确定性与随机优化对比', 'Color', 'w', 'Position', [100, 100, 600, 400]);
    
    yyaxis left
    bar([1, 2], [cost_det, cost_stoch], 0.4, 'FaceColor', [0.2 0.6 0.8]);
    ylabel('运行成本 (元)', 'FontSize', 15); 
    ylim([min(cost_det, cost_stoch)*0.8, max(cost_det, cost_stoch)*1.1]);
    
    text(1, cost_det, sprintf('%.0f', cost_det), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14, 'FontWeight', 'bold');
    text(2, cost_stoch, sprintf('%.0f', cost_stoch), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14, 'FontWeight', 'bold');
    
    yyaxis right
    bar([1.4, 2.4], [cvar_det_cost, cvar_stoch], 0.4, 'FaceColor', [0.8 0.3 0.3]);
    ylabel('CVaR 风险成本 (元)', 'FontSize', 15); 
    ylim([0, max(cvar_det_cost, cvar_stoch)*1.5]); 
    
    text(1.4, cvar_det_cost, sprintf('%.2f', cvar_det_cost), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14, 'Color', 'r');
    text(2.4, cvar_stoch, sprintf('%.2f', cvar_stoch), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14, 'Color', 'r');
    
    set(gca, 'XTick', [1.2, 2.2], 'XTickLabel', {'确定性优化', '随机优化'}, 'FontSize', 14); 
    legend({'运行成本', 'CVaR 风险成本'}, 'Location', 'north', 'Orientation', 'horizontal', 'FontSize', 14);
    
    print(fig, '场景G_确定性与随机优化对比.png', '-dpng', '-r600');
    fprintf('  > 已保存对比图: 场景G_确定性与随机优化对比.png\n');
    
    % ================= [图2和图3保持不变] =================
    t_vector = 8 : dt : (8 + (T_steps-1)*dt);
    dir_sign = sign(direction_signal);
    dir_sign(dir_sign == 0) = 1;
    P_Reg_AC_Det = st_det.P_AC .* dir_sign;
    P_Reg_EV_Det = st_det.P_EV .* dir_sign;
    P_Reg_Sum_Det = P_Reg_AC_Det + P_Reg_EV_Det;
    
    P_Reg_AC_Stoch = st_stoch.P_AC .* dir_sign;
    P_Reg_EV_Stoch = st_stoch.P_EV .* dir_sign;
    P_Reg_Sum_Stoch = P_Reg_AC_Stoch + P_Reg_EV_Stoch;
    P_Tot_Sum_Det = (Reliable_AC_Base + Reliable_EV_Base) + P_Reg_Sum_Det;
    P_Tot_Sum_Stoch = (Reliable_AC_Base + Reliable_EV_Base) + P_Reg_Sum_Stoch;
    P_Tot_Base = Reliable_AC_Base + Reliable_EV_Base;
    
    % --- 图 2：AC与EV聚合体总功率对比 ---
    fig2 = figure('Name', 'AC与EV聚合体总功率对比', 'Color', 'w', 'Position', [150, 150, 700, 400]);
    hold on;
    plot(t_vector, P_Tot_Base, 'k--', 'LineWidth', 1.5, 'DisplayName', '基线总功率 (无调节)');
    plot(t_vector, P_Tot_Sum_Det, 'b-', 'LineWidth', 1.5, 'DisplayName', '总功率 (确定性优化)');
    plot(t_vector, P_Tot_Sum_Stoch, 'r-', 'LineWidth', 2.0, 'DisplayName', '总功率 (随机优化)');
    
    ylabel('总功率 (MW)', 'FontName', 'Microsoft YaHei', 'FontSize', 19); 
    xlabel('时刻', 'FontName', 'Microsoft YaHei', 'FontSize', 19);
    xlim([8, 32]);
    set(gca, 'XTick', [8, 12, 16, 20, 24, 28, 32], ...
             'XTickLabel', {'08:00', '12:00', '16:00', '20:00', '00:00', '04:00', '08:00'}, ...
             'FontName', 'Microsoft YaHei', 'FontSize', 16); 
    legend('Location', 'best', 'FontName', 'Microsoft YaHei', 'FontSize', 15);
    print(fig2, '场景G_AC与EV聚合体总功率对比.png', '-dpng', '-r600');
    fprintf('  > 已保存对比图: 场景G_AC与EV聚合体总功率对比.png\n');
    
    % --- 图 3：AC与EV聚合体功率调节量对比 ---
    fig3 = figure('Name', 'AC与EV聚合体功率调节量对比', 'Color', 'w', 'Position', [200, 200, 700, 400]);
    hold on;
    yline(0, 'k-', 'HandleVisibility', 'off');
    
    plot(t_vector, P_Reg_Sum_Det, 'b--', 'LineWidth', 1.5, 'DisplayName', '调节量 (确定性优化)');
    plot(t_vector, P_Reg_Sum_Stoch, 'r-', 'LineWidth', 2.0, 'DisplayName', '调节量 (随机优化)');
    
    x_fill = [t_vector, fliplr(t_vector)];
    y_fill = [P_Reg_Sum_Det', fliplr(P_Reg_Sum_Stoch')];
    fill(x_fill, y_fill, [0.8 0.8 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '策略差异');
    ylabel('功率调节量 (MW)', 'FontName', 'Microsoft YaHei', 'FontSize', 19);
    xlabel('时刻', 'FontName', 'Microsoft YaHei', 'FontSize', 19);
    xlim([8, 32]);
    set(gca, 'XTick', [8, 12, 16, 20, 24, 28, 32], ...
             'XTickLabel', {'08:00', '12:00', '16:00', '20:00', '00:00', '04:00', '08:00'}, ...
             'FontName', 'Microsoft YaHei', 'FontSize', 16); 
    legend('Location', 'best', 'FontName', 'Microsoft YaHei', 'FontSize', 15); 
    print(fig3, '场景G_AC与EV聚合体功率调节量对比.png', '-dpng', '-r600');
    fprintf('  > 已保存对比图: 场景G_AC与EV聚合体功率调节量对比.png\n');

    % ================= [图4：修改了图例大小和位置] =================
    % --- 图 4：成本结构分项对比 (新增) ---
    fig4 = figure('Name', '成本结构分项对比', 'Color', 'w', 'Position', [250, 250, 700, 400]);
    
    % 准备堆积图数据: [基础成本, 备用成本, 惩罚+风险成本]
    y_data = [det_breakdown.base, det_breakdown.res, total_penalty_det; 
              stoch_breakdown.base, stoch_breakdown.res, total_penalty_stoch];
    
    b_stack = bar(y_data, 'stacked');
    
    % 设置颜色 (基础:蓝, 备用:橙, 惩罚:红)
    b_stack(1).FaceColor = [0.2 0.6 0.8];
    b_stack(2).FaceColor = [0.9 0.6 0.2];
    b_stack(3).FaceColor = [0.8 0.3 0.3];
    
    ylabel('成本 (元)', 'FontName', 'Microsoft YaHei', 'FontSize', 19);
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'确定性优化', '随机优化'}, ...
             'FontName', 'Microsoft YaHei', 'FontSize', 16);
    xlim([0.5, 2.5]);
    
    % 设置 Y 轴上限，留出一点空间 (从1.3倍改为1.15倍，因为图例不占顶部了)
    ymax = max(sum(y_data, 2)) * 1.15;
    ylim([0, ymax]);
    
    % [修改点] 缩小字体至12，位置改为 'best' 让其自动寻找合适位置
    legend({'基础调度成本 (AC+EV)', '备用/火电成本', '惩罚与风险成本'}, ...
           'Location', 'best', 'FontName', 'Microsoft YaHei', 'FontSize', 12);
       
    % 添加数值标签
    for i = 1:2
        % 计算每个堆积块的中心位置
        h1 = y_data(i, 1);
        h2 = y_data(i, 2);
        h3 = y_data(i, 3);
        
        % 仅当成本大于0时显示标签，避免拥挤
        if h1 > 100, text(i, h1/2, sprintf('%.0f', h1), 'HorizontalAlignment', 'center', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold'); end
        if h2 > 100, text(i, h1 + h2/2, sprintf('%.0f', h2), 'HorizontalAlignment', 'center', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold'); end
        if h3 > 100, text(i, h1 + h2 + h3/2, sprintf('%.0f', h3), 'HorizontalAlignment', 'center', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold'); end
        
        % 总成本标签
        text(i, h1+h2+h3, sprintf('%.0f', h1+h2+h3), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14, 'FontWeight', 'bold');
    end
    
    print(fig4, '场景G_成本结构分项对比.png', '-dpng', '-r600');
    fprintf('  > 已保存对比图: 场景G_成本结构分项对比.png\n');
end

%% --- 辅助函数：计算分项运行成本 (修改版) ---
function [total_c, breakdown] = calculate_operating_cost_breakdown(st, p, dt)
    % 兼容不同格式的成本系数（单值或向量）
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
    
    % 分项求和
    breakdown.base = sum((term_ac + term_ev) * dt);
    breakdown.res = sum(term_gen * dt);
    breakdown.shed = sum(term_shed * dt);
    
    total_c = breakdown.base + breakdown.res + breakdown.shed;
end

%% --- 辅助函数：事后评估风险 (保持不变) ---
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