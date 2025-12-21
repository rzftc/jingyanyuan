function run_scenario_G_comparison(beta_val, Max_Iter, N_scenarios, N_bus, N_line, dt, ...
    P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
    Reliable_AC_Up, Reliable_EV_Up, ... 
    R_Gen_Max, R_Shed_Max, cost_params, net_params, direction_signal, options)
    fprintf('\n==========================================================\n');
    fprintf('>>> 场景 G: 确定性 vs 随机优化 效益对比分析 <<<\n');
    fprintf('==========================================================\n');
    T_steps = length(P_grid_demand);
    
    %% 1. 运行确定性优化 (Deterministic Dispatch)
    % 逻辑：只把“可靠调节域”当作唯一的真理，完全不考虑具体的场景分布
    fprintf('  正在运行: 确定性优化 (Deterministic) ...\n');
    st_det = solve_deterministic_dispatch(P_grid_demand, direction_signal, ...
        Reliable_AC_Up, Reliable_EV_Up, Reliable_AC_Up, Reliable_EV_Up, ... % 上下限都传可靠边界
        R_Gen_Max, R_Shed_Max, cost_params, dt, options, T_steps);
    
    % --- 关键：事后评估确定性策略的“真实风险” ---
    % 虽然确定性模型自己觉得风险是0，但我们用蒙特卡洛场景去检验它，会发现它其实有风险
    [cvar_det, risk_cost_det] = evaluate_risk_post_hoc(st_det.P_AC, st_det.P_EV, ...
        Scenarios_AC_Up, Scenarios_EV_Up, beta_val, 0.95);
    
    % 计算确定性策略的运行成本
    cost_det = calculate_operating_cost(st_det, cost_params, dt);
    
    
    %% 2. 运行随机优化 (Stochastic Dispatch)
    % 逻辑：考虑所有场景的分布，并显式优化 CVaR
    fprintf('  正在运行: 随机优化 (Stochastic, Beta=%d) ...\n', beta_val);
    
    risk_p.beta = beta_val;
    risk_p.confidence = 0.95;
    risk_p.rho_pen = 5000; 
    
    % 复用场景 B/F 的逻辑，不需要 SDCI/Rho (设为0)，纯粹对比随机性
    % 注意：为了公平对比，这里也把上限设为 Reliable_AC_Up (稳健策略)，看它是否会进一步回撤
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
        
        % 提取 CVaR
        eta_val = x_opt(info.idx_eta);
        z_val   = x_opt(info.idx_z);
        cvar_stoch = eta_val + (1 / (N_scenarios * (1 - risk_p.confidence))) * sum(z_val);
        risk_cost_stoch = cvar_stoch; % 这里直接用 CVaR 值作为风险度量
        
        cost_stoch = calculate_operating_cost(st_stoch, cost_params, dt);
    else
        warning('随机优化求解失败');
        cost_stoch = NaN; cvar_stoch = NaN;
    end
    
    %% 3. 结果对比与输出
    fprintf('\n----------------------------------------------------------\n');
    fprintf('指标对比 \t\t| 确定性优化 (Deterministic) \t| 随机优化 (Stochastic)\n');
    fprintf('----------------------------------------------------------\n');
    fprintf('运行成本 (元) \t| %.2f \t\t\t| %.2f \t(+%.2f%%)\n', ...
        cost_det, cost_stoch, (cost_stoch - cost_det)/cost_det*100);
    fprintf('CVaR 风险 (MW) \t| %.4f \t\t\t| %.4f \t(-%.2f%%)\n', ...
        cvar_det, cvar_stoch, (cvar_det - cvar_stoch)/cvar_det*100);
    fprintf('----------------------------------------------------------\n');
    
   %% 4. 绘图：双柱状图对比 (带数值标注版)
    fig = figure('Name', '确定性与随机优化对比', 'Color', 'w', 'Position', [100, 100, 600, 400]);
    
    % --- 左轴：绘制运行成本 ---
    yyaxis left
    % 绘制柱状图 (位置在 x=1 和 x=2)
    b1 = bar([1, 2], [cost_det, cost_stoch], 0.4, 'FaceColor', [0.2 0.6 0.8]);
    ylabel('运行成本 (元)', 'FontSize', 11);
    % 设置左轴范围，留出顶部空间给文字
    ylim([min(cost_det, cost_stoch)*0.9, max(cost_det, cost_stoch)*1.1]);
    
    % [新增] 添加成本数值标注 (显示整数)
    text(1, cost_det, sprintf('%.0f', cost_det), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
    text(2, cost_stoch, sprintf('%.0f', cost_stoch), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
    
    % --- 右轴：绘制风险 ---
    yyaxis right
    % 绘制柱状图 (位置错开，在 x=1.4 和 x=2.4)
    b2 = bar([1.4, 2.4], [cvar_det, cvar_stoch], 0.4, 'FaceColor', [0.8 0.3 0.3]);
    ylabel('CVaR 风险 (MW)', 'FontSize', 11);
    % 设置右轴范围
    ylim([0, max(cvar_det, cvar_stoch)*1.2]); 
    
    % [新增] 添加风险数值标注 (保留4位小数)
    % 注意：必须在 yyaxis right 激活状态下绘制，坐标才对
    text(1.4, cvar_det, sprintf('%.4f', cvar_det), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold');
    text(2.4, cvar_stoch, sprintf('%.4f', cvar_stoch), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold');
    
    % --- 公共设置 ---
    set(gca, 'XTick', [1.2, 2.2], 'XTickLabel', {'确定性优化', '随机优化'}, 'FontSize', 11);
    
    % 图例设置 (中文，上方居中)
    legend({'运行成本 (左轴)', 'CVaR 风险 (右轴)'}, ...
        'Location', 'north', 'Orientation', 'horizontal', 'Box', 'off');
    
    grid on; box on;
    
    % 增加数值标签
    text(1, cost_det, sprintf('%.0f', cost_det), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(2, cost_stoch, sprintf('%.0f', cost_stoch), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    
    % 保存图片（中文文件名）
    print(fig, '场景G_确定性与随机优化对比.png', '-dpng', '-r600');
    print(fig, '场景G_确定性与随机优化对比.emf', '-dmeta');
    fprintf('  > 已保存对比图: 场景G_确定性与随机优化对比.png\n');
end
%% --- 辅助函数：计算运行成本 ---
function c = calculate_operating_cost(st, p, dt)
    term_ac = (p.c1_ac .* st.P_AC + p.c2_ac .* st.P_AC.^2);
    term_ev = (p.c1_ev .* st.P_EV + p.c2_ev .* st.P_EV.^2);
    term_gen = (p.c1_gen .* st.P_Gen + p.c2_gen .* st.P_Gen.^2);
    term_shed = (p.c1_shed .* st.P_Shed);
    c = sum((term_ac + term_ev + term_gen + term_shed) * dt);
end
%% --- 辅助函数：事后评估风险 (核心) ---
function [cvar, risk_val] = evaluate_risk_post_hoc(P_AC, P_EV, S_AC, S_EV, beta, conf)
    % 这是一个"照妖镜"函数
    % 确定性策略以为自己没风险，我们用几千个场景去照它，看看它到底违约了多少
    
    [T, N_scenarios] = size(S_AC);
    total_sched = P_AC + P_EV;      % T x 1
    total_cap_s = S_AC + S_EV;      % T x N
    
    % 计算每个场景、每个时刻的违约量
    % Violation = max(0, 调度计划 - 实际能力)
    violations = zeros(N_scenarios, 1);
    
    for s = 1:N_scenarios
        % 这里采用"最严重时刻"逻辑，与随机优化中的约束一致
        % 即：只要这一天里有一个时刻违约了，这一天的场景就算违约
        v_t = max(0, total_sched - total_cap_s(:, s));
        violations(s) = max(v_t); 
    end
    
    % 计算 CVaR (VaR + 尾部均值)
    % 1. 排序违约量
    sorted_v = sort(violations);
    % 2. 找到 VaR (95% 分位数)
    idx_var = ceil(conf * N_scenarios);
    VaR = sorted_v(idx_var);
    % 3. 尾部均值 (超过 VaR 的部分的平均)
    % 注意：标准的 CVaR 定义
    tail_losses = sorted_v(idx_var:end);
    cvar = mean(tail_losses);
    
    risk_val = cvar;
end