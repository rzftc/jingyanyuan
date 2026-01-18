function run_scenario_H_robust_comparison(beta_val, N_scenarios, N_bus, N_line, dt, ...
    P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
    R_Gen_Max, R_Shed_Max, cost_params, net_params, direction_signal, options)

    fprintf('\n==========================================================\n');
    fprintf('>>> 场景 H: 传统鲁棒优化(Robust)与CVaR随机优化对比验证 (修正版) <<<\n');
    fprintf('==========================================================\n');
    
    T_steps = length(P_grid_demand);
    
    %% 1. 运行鲁棒优化 (Robust Optimization - Worst Case Benchmark)
    % 修正说明：为了公平对比，鲁棒优化必须使用与CVaR相同的物理模型（含PTDF、爬坡、状态方程），
    % 区别在于鲁棒优化将物理边界（ub）强制限制在“最恶劣场景（最小值）”，且不考虑风险惩罚项。
    
    fprintf('  正在运行: 传统鲁棒优化 (基于最恶劣场景边界 + 全物理约束) ...\n');
    
    % 1.1 构造鲁棒参数
    rob_risk_p.beta = 0;        % 鲁棒优化不含CVaR目标惩罚
    rob_risk_p.confidence = 0.99; 
    rob_risk_p.rho_pen = 0;     % 不需要的惩罚设为0
    
    % 1.2 计算鲁棒边界：取所有场景在各时刻的最小值 (Box Uncertainty Set)
    Robust_AC_Bound = min(Scenarios_AC_Up, [], 2); 
    Robust_EV_Bound = min(Scenarios_EV_Up, [], 2);
    
    % 1.3 构建 QP 问题 (使用全约束模型)
    % 关键：将 R_AC/R_EV (即物理ub) 替换为 Robust_AC_Bound
    [H_r, f_r, A_r, b_r, Aeq_r, beq_r, lb_r, ub_r, info_r] = construct_risk_constrained_qp_fast_ramp_tly(...
        P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
        Robust_AC_Bound, Robust_EV_Bound, ... % <--- 强制约束为最小值
        R_Gen_Max, R_Shed_Max, cost_params, rob_risk_p, net_params);
    
    % 1.4 处理方向信号 (与主流程一致)
    % start_row_net = 2 * N_scenarios;
    start_row_net = (T_steps + 1) * N_scenarios;
    for t = 1:T_steps
        if direction_signal(t) == 1
            rows_t = start_row_net + (t-1)*2*N_line + (1 : 2*N_line);
            A_r(rows_t, info_r.idx_P_AC(t)) = -A_r(rows_t, info_r.idx_P_AC(t));
            A_r(rows_t, info_r.idx_P_EV(t)) = -A_r(rows_t, info_r.idx_P_EV(t));
            A_r(rows_t, info_r.idx_P_Gen(t)) = -A_r(rows_t, info_r.idx_P_Gen(t));
        end
    end
    
    % 1.5 求解鲁棒模型
    % [x_rob, ~, exitflag_rob] = quadprog(H_r, f_r, A_r, b_r, Aeq_r, beq_r, lb_r, ub_r, [], options);
    [x_rob, ~, exitflag_rob] = cplexqp(H_r, f_r, A_r, b_r, Aeq_r, beq_r, lb_r, ub_r, [], options);
    if exitflag_rob > 0
        st_rob.P_AC = x_rob(info_r.idx_P_AC);
        st_rob.P_EV = x_rob(info_r.idx_P_EV);
        st_rob.P_Gen = x_rob(info_r.idx_P_Gen);
        st_rob.P_Shed = x_rob(info_r.idx_P_Shed);
        
        cost_rob = calculate_operating_cost(st_rob, cost_params, dt);
        [~, viol_rob] = evaluate_risk_post_hoc(st_rob.P_AC, st_rob.P_EV, Scenarios_AC_Up, Scenarios_EV_Up);
    else
        warning('鲁棒优化求解失败 (可能是最恶劣场景下无解)');
        cost_rob = NaN; viol_rob = NaN; st_rob.P_AC = zeros(T_steps,1); st_rob.P_EV = zeros(T_steps,1);
    end

    %% 2. 运行 CVaR 随机优化 (Stochastic Optimization)
    fprintf('  正在运行: CVaR 随机优化 (Beta=%d, 全物理约束) ...\n', beta_val);
    
    % 2.1 构造CVaR参数
    cvar_risk_p.beta = beta_val;
    cvar_risk_p.confidence = 0.95;
    cvar_risk_p.rho_pen = 5000; 
    
    % 2.2 物理边界使用最大值 (允许在风险约束下探索更大空间)
    Bound_AC_Ref = max(Scenarios_AC_Up, [], 2); 
    Bound_EV_Ref = max(Scenarios_EV_Up, [], 2);

    % 2.3 构建 QP 问题
    [H_c, f_c, A_c, b_c, Aeq_c, beq_c, lb_c, ub_c, info_c] = construct_risk_constrained_qp_fast_ramp_tly(...
        P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
        Bound_AC_Ref, Bound_EV_Ref, ... 
        R_Gen_Max, R_Shed_Max, cost_params, cvar_risk_p, net_params);
    
    % 2.4 处理方向信号
    for t = 1:T_steps
        if direction_signal(t) == 1
            rows_t = start_row_net + (t-1)*2*N_line + (1 : 2*N_line);
            A_c(rows_t, info_c.idx_P_AC(t)) = -A_c(rows_t, info_c.idx_P_AC(t));
            A_c(rows_t, info_c.idx_P_EV(t)) = -A_c(rows_t, info_c.idx_P_EV(t));
            A_c(rows_t, info_c.idx_P_Gen(t)) = -A_c(rows_t, info_c.idx_P_Gen(t));
        end
    end
    
    % 2.5 求解CVaR模型
    [x_cvar, ~, exitflag_cvar] = quadprog(H_c, f_c, A_c, b_c, Aeq_c, beq_c, lb_c, ub_c, [], options);
    
    if exitflag_cvar > 0
        st_cvar.P_AC = x_cvar(info_c.idx_P_AC);
        st_cvar.P_EV = x_cvar(info_c.idx_P_EV);
        st_cvar.P_Gen = x_cvar(info_c.idx_P_Gen);
        st_cvar.P_Shed = x_cvar(info_c.idx_P_Shed);
        
        cost_cvar = calculate_operating_cost(st_cvar, cost_params, dt);
        [~, viol_cvar] = evaluate_risk_post_hoc(st_cvar.P_AC, st_cvar.P_EV, Scenarios_AC_Up, Scenarios_EV_Up);
    else
        warning('CVaR 优化求解失败');
        cost_cvar = NaN; viol_cvar = NaN;
    end
    
    %% 3. 结果数据对比与打印
    fprintf('\n----------------------------------------------------------\n');
    fprintf('对比指标 (场景 H 修正版) \t| 鲁棒优化 (Robust) \t| CVaR 随机优化\n');
    fprintf('----------------------------------------------------------\n');
    fprintf('运行成本 (元) \t\t| %.2f \t\t| %.2f \n', cost_rob, cost_cvar);
    
    cost_saving = (cost_rob - cost_cvar)/cost_rob*100;
    fprintf('成本节省率 (%%) \t\t| -- \t\t\t| %.2f%% \n', cost_saving);
    fprintf('平均潜在违约 (MW)\t\t| %.4f \t\t| %.4f \n', viol_rob, viol_cvar);
    fprintf('----------------------------------------------------------\n');
    
    %% 4. 绘图 (无标题，中文保存)
    
    % --- 图 H1：成本对比柱状图 ---
    fig_h1 = figure('Name', 'Robust_vs_CVaR_Cost', 'Color', 'w', 'Position', [100, 100, 500, 400]);
    b = bar([cost_rob, cost_cvar], 0.5);
    b.FaceColor = 'flat';
    b.CData(1,:) = [0.6 0.6 0.6]; % 灰色代表保守的鲁棒
    b.CData(2,:) = [0.2 0.6 0.8]; % 蓝色代表本文方法
    ylabel('系统运行成本 (元)', 'FontName', 'Microsoft YaHei', 'FontSize', 12);
    set(gca, 'XTickLabel', {'鲁棒优化', 'CVaR调度'}, 'FontName', 'Microsoft YaHei', 'FontSize', 12);
    grid on;
    % 添加数值标签
    if ~isnan(cost_rob), text(1, cost_rob, sprintf('%.0f', cost_rob), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 11); end
    if ~isnan(cost_cvar), text(2, cost_cvar, sprintf('%.0f', cost_cvar), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 11); end
    
    print(fig_h1, '场景H_鲁棒与CVaR成本对比.png', '-dpng', '-r600');
    
    % --- 图 H2：调度策略与边界对比 (修正版：真实正负调度量 + 08:00起止) ---
    if ~isnan(cost_rob) && ~isnan(cost_cvar)
        fig_h2 = figure('Name', 'Robust_vs_CVaR_Dispatch', 'Color', 'w', 'Position', [650, 100, 700, 400]);
        hold on;
        
        % 1. 构建时间轴 (08:00 -> 次日 08:00)
        % 假设 dt 是以小时为单位 (例如 0.25h), 起始点是 8
        t_vector = 8 : dt : (8 + (T_steps-1)*dt);
        
        % 2. 处理正负方向
        % 原始 P_AC, P_EV 均为非负的“调节幅度”
        % 真实调度量 = 调节幅度 * 方向信号 (+1/-1)
        dir_sign = sign(direction_signal);
        dir_sign(dir_sign == 0) = 1; % 默认为正

        P_Real_Rob  = (st_rob.P_AC + st_rob.P_EV) .* dir_sign;
        P_Real_CVaR = (st_cvar.P_AC + st_cvar.P_EV) .* dir_sign;
        
        % 鲁棒边界也需要根据方向翻转，展示“调节包络”
        Bound_Real = (Robust_AC_Bound + Robust_EV_Bound) .* dir_sign;

        % 3. 绘图
        % 绘制最恶劣边界 (鲁棒基准)
        h1 = plot(t_vector, Bound_Real, 'r:', 'LineWidth', 1.5, 'DisplayName', '最恶劣场景边界');

        % 绘制调度曲线
        h2 = plot(t_vector, P_Real_Rob, 'k--', 'LineWidth', 1.5, 'DisplayName', '鲁棒优化调度');
        h3 = plot(t_vector, P_Real_CVaR, 'b-', 'LineWidth', 2, 'DisplayName', 'CVaR随机优化调度');

        % 填充差异区域
        x_fill = [t_vector, fliplr(t_vector)];
        y_fill = [P_Real_Rob', fliplr(P_Real_CVaR')];
        fill(x_fill, y_fill, [0.2 0.6 0.8], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', '释放的调节潜力');

        ylabel('VPP 实际调度功率 (MW)', 'FontName', 'Microsoft YaHei', 'FontSize', 12);
        xlabel('时刻', 'FontName', 'Microsoft YaHei', 'FontSize', 12);
        
        % 设置 X 轴刻度 (8:00 到 次日 8:00)
        xlim([8, 32]);
        set(gca, 'XTick', [8, 12, 16, 20, 24, 28, 32]);
        set(gca, 'XTickLabel', {'08:00', '12:00', '16:00', '20:00', '00:00', '04:00', '08:00'}, 'FontSize', 11);
        
        legend([h1, h2, h3], 'Location', 'best', 'FontName', 'Microsoft YaHei');
        grid on;
        
        print(fig_h2, '场景H_调度策略与边界对比.png', '-dpng', '-r600');
    end
    
    fprintf('  > 已保存对比图: 场景H_鲁棒与CVaR成本对比.png, 场景H_调度策略与边界对比.png\n');
end

%% --- 辅助函数：计算运行成本 ---
function c = calculate_operating_cost(st, p, dt)
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

%% --- 辅助函数：事后评估违约 ---
function [cvar, risk_val] = evaluate_risk_post_hoc(P_AC, P_EV, S_AC, S_EV)
    [~, N_scenarios] = size(S_AC);
    total_sched = P_AC + P_EV;
    total_cap_s = S_AC + S_EV; 
    
    violations = zeros(N_scenarios, 1);
    for s = 1:N_scenarios
        violations(s) = sum(max(0, total_sched - total_cap_s(:, s))); 
    end
    
    risk_val = mean(violations); % 平均违约量
    cvar = max(violations);      % 最大违约量
end