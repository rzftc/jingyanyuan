function run_test3_risk_analysis(t_axis, P_req_pu, dir_vec, net, cost_p, ...
    S_AC_U, S_EV_U, S_AC_D, S_EV_D, R_AC_U, R_EV_U, R_AC_D, R_EV_D, baseMVA)

    fprintf('  - 正在对比风险偏好...\n');
    T = length(t_axis);
    
    % 准备数据: 根据 dir_vec 组装当前时刻的场景
    S_AC_curr = zeros(T, size(S_AC_U, 2));
    S_EV_curr = zeros(T, size(S_EV_U, 2));
    R_AC_curr = zeros(T, 1);
    R_EV_curr = zeros(T, 1);
    
    for t = 1:T
        if dir_vec(t) == 1
            S_AC_curr(t,:) = S_AC_U(t,:); S_EV_curr(t,:) = S_EV_U(t,:);
            R_AC_curr(t) = R_AC_U(t);     R_EV_curr(t) = R_EV_U(t);
        else
            S_AC_curr(t,:) = S_AC_D(t,:); S_EV_curr(t,:) = S_EV_D(t,:);
            R_AC_curr(t) = R_AC_D(t);     R_EV_curr(t) = R_EV_D(t);
        end
    end
    
    % 策略A: 风险中性 (Beta=0)
    risk_A.beta = 0; risk_A.confidence = 0.95; risk_A.rho_pen = 1000;
    [P_AC_A, P_EV_A, P_sl_A] = solve_dist_optimization_corrected(P_req_pu, net, cost_p, risk_A, S_AC_curr, S_EV_curr, R_AC_curr, R_EV_curr, 0, 0, dir_vec);
    
    % 策略B: 风险规避 (Beta=10)
    risk_B.beta = 10; risk_B.confidence = 0.95; risk_B.rho_pen = 1000;
    [P_AC_B, P_EV_B, P_sl_B] = solve_dist_optimization_corrected(P_req_pu, net, cost_p, risk_B, S_AC_curr, S_EV_curr, R_AC_curr, R_EV_curr, 0, 0, dir_vec);
    
    % 绘图
    fig = figure('Name', 'Test3_Risk', 'Color', 'w', 'Position', [100, 100, 800, 600]);
    subplot(2,1,1); hold on;
    % 绘制带符号的总出力 (Up为正, Down为负以区分显示)
    P_total_A = (P_AC_A + P_EV_A) .* dir_vec * baseMVA;
    P_total_B = (P_AC_B + P_EV_B) .* dir_vec * baseMVA;
    P_target  = P_req_pu .* dir_vec * baseMVA;
    
    plot(t_axis, P_target, 'k:', 'LineWidth', 2);
    plot(t_axis, P_total_A, 'b--', 'LineWidth', 1.5);
    plot(t_axis, P_total_B, 'r-', 'LineWidth', 1.5);
    ylabel('调节功率 (MW)'); legend('指令', '中性策略', '规避策略'); grid on;
    title('调度计划追踪对比 (正:上调, 负:下调)');
    
    subplot(2,1,2); 
    bar(t_axis, [P_sl_A, P_sl_B]*baseMVA);
    legend('中性切负荷', '规避切负荷'); ylabel('缺额 (MW)');
    title('违约风险对比');
    
    print(fig, 'Test3_Risk.png', '-dpng', '-r300');
end