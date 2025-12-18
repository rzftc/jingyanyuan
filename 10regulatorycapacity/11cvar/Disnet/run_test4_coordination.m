function run_test4_coordination(t_axis, P_req_pu, dir_vec, net, cost_p, ...
    R_AC_U, R_EV_U, R_AC_D, R_EV_D, baseMVA)

    fprintf('  - 正在分析协同效果...\n');
    T = length(t_axis);
    
    % 使用可靠域作为确定性场景进行协同演示
    S_AC_mean = zeros(T, 1); S_EV_mean = zeros(T, 1);
    R_AC_curr = zeros(T, 1); R_EV_curr = zeros(T, 1);
    
    for t = 1:T
        if dir_vec(t) == 1
            S_AC_mean(t) = R_AC_U(t); S_EV_mean(t) = R_EV_U(t);
            R_AC_curr(t) = R_AC_U(t); R_EV_curr(t) = R_EV_U(t);
        else
            S_AC_mean(t) = R_AC_D(t); S_EV_mean(t) = R_EV_D(t);
            R_AC_curr(t) = R_AC_D(t); R_EV_curr(t) = R_EV_D(t);
        end
    end
    
    risk_p.beta = 1; risk_p.confidence = 0.95; risk_p.rho_pen = 100;

    % 无协同
    [P_AC_NC, P_EV_NC] = solve_dist_optimization_corrected(P_req_pu, net, cost_p, risk_p, S_AC_mean, S_EV_mean, R_AC_curr, R_EV_curr, 0, 0, dir_vec);
    
    % 有协同
    lambda = 2000;
    [P_AC_C, P_EV_C] = solve_dist_optimization_corrected(P_req_pu, net, cost_p, risk_p, S_AC_mean, S_EV_mean, R_AC_curr, R_EV_curr, lambda, lambda, dir_vec);
    
    % 绘图
    fig = figure('Name', 'Test4_Coord', 'Color', 'w', 'Position', [100, 100, 1000, 400]);
    
    subplot(1,2,1); 
    % 使用 area 绘图需要处理负值，这里只画绝对值堆叠
    area(t_axis, [P_AC_NC, P_EV_NC]*baseMVA);
    title('无协同'); legend('AC', 'EV'); ylabel('MW'); grid on; ylim([0, 3]);
    
    subplot(1,2,2); 
    area(t_axis, [P_AC_C, P_EV_C]*baseMVA);
    title('协同优化'); legend('AC', 'EV'); grid on; ylim([0, 3]);
    
    print(fig, 'Test4_Coordination.png', '-dpng', '-r300');
end