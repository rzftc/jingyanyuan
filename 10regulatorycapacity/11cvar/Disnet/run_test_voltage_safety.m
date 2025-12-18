function run_test_voltage_safety(t_axis, net, cost_p, R_AC_Up, R_EV_Up, baseMVA)
% run_test_voltage_safety
% 验证 LinDistFlow 约束是否有效防止电压越限
% 场景：施加极端的“上调”指令（增加负荷），这是造成配网电压跌落的最恶劣工况

    fprintf('  - 正在进行电压安全验证 (极限增负荷测试)...\n');
    T = length(t_axis);
    
    % 构造极限压力测试指令：全时段按最大能力的 90% 进行增负荷
    P_stress_pu = 0.9 * (R_AC_Up + R_EV_Up) / baseMVA;
    dir_vec = ones(T, 1); % 全程上调
    
    % 场景设定
    S_AC = repmat(R_AC_Up, 1, 1);
    S_EV = repmat(R_EV_Up, 1, 1);
    risk_p.beta = 10; risk_p.confidence = 0.95; risk_p.rho_pen = 100;
    
    % 求解
    [~, ~, ~, U_opt, ~] = solve_dist_optimization_corrected(P_stress_pu, net, cost_p, risk_p, S_AC, S_EV, R_AC_Up, R_EV_Up, 0, 0, dir_vec);
    
    V_opt = sqrt(U_opt); % 转换为电压幅值
    
    % 提取关键节点 (末端节点最易越限)
    idx_comm = 18; % 商业区末端
    idx_res  = 33; % 居民区末端
    
    fig = figure('Name', 'VoltageSafety', 'Color', 'w', 'Position', [100, 100, 800, 600]);
    
    subplot(2,1,1);
    [X,Y] = meshgrid(1:net.N_bus, t_axis);
    surf(X, Y, V_opt); shading interp; colorbar;
    xlabel('Bus ID'); ylabel('Time (h)'); zlabel('Voltage (p.u.)');
    title('全网节点电压分布 (极限负荷工况)');
    zlim([0.94, 1.05]);
    
    subplot(2,1,2); hold on;
    plot(t_axis, V_opt(:, idx_comm), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Bus 18 (Comm End)');
    plot(t_axis, V_opt(:, idx_res), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Bus 33 (Res End)');
    yline(0.95, 'k--', 'Limit 0.95');
    ylabel('Voltage (p.u.)'); legend('Location', 'best'); grid on;
    title('关键末端节点电压曲线');
    
    min_v = min(V_opt(:));
    if min_v >= 0.949
        fprintf('  > 电压验证通过! 最低电压: %.4f p.u.\n', min_v);
    else
        fprintf('  > 警告: 存在电压越限, 最低: %.4f p.u.\n', min_v);
    end
    
    print(fig, 'Test_Voltage_Safety.png', '-dpng', '-r300');
end