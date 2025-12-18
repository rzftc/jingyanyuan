function run_test2_uncertainty(t_axis, S_AC_U, S_EV_U, S_AC_D, S_EV_D, R_AC_U, R_EV_U, R_AC_D, R_EV_D)
% 绘制调节能力云图和可靠边界 (含上下调)

    fig = figure('Name', 'Test2_Uncertainty', 'Color', 'w', 'Position', [100, 100, 1000, 800]);
    
    % 抽样绘图索引
    [~, N] = size(S_AC_U);
    idx = 1:10:N; 
    
    % --- 子图1: AC 上调 ---
    subplot(2, 2, 1); hold on;
    plot(t_axis, S_AC_U(:, idx), 'Color', [0.8 0.8 1]);
    plot(t_axis, R_AC_U, 'b-', 'LineWidth', 2);
    title('AC 上调能力 (增加负荷)'); ylabel('MW'); xlim([8, 32]); grid on;
    
    % --- 子图2: EV 上调 ---
    subplot(2, 2, 2); hold on;
    plot(t_axis, S_EV_U(:, idx), 'Color', [1 0.8 0.8]);
    plot(t_axis, R_EV_U, 'r-', 'LineWidth', 2);
    title('EV 上调能力 (增加负荷)'); xlim([8, 32]); grid on;
    
    % --- 子图3: AC 下调 ---
    subplot(2, 2, 3); hold on;
    plot(t_axis, S_AC_D(:, idx), 'Color', [0.8 0.8 1]);
    plot(t_axis, R_AC_D, 'b-', 'LineWidth', 2);
    title('AC 下调能力 (减少负荷)'); ylabel('MW'); xlim([8, 32]); grid on;
    
    % --- 子图4: EV 下调 ---
    subplot(2, 2, 4); hold on;
    plot(t_axis, S_EV_D(:, idx), 'Color', [1 0.8 0.8]);
    plot(t_axis, R_EV_D, 'r-', 'LineWidth', 2);
    title('EV 下调能力 (减少负荷)'); xlim([8, 32]); grid on;
    
    print(fig, 'Test2_Uncertainty.png', '-dpng', '-r300');
end