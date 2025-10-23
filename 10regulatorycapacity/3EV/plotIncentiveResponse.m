function plotIncentiveResponse()
    % 参数设置（示例值，可根据实际情况调整）
    p_min = 15;       % 原始电价下限
    p_max = 50;       % 原始电价上限
    p_min_prime = 10; % 修正后电价下限
    p_max_prime = 35; % 修正后电价上限
    E_tar_max = 100;  % 最大能量调节量
    
    % 生成激励电价范围（覆盖所有关键点）
    p = linspace(0, 60, 500); % 扩展到60以显示饱和区
    
    % 调用激励响应函数
    [E_up, E_down, ~] = incentiveTempEV(p, p_min, p_max, p_min_prime, p_max_prime, E_tar_max*ones(size(p)));
    
    % 创建专业可视化图形
    figure('Position', [100 100 900 500], 'Color', 'w')
    
    % 绘制主曲线
    plot(p, E_up,  'LineWidth', 2.5, 'Color', [0.2 0.6 0.8], 'DisplayName', '上调潜力 E_{tar}^{up}')
    hold on
    plot(p, E_down, 'LineWidth', 2.5, 'Color', [0.8 0.4 0.1], 'LineStyle', '--', 'DisplayName', '下调潜力 E_{tar}^{down}')
    
    % 标注关键参数点
    plotCriticalPoints(p_min, p_max, p_min_prime, p_max_prime, E_tar_max)
    
    % 图形美化
    xlabel('激励电价 p (元)', 'FontSize', 12, 'FontWeight', 'bold')
    ylabel('能量调节潜力 (kWh)', 'FontSize', 12, 'FontWeight', 'bold')
    title('虚拟电厂调节潜力特性曲线', 'FontSize', 14, 'FontWeight', 'bold')
    legend('Location', 'northwest')
    grid on
    axis tight
    set(gca, 'FontSize', 11, 'LineWidth', 1.2)
    
    % 添加理论公式标注
    annotation('textbox', [0.15 0.7 0.2 0.15], 'String', ...
        {'$$\Delta Z_{set}^{up} = \frac{p - p''_{min}}{p''_{max} - p''_{min}} \Delta Z_{max}$$', ...
         '$$\Delta Z_{set}^{down} = \frac{p - p_{min}}{p_{max} - p_{min}} \Delta Z_{max}$$'}, ...
         'Interpreter', 'latex', 'FontSize', 12, 'EdgeColor', 'none')
end

function plotCriticalPoints(p_min, p_max, p_min_prime, p_max_prime, E_max)
    % 绘制关键转折点
    xline(p_min_prime, 'k:', 'LineWidth', 1.5, 'DisplayName', 'p''_{min}')
    xline(p_max_prime, 'k:', 'LineWidth', 1.5, 'DisplayName', 'p''_{max}')
    xline(p_min, 'k-.', 'LineWidth', 1.5, 'DisplayName', 'p_{min}')
    xline(p_max, 'k-.', 'LineWidth', 1.5, 'DisplayName', 'p_{max}')
    
    % 标注饱和区域
    text(p_max_prime+2, E_max*0.95, '饱和区', ...
        'FontSize', 10, 'Color', [0.5 0.5 0.5], 'Rotation', 90)
    text(p_max+2, E_max*0.95, '饱和区', ...
        'FontSize', 10, 'Color', [0.5 0.5 0.5], 'Rotation', 90)
    
    % 添加特征点标记
    scatter([p_min_prime, p_max_prime], [0, E_max], 80, ...
        'MarkerEdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.5, ...
        'DisplayName', '特征点')
end
