function plotSensitivityCurves(p, mUp, sUp, mDown, sDown,...
    p_min, p_max, p_min_prime, p_max_prime, cIdx, bUp, ~)
    
    %% 维度一致性处理
    p = p(:).';  % 强制转换为行向量
    mUp = mUp(:).'; sUp = sUp(:).';
    mDown = mDown(:).'; sDown = sDown(:).';
    
    %% 上调潜力图绘制与保存
    figure('Position', [200 200 1400 700]);
    hold on;
    
    % 置信区间填充
    upper = mUp + sUp;
    lower = mUp - sUp;
    fill([p, fliplr(p)], [upper, fliplr(lower)], [0.8 0.8 1], 'EdgeColor','none');
    
    % 均值曲线和临界点
    plot(p, mUp, 'b-', 'LineWidth', 2);
    scatter(p(cIdx), mUp(cIdx), 100, 'r', 'filled');
    
    % 右侧经济收益曲线
    yyaxis right;
    plot(p, bUp, 'g--', 'LineWidth', 1.5);
    
    % 坐标轴与字体设置（加大字体）
    yyaxis left;
    xlabel('激励电价 (元)', 'FontSize', 16);
    ylabel('灵敏度 (kW/元)', 'FontSize', 16);
    set(gca, 'FontSize', 16);  % 刻度字体
    grid on;
    
    yyaxis right;
    ylabel('经济收益 (元/kW)', 'FontSize', 16);
    set(gca, 'FontSize', 16);
    
    % 理论边界线
    xline([p_min, p_max, p_min_prime, p_max_prime], 'LineStyle', '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
    
    % 图例（加大字体）
    legend('置信区间', '灵敏度均值', '临界点', 'Location', 'northwest', 'FontSize', 16);
    
    % 边界说明注释
    annotation('textbox', [0.15 0.15 0.2 0.1],...
        'String', {sprintf('原始范围: %.1f-%.1f', p_min, p_max),...
        sprintf('调整范围: %.1f-%.1f', p_min_prime, p_max_prime)},...
        'FitBoxToText', 'on',...
        'BackgroundColor', [1 1 1 0.8]);
    
    % 保存400dpi PNG
    print('上调灵敏度特性', '-dpng', '-r400');
    
    %% 下调潜力图绘制与保存
    figure('Position', [200 200 1400 700]);
    hold on;
    
    % 置信区间填充
    upper_down = mDown + sDown;
    lower_down = mDown - sDown;
    fill([p, fliplr(p)], [upper_down, fliplr(lower_down)], [1 0.8 0.8], 'EdgeColor','none');
    
    % 均值曲线和临界点
    plot(p, mDown, 'r-', 'LineWidth', 2);
    scatter(p(cIdx), mDown(cIdx), 100, 'b', 'filled');
    
    % 坐标轴与字体设置（加大字体）
    yyaxis left;
    xlabel('激励电价 (元)', 'FontSize', 16);
    ylabel('灵敏度 (kW/元)', 'FontSize', 16);
    set(gca, 'FontSize', 12);  % 刻度字体
    grid on;
    
    yyaxis right;
    ylabel('经济收益 (元/kW)', 'FontSize', 16);
    set(gca, 'FontSize', 16);
    
    % 理论边界线
    xline([p_min, p_max, p_min_prime, p_max_prime], 'LineStyle', '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
    
    % 图例（加大字体）
    legend('置信区间', '灵敏度均值', '临界点', 'Location', 'northwest', 'FontSize', 16);
    
    % 边界说明注释
    annotation('textbox', [0.15 0.15 0.2 0.1],...
        'String', {sprintf('原始范围: %.1f-%.1f', p_min, p_max),...
        sprintf('调整范围: %.1f-%.1f', p_min_prime, p_max_prime)},...
        'FitBoxToText', 'on',...
        'BackgroundColor', [1 1 1 0.8]);
    
    % 保存400dpi PNG
    print('下调灵敏度特性', '-dpng', '-r400');
end
