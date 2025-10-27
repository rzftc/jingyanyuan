
function plotLambdaAndAggSOC(results, dt_short)
    % 图3：lambda与聚合SOC协同分析
    time_hours = (0:length(results.lambda)-1)*dt_short/60;
    
    fig = figure('Name','Lambda与聚合SOC协同',...
        'Position',[100 100 1000 400], 'NumberTitle','off');
    
    % ===== 主曲线 =====
    % 左侧坐标轴（聚合SOC）
    yyaxis left
    main_soc = plot(time_hours, results.S_agg,...
        'LineWidth',1.2, 'Color',[0.1 0.5 0.2], 'DisplayName','聚合SOC');
    ylabel('聚合SOC (-1~1)', 'FontSize',16, 'Color',[0.1 0.5 0.2])
    ylim([-1.1 1.1])
    set(gca, 'YColor', [0.1 0.5 0.2])
    
    % 右侧坐标轴（Lambda）
    yyaxis right
    main_lambda = plot(time_hours, results.lambda,...
        'LineWidth',1.2, 'Color',[0.2 0.4 0.8], 'DisplayName','\lambda^*');
    ylabel('\lambda^*', 'FontSize',16, 'Color',[0.2 0.4 0.8])
    ylim([floor(min(results.lambda)*2)/2, ceil(max(results.lambda)*2)/2])
    
    % ===== 趋势线 =====
    hold on
    % SOC趋势（7点移动平均）
    trend_soc = plot(time_hours, movmean(results.S_agg,7),...
        '-.', 'LineWidth',2, 'Color',[0 0.3 0], 'DisplayName','SOC趋势');
    
    % Lambda趋势（7点移动平均）
    trend_lambda = plot(time_hours, movmean(results.lambda,7),...
        ':', 'LineWidth',2, 'Color',[0 0 0.6], 'DisplayName','\lambda趋势');
    hold off
    
    % ===== 公共设置 =====
    xlabel('时间 (小时)', 'FontSize',14)
    set(gca, 'FontSize',12) % 坐标刻度字号12
    xlim([0 24])
    grid on
    legend([main_soc, main_lambda, trend_soc, trend_lambda],...
        'Location','best', 'FontSize',12)
    
    % ===== 保存图像 =====
    print(fig, 'Lambda_AggSOCv3.png', '-dpng', '-r600')
end