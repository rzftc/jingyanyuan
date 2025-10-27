
function plotPowerComparison(results, dt_short)
    % 图2：功率对比分析
    time_hours = (0:length(results.P_agg)-1)*dt_short/60;
    
    fig = figure('Name','功率跟踪效果分析',...
        'Position',[100 100 1000 400], 'NumberTitle','off');
    
    % ===== 主曲线 =====
    % 目标功率（橙色半透明区域）
    area(time_hours, results.P_tar,...
        'FaceColor',[1 0.6 0], 'FaceAlpha',0.3, 'EdgeColor','none',...
        'DisplayName','目标功率')
    hold on
    
    % 实际功率（深绿色实线）
    main_agg = plot(time_hours, results.P_agg,...
        'LineWidth',1.5, 'Color',[0 0.5 0], 'DisplayName','实际功率');
    
    % ===== 趋势线 =====
    % 实际功率趋势（15点移动平均）
    trend_agg = plot(time_hours, movmean(results.P_agg,15),...
        '-', 'LineWidth',2, 'Color',[0 0.2 0], 'DisplayName','功率趋势');
    hold off
    
    % ===== 坐标轴设置 =====
    xlabel('时间 (小时)', 'FontSize',14)
    set(gca, 'FontSize',12) % 坐标刻度字号12
    ylabel('功率 (kW)', 'FontSize',14)
    xlim([0 24])
    ylim([0 max([results.P_agg, results.P_tar])*1.1])
    grid on
    legend('Location','northeast', 'FontSize',12)
    
    % ===== 保存图像 =====
    print(fig, 'PowerComparison3.png', '-dpng', '-r600')
end
