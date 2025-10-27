function plotEVStatusEnhanced(results, dt_short, ev_idx)
    % 时间轴生成（统一小时单位）
    time_hours = (0:length(results.lambda)-1)*dt_short/60; 
    
    % 创建三轴坐标系
    figure('Name',sprintf('EV%d协同分析',ev_idx), 'Position',[100 100 1200 500])
    
    % 主坐标轴（左侧）- SOC原始值
    yyaxis left
    soc_plot = plot(time_hours, results.EV_S_original(ev_idx,:),...
        'LineWidth',1.8, 'Color',[0.85 0.2 0.2], 'DisplayName','EV SOC');
    ylabel('EV SOC (-1~1)', 'FontSize',12, 'Color',[0.85 0.2 0.2])
    ylim([-1.1 1.1])
    set(gca, 'YColor', [0.85 0.2 0.2])
    
    % 右侧第一坐标轴 - 聚合SOC
    yyaxis right
    ax = gca;
    ax.YAxis(1).Visible = 'off'; % 隐藏默认右侧轴
    
    % 创建第二右侧坐标轴
    ax_pos = ax.Position;
    ax2 = axes('Position',ax_pos,...
        'YAxisLocation','right',...
        'Color','none');
    linkprop([ax ax2], {'Position','XLim'}); % 位置联动
    
    % 绘制聚合SOC
    hold(ax2, 'on')
    agg_plot = plot(ax2, time_hours, results.S_agg,...
        '-.', 'LineWidth',1.5, 'Color',[0.1 0.6 0.3], 'DisplayName','聚合SOC');
    ylabel(ax2, '聚合SOC', 'FontSize',12, 'Color',[0.1 0.6 0.3])
    set(ax2, 'YColor', [0.1 0.6 0.3])
    
    % 绘制Lambda（共用第二右侧轴）
    lambda_plot = plot(ax2, time_hours, results.lambda,...
        '--', 'LineWidth',1.2, 'Color',[0.2 0.4 0.8], 'DisplayName','λ^*');
    
    % 坐标轴设置
    xlabel('时间 (小时)', 'FontSize',12)
    xlim([0 24])
    title(sprintf('EV%d状态协同分析',ev_idx), 'FontSize',14)
    grid on
    
    % 动态刻度
    y_range = [floor(min([results.lambda, results.S_agg])*2)/2,...
              ceil(max([results.lambda, results.S_agg])*2)/2];
    ylim(ax2, y_range)
    
    % 图例整合
    leg = legend([soc_plot, agg_plot, lambda_plot],...
        'Location','northwest');
    set(leg, 'FontSize',10)
   
end
