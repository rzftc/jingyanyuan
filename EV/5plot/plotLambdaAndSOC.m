
function plotLambdaAndSOC(results, dt_short, ev_idx)
    % 图1：lambda与单台EV的SOC协同分析
    time_hours = (0:length(results.lambda)-1)*dt_short/60;
    
    fig = figure('Name',sprintf('EV%d-Lambda&SOC协同',ev_idx),...
        'Position',[100 100 1000 400], 'NumberTitle','off');
    
    % ===== 主曲线 =====
    % 左侧坐标轴（SOC）
    yyaxis left
    main_soc = plot(time_hours, results.EV_S_original(ev_idx,:),...
        'LineWidth',1.2, 'Color',[0.8 0.2 0.2], 'DisplayName','SOC原始值');
    ylabel('SOC (-1~1)', 'FontSize',16, 'Color',[0.8 0.2 0.2])
    ylim([-1.1 1.1])
    set(gca, 'YColor', [0.8 0.2 0.2])
    
    % 右侧坐标轴（Lambda）
    yyaxis right
    main_lambda = plot(time_hours, results.lambda,...
        'LineWidth',1.2, 'Color',[0.2 0.4 0.8], 'DisplayName','\lambda^*');
    ylabel('\lambda^*', 'FontSize',16, 'Color',[0.2 0.4 0.8])
    ylim([floor(min(results.lambda)*2)/2, ceil(max(results.lambda)*2)/2])
    
    % ===== 趋势线 =====
    hold on
    % SOC趋势（7点移动平均）
    trend_soc = plot(time_hours, movmean(results.EV_S_original(ev_idx,:),7),...
        '-.', 'LineWidth',2, 'Color',[0.5 0 0], 'DisplayName','SOC趋势');
    
    % Lambda趋势（7点移动平均）
    trend_lambda = plot(time_hours, movmean(results.lambda,7),...
        ':', 'LineWidth',2, 'Color',[0 0 0.6], 'DisplayName','\lambda趋势');
    hold off
    
    % ===== 公共设置 =====
    xlabel('时间 (小时)', 'FontSize',16)
    set(gca, 'FontSize',12) % 坐标刻度字号12
    xlim([0 24])
    grid on
    legend([main_soc, main_lambda, trend_soc, trend_lambda],...
        'Location','northwest', 'FontSize',14)
    
    % ===== 保存图像 =====
    print(fig, 'EV%d_LambdaSOC3.png', '-dpng', '-r600')
end
