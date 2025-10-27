function plotEnhancedResults(results, dt_short)
    % 创建时间轴（单位：分钟）
    time_min = (0:length(results.lambda)-1)*dt_short;
    
    % 创建双轴画布
    figure('Position', [100 100 1200 600])
    yyaxis left
    
    %% 绘制所有EV的SOC曲线（半透明处理）
    hold on
    for ev = 1:size(results.EV_S_original,1)
        plot(time_min, results.EV_S_original(ev,:),...
            'Color',[0.2 0.4 0.8 0.05], 'LineWidth',0.3) % 半透明浅蓝色
    end
    
    %% 绘制SOC统计特征
    soc_median = median(results.EV_S_original, 1);
    soc_upper = quantile(results.EV_S_original, 0.75, 1);
    soc_lower = quantile(results.EV_S_original, 0.25, 1);
    
    % 绘制四分位区域
    fill([time_min fliplr(time_min)],...
         [soc_upper fliplr(soc_lower)],...
         [0.4 0.6 1], 'FaceAlpha',0.2, 'EdgeColor','none')
    
    % 绘制中位线
    plot(time_min, soc_median, 'b-', 'LineWidth',1.5)
    
    %% 设置左轴属性
    ylabel('SOC')
    ylim([0 1])
    grid on
    
    %% 绘制Lambda曲线
    yyaxis right
    plot(time_min, results.lambda, 'r--', 'LineWidth',1.5)
    ylabel('\lambda^*')
    ylim([0 1.1])
    
    %% 通用设置
    xlabel('时间 (分钟)')
    title('电动汽车SOC分布与最优Lambda参数')
    legend('各车SOC', '25-75%分位区', 'SOC中位数', '\lambda^*',...
           'Location','southeast')
    
    set(gcf, 'Color','w', 'InvertHardcopy','off')
    set(gca, 'FontSize',10)
    box on
end
