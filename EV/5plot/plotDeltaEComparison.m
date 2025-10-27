function plotDeltaEComparison(EVs)
    % 创建新画布
    fig = figure('Name','Delta_E对比分析', 'Position',[200 200 1200 800]);
    
    % 参数设置（仅取前10辆EV）
    num_plot = min(10, length(EVs)); % 关键修改点
    selected_EVs = EVs(1:num_plot);
    num_points = 50;
    line_colors = jet(num_plot); % 根据实际数量生成颜色
    
    % 主绘图循环
    for ev_idx = 1:num_plot
        ev = selected_EVs(ev_idx);
        P_range = linspace(ev.P_l_min*0.9, ev.P_h_max*1.1, num_points);
        Delta_E = arrayfun(@(p) calculateDeltaE(ev,p), P_range);
        
        % 绘制曲线（带透明度）
        plot(P_range, Delta_E, 'LineWidth',1.8,...
             'Color',[line_colors(ev_idx,:) 0.7],... % 添加透明度
             'DisplayName',sprintf('EV%d',ev_idx))
        hold on
        
        % 标记特征点
        plot([ev.P_l_min, ev.P_0, ev.P_h_max],...
             [ev.Delta_E_q_max, 0, ev.Delta_E_h_max],...
             '^','MarkerSize',8,...
             'MarkerFaceColor',line_colors(ev_idx,:),...
             'MarkerEdgeColor','k')
    end
    
    % 图形标注
    hold off
    xlabel('实时电价 p_{real} (元/℃)','FontSize',24)
    ylabel('\Delta T (℃)','FontSize',24)
    legend('show','Location','eastoutside','FontSize',10)
    grid on
    
    % 设置坐标轴范围
    x_limits = [min([selected_EVs.P_l_min])*0.9,...
               max([selected_EVs.P_h_max])*1.1];
    y_limits = [min([selected_EVs.Delta_E_q_max])*1.1,...
               max([selected_EVs.Delta_E_h_max])*1.1];
    xlim(x_limits)
    ylim(y_limits)
    ax = gca;
    ax.FontSize = 18;       % 统一设置刻度字体
    ax.LabelFontSizeMultiplier = 1.3; % 标签字体放大系数

    % 保存图像
end
