function plotResults(results, dt_short)
    time_hours = (0:dt_short:length(results.lambda)-1)/60;
    
    figure('Position', [200 200 1200 600]);
    
    % 主图：功率跟踪与λ*
    subplot(2,1,1);
    yyaxis left
    plot(time_hours, results.P_agg, 'b', 'LineWidth',1.5);
    hold on;
    stairs(time_hours(1:10:end), results.P_base(1:10:end), 'r--', 'LineWidth',1.2);
    ylabel('功率 (kW)');
    
    yyaxis right
    plot(time_hours, results.lambda, 'g', 'LineWidth',1.2);
    ylabel('λ*');
    legend('实际功率','基准功率','λ*');
    title('功率跟踪性能');
    
    % 子图：SOC动态
    subplot(2,1,2);
    plot(time_hours, results.S_agg, 'k', 'LineWidth',1.5);
    hold on;
    plot(time_hours, results.EV_S_original, 'm:', 'LineWidth',1.2);
    legend('聚合虚拟SOC','EV5原始SOC');
    title('SOC动态对比');
    xlabel('时间 (小时)');
end
