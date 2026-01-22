function run_scenario_C(strategies, t_axis, P_grid_demand, direction_signal, ...
    Reliable_AC_Up, Reliable_EV_Up, Reliable_AC_Down, Reliable_EV_Down, ...
    x_ticks_set, x_labels_set)

    fprintf('\n>>> 绘制调度方案 (Beta=10) <<<\n');
    idx_plot = 2; 
    T_steps = length(t_axis);
    
    if ~isempty(strategies{idx_plot})
        st = strategies{idx_plot};
        dir_sign = sign(direction_signal); dir_sign(dir_sign==0)=1;
        
        Y_Stack_Plot = [st.P_AC, st.P_EV, st.P_Gen, st.P_Shed] .* dir_sign;
        P_Demand_plot = P_grid_demand .* dir_sign;
        
        % --- 图 1: 功率堆叠 ---
        fig_stack = figure('Name', '多源协同调度堆叠图', 'Color', 'w', 'Position', [150, 150, 1000, 600]);
        hold on;
        h_area = area(t_axis, Y_Stack_Plot);
        h_area(1).FaceColor = [0.00, 0.45, 0.74]; h_area(1).EdgeColor = 'none'; % AC
        h_area(2).FaceColor = [0.47, 0.67, 0.19]; h_area(2).EdgeColor = 'none'; % EV
        h_area(3).FaceColor = [0.92, 0.69, 0.13]; h_area(3).EdgeColor = 'none'; % Gen
        h_area(4).FaceColor = [0.85, 0.33, 0.10]; h_area(4).EdgeColor = 'none'; h_area(4).FaceAlpha = 0.8; % Shed
        plot(t_axis, P_Demand_plot, 'k--', 'LineWidth', 2.0, 'DisplayName', '电网总需求');
        
        % [修改] 坐标轴标签大字号加粗
        ylabel('功率 (MW) [正=增加用电, 负=减少用电]', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Microsoft YaHei');
        % [修改] 图例大字号
        legend([h_area(1), h_area(2), h_area(3), h_area(4)], ...
               {'空调 (AC)', '电动汽车 (EV)', '火电调节 (Gen)', '切负荷/弃光(Shed)'}, ...
               'Location', 'best', 'FontSize', 18, 'FontName', 'Microsoft YaHei');
        
        % [修改] 去除网格线
        % grid on; % 已注释
        xlim([8, 32]); 
        % [修改] 坐标轴刻度大字号，框线加粗
        set(gca, 'XTick', x_ticks_set, 'XTickLabel', x_labels_set, 'FontSize', 18, 'LineWidth', 1.5, 'FontName', 'Microsoft YaHei');
        print(fig_stack, '多源协同调度堆叠图.png', '-dpng', '-r600');
        
        % --- 图 2: AC与EV时序出力对比 ---
        fprintf('  >>> 正在绘制 AC与EV时序调度量对比图...\n');
        fig_comp = figure('Name', 'AC与EV时序出力对比', 'Color', 'w', 'Position', [200, 200, 1000, 600]);
        hold on;
        
        Reliable_AC_plot = zeros(T_steps, 1);
        Reliable_EV_plot = zeros(T_steps, 1);
        for t = 1:T_steps
            if dir_sign(t) >= 0
                Reliable_AC_plot(t) = Reliable_AC_Up(t);
                Reliable_EV_plot(t) = Reliable_EV_Up(t);
            else
                Reliable_AC_plot(t) = -Reliable_AC_Down(t); 
                Reliable_EV_plot(t) = -Reliable_EV_Down(t);
            end
        end
        P_AC_line = st.P_AC .* dir_sign;
        P_EV_line = st.P_EV .* dir_sign;

        % [修改] 曲线线宽统一为 2.0
        plot(t_axis, Reliable_AC_plot, 'b:', 'LineWidth', 2.0, 'DisplayName', 'AC 可靠边界');
        plot(t_axis, P_AC_line, 'b-', 'LineWidth', 2.0, 'DisplayName', 'AC 调度');
        plot(t_axis, Reliable_EV_plot, 'g:', 'LineWidth', 2.0, 'DisplayName', 'EV 可靠边界');
        plot(t_axis, P_EV_line, 'g-', 'LineWidth', 2.0, 'DisplayName', 'EV 调度');
        yline(0, 'k-', 'HandleVisibility', 'off', 'LineWidth', 1.5);
        
        % [修改] 坐标轴标签大字号加粗
        ylabel('调节功率 (MW)', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Microsoft YaHei'); 
        xlabel('时间', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Microsoft YaHei');
        % [修改] 图例大字号
        legend('Location', 'best', 'FontSize', 18, 'FontName', 'Microsoft YaHei'); 
        
        % [修改] 去除网格线
        % grid on; % 已注释
        xlim([8, 32]); 
        % [修改] 坐标轴刻度大字号，框线加粗
        set(gca, 'XTick', x_ticks_set, 'XTickLabel', x_labels_set, 'FontSize', 18, 'LineWidth', 1.5, 'FontName', 'Microsoft YaHei');
        print(fig_comp, 'AC与EV时序调度量对比_Optimized.png', '-dpng', '-r600');
        
        % --- 指标迭代对比图 ---
        if isfield(st, 'SDCI_History')
            fig_sdci = figure('Name', 'SDCI 迭代对比', 'Color', 'w', 'Position', [300, 300, 600, 400]);
            sdci_vals = [st.SDCI_History(1), st.SDCI_History(end)];
            bar(sdci_vals, 0.4, 'FaceColor', [0.2 0.6 0.8]);
            
            % [修改] 坐标轴刻度大字号，框线加粗
            set(gca, 'XTickLabel', {'初始迭代', '最终迭代'}, 'FontSize', 18, 'LineWidth', 1.5, 'FontName', 'Microsoft YaHei');
            % [修改] 坐标轴标签大字号加粗
            ylabel('SDCI 指标', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Microsoft YaHei');
            % [修改] 数值标签字号调整为 16
            text(1:2, sdci_vals, num2str(sdci_vals', '%.4f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 16);
            print(fig_sdci, 'SDCI对比.png', '-dpng', '-r600');
            
            fig_rho = figure('Name', 'Rho 迭代对比', 'Color', 'w', 'Position', [400, 400, 600, 400]);
            rho_vals = [st.Rho_History(1), st.Rho_History(end)];
            bar(rho_vals, 0.4, 'FaceColor', [0.8 0.4 0.2]);
            
            % [修改] 坐标轴刻度大字号，框线加粗
            set(gca, 'XTickLabel', {'初始迭代', '最终迭代'}, 'FontSize', 18, 'LineWidth', 1.5, 'FontName', 'Microsoft YaHei');
            % [修改] 坐标轴标签大字号加粗
            ylabel('Rho 指标', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Microsoft YaHei');
            % [修改] 数值标签字号调整为 16
            text(1:2, rho_vals, num2str(rho_vals', '%.4f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 16);
            print(fig_rho, 'Rho对比.png', '-dpng', '-r600');
        end
    end
end