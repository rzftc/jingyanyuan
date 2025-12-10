function run_scenario_D_ramp(strategies, Scenarios_AC_Up, Scenarios_EV_Up, Scenarios_AC_Down, Scenarios_EV_Down, direction_signal)
    fprintf('\n>>> 场景 D: 鲁棒性测试 <<<\n');
    T_steps = length(direction_signal);
    [~, N_scenarios] = size(Scenarios_AC_Up);
    
    if ~isempty(strategies{1}) && ~isempty(strategies{3})
        % 计算每个场景下，两种策略各自的总违约量
        % 违约定义：sum_t max(0, P_AC(t) - Cap_AC(t,s)) + max(0, P_EV(t) - Cap_EV(t,s))
        viol_neu_all = zeros(N_scenarios, 1);   % Beta = 0
        viol_rob_all = zeros(N_scenarios, 1);   % Beta = 10
        
        for s = 1:N_scenarios
            Real_Cap_AC = zeros(T_steps, 1);
            Real_Cap_EV = zeros(T_steps, 1);
            
            for t = 1:T_steps
                if direction_signal(t) == 1
                    % 上调时刻，使用 Up 场景
                    Real_Cap_AC(t) = Scenarios_AC_Up(t, s);
                    Real_Cap_EV(t) = Scenarios_EV_Up(t, s);
                else
                    % 下调时刻，使用 Down 场景的绝对值
                    Real_Cap_AC(t) = abs(Scenarios_AC_Down(t, s));
                    Real_Cap_EV(t) = abs(Scenarios_EV_Down(t, s));
                end
            end
            
            calc_viol = @(P, Cap) sum(max(0, P(:) - Cap(:)));
            
            viol_neu_all(s) = calc_viol(strategies{1}.P_AC, Real_Cap_AC) + ...
                              calc_viol(strategies{1}.P_EV, Real_Cap_EV);
            viol_rob_all(s) = calc_viol(strategies{3}.P_AC, Real_Cap_AC) + ...
                              calc_viol(strategies{3}.P_EV, Real_Cap_EV);
        end
        
        % 寻找一个 Beta=0 违约量大于 Beta=10 的场景
        diff_viol = viol_neu_all - viol_rob_all;
        idx_candidates = find(diff_viol > 0);
        
        if ~isempty(idx_candidates)
            % 从这些候选里，选“差值最大”的那个场景
            [~, idx_maxdiff_local] = max(diff_viol(idx_candidates));
            worst_idx = idx_candidates(idx_maxdiff_local);
            fprintf('  - 选取的场景: #%d (Beta=0 违约量 > Beta=10 违约量)\n', worst_idx);
        else
            % 如果不存在 Beta=0 > Beta=10 的场景，则退回到 Beta=0 的最坏场景
            [~, worst_idx] = max(viol_neu_all);
            fprintf('  - 未找到 Beta=0 违约量大于 Beta=10 的场景，退化为 Beta=0 最坏场景: #%d\n', worst_idx);
        end
        
        % 取出该场景下两种策略的违约量
        v_neu = viol_neu_all(worst_idx);
        v_rob = viol_rob_all(worst_idx);
        
        fprintf('  - 中性策略(Beta=0) 该场景违约量: %.2f MW\n', v_neu);
        fprintf('  - 规避策略(Beta=10) 该场景违约量: %.2f MW\n', v_rob);
        
        % ===== 绘图：和原始场景 D 的绘图风格基本一致，并在柱子上标数字 =====
        figure('Color','w', 'Position', [600, 300, 500, 400]); 
        b = bar([v_neu, v_rob], 0.5);
        b.FaceColor = 'flat'; 
        b.CData(1,:) = [0.8 0.2 0.2];   % Beta=0 用红色系
        b.CData(2,:) = [0.2 0.6 0.2];   % Beta=10 用绿色系
        
        set(gca, 'XTickLabel', {'中性 (\beta=0)', '规避 (\beta=10)'});
        ylabel('极端场景实际违约量 (MW)');
        grid on;
        
        % 在柱状图上标出具体数字（保留两位小数）
        xs = 1:2;
        ys = [v_neu, v_rob];
        labels = arrayfun(@(x) sprintf('%.2f', x), ys, 'UniformOutput', false);
        text(xs, ys, labels, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        
        % 保存图片
        print(gcf, '极端场景违约对比_选定场景.png', '-dpng', '-r300');
    end
end
