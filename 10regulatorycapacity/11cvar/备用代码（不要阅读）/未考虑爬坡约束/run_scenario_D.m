function run_scenario_D(strategies, Scenarios_AC_Up, Scenarios_EV_Up, Scenarios_AC_Down, Scenarios_EV_Down, direction_signal)
    fprintf('\n>>> 场景 D: 鲁棒性测试 <<<\n');
    T_steps = length(direction_signal);
    
    if ~isempty(strategies{1}) && ~isempty(strategies{3})
        [~, worst_idx] = min(sum(Scenarios_AC_Up + Scenarios_EV_Up));
        fprintf('  - 最恶劣场景: #%d\n', worst_idx);
        
        Real_Cap_AC = zeros(T_steps, 1); Real_Cap_EV = zeros(T_steps, 1);
        for t=1:T_steps
            if direction_signal(t)==1, Real_Cap_AC(t)=Scenarios_AC_Up(t,worst_idx); Real_Cap_EV(t)=Scenarios_EV_Up(t,worst_idx);
            else, Real_Cap_AC(t)=abs(Scenarios_AC_Down(t,worst_idx)); Real_Cap_EV(t)=abs(Scenarios_EV_Down(t,worst_idx)); end
        end
        
        calc_viol = @(P, Cap) sum(max(0, P - Cap));
        v_neu = calc_viol(strategies{1}.P_AC, Real_Cap_AC) + calc_viol(strategies{1}.P_EV, Real_Cap_EV);
        v_rob = calc_viol(strategies{3}.P_AC, Real_Cap_AC) + calc_viol(strategies{3}.P_EV, Real_Cap_EV);
        
        fprintf('  - 中性策略(Beta=0) 违约量: %.2f MW\n', v_neu);
        fprintf('  - 规避策略(Beta=10) 违约量: %.2f MW\n', v_rob);
        
        figure('Color','w', 'Position', [600, 300, 500, 400]); 
        b = bar([v_neu, v_rob], 0.5);
        b.FaceColor = 'flat'; b.CData(1,:) = [0.8 0.2 0.2]; b.CData(2,:) = [0.2 0.6 0.2];
        set(gca, 'XTickLabel', {'中性 (\beta=0)', '规避 (\beta=10)'});
        ylabel('极端场景实际违约量 (MW)');
        grid on;
        print(gcf, '极端场景鲁棒性测试.png', '-dpng', '-r300');
    end
end
