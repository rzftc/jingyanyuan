function run_scenario_E_expand(P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
    Physical_AC_Up, Physical_EV_Up, R_Gen_Max, R_Shed_Max, ...
    cost_params, net_params, direction_signal, dt, options, N_scenarios, N_line, N_bus)

    fprintf('\n>>> 场景 E: 置信水平对经济性的影响测试 (Reliability vs Cost) <<<\n');
    confs = [0.85, 0.90, 0.95, 0.99];
    num_conf = length(confs);
    e_total_cost = zeros(1, num_conf);
    e_slack_sum = zeros(1, num_conf);
    T_steps = length(P_grid_demand);

    risk_E.beta = 1;
    % risk_E.rho_pen = 10000;
    risk_E.rho_pen = 1000; 
    risk_E.tight_factor = 0.9;
    for k = 1:num_conf
        curr_alpha = confs(k);
        risk_E.confidence = curr_alpha;
        fprintf('  [测试 %d/%d] 置信水平 alpha = %.2f ... ', k, num_conf, curr_alpha);
        
        net_params_safe = net_params;
        net_params_safe.ShedDist = zeros(N_bus, 1);

        [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast_ramp(...
                P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
                Physical_AC_Up, Physical_EV_Up, R_Gen_Max, R_Shed_Max, ...
                cost_params, risk_E, net_params_safe);
                
        start_row_net = 2 * N_scenarios; 
        % 位于 run_scenario_B.m 的循环内
        for t = 1:T_steps
            if direction_signal(t) == 1
                rows_t = start_row_net + (t-1)*2*N_line + (1 : 2*N_line);

                % 1. 翻转 AC 和 EV（代表增加负荷/吸电）—— 原有逻辑，正确
                A(rows_t, info.idx_P_AC(t)) = -A(rows_t, info.idx_P_AC(t));
                A(rows_t, info.idx_P_EV(t)) = -A(rows_t, info.idx_P_EV(t));

                % 2. 【新增】翻转 Gen（代表减少发电/等效吸电）—— 修复逻辑
                % 只有翻转后，P_Gen > 0 才代表“火电出力下降”，对潮流表现为负贡献
                A(rows_t, info.idx_P_Gen(t)) = -A(rows_t, info.idx_P_Gen(t));
            end
        end
        
        idx_pow = [info.idx_P_AC, info.idx_P_EV, info.idx_P_Gen, info.idx_P_Shed];
        for idx = idx_pow, H(idx, idx) = H(idx, idx) * dt; end
        f(idx_pow) = f(idx_pow) * dt;
        
        [x_opt, ~, exitflag] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
        
        if exitflag > 0
            P_AC_curr = x_opt(info.idx_P_AC); P_EV_curr = x_opt(info.idx_P_EV);
            P_Gen_curr = x_opt(info.idx_P_Gen); P_Shed_curr = x_opt(info.idx_P_Shed);
            
            real_cost = sum((cost_params.c1_ac*P_AC_curr + cost_params.c2_ac*P_AC_curr.^2)*dt + ...
                            (cost_params.c1_ev*P_EV_curr + cost_params.c2_ev*P_EV_curr.^2)*dt + ...
                            (cost_params.c1_gen*P_Gen_curr + cost_params.c2_gen*P_Gen_curr.^2)*dt + ...
                            (cost_params.c1_shed*P_Shed_curr)*dt);
            e_total_cost(k) = real_cost;
            e_slack_sum(k) = sum(P_Gen_curr + P_Shed_curr) * dt;
            fprintf('成功。总成本: %.2f 元, 昂贵资源调用: %.2f MWh\n', real_cost, e_slack_sum(k));
        else
            fprintf('失败。\n');
        end
    end

    if any(e_total_cost > 0)
        fig_conf = figure('Name', '场景E_置信度影响', 'Color', 'w', 'Position', [600, 100, 700, 450]);
        yyaxis left;
        b = bar(categorical(confs), e_total_cost, 0.5, 'FaceColor', [0.2 0.6 0.8]);
        ylabel('系统总运行成本 (元)', 'FontSize', 12, 'FontName', 'Microsoft YaHei');
        ylim([min(e_total_cost)*0.9, max(e_total_cost)*1.1]);
        for i = 1:length(e_total_cost)
            text(i, e_total_cost(i), sprintf('%.0f', e_total_cost(i)), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                 'FontSize', 10, 'Color', 'b');
        end
        yyaxis right;
        plot(1:num_conf, e_slack_sum, 'r-^', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        ylabel('火电与切负荷调用量 (MWh)', 'FontSize', 12, 'FontName', 'Microsoft YaHei');
        ax = gca; ax.YColor = 'r';
        xlabel('置信水平 \alpha (可靠性要求)', 'FontSize', 12, 'FontName', 'Microsoft YaHei');
        legend({'运行成本', '备用资源调用'}, 'Location', 'northwest', 'FontSize', 11);
        grid on;
        print(fig_conf, '置信水平经济性分析.png', '-dpng', '-r300');
    end
end