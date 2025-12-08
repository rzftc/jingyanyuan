function strategies = run_scenario_B(beta_values, Max_Iter, N_scenarios, N_bus, N_line, dt, ...
    P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, Physical_AC_Up, Physical_EV_Up, ...
    R_Gen_Max, R_Shed_Max, cost_params, net_params, direction_signal, lambda_SDCI, lambda_Rho, options)

    fprintf('\n>>> 场景 B: 风险偏好灵敏度分析 <<<\n');
    b_run_cost = nan(1, length(beta_values)); 
    b_slack_sum = nan(1, length(beta_values));
    b_risk_val = nan(1, length(beta_values));
    strategies = cell(1, length(beta_values));
    T_steps = length(P_grid_demand);

    for i = 1:length(beta_values)
        beta = beta_values(i);
        fprintf('  工况 %d (Beta=%d): \n', i, beta); 
        
        risk_p.beta = beta;
        risk_p.confidence = 0.95;
        risk_p.rho_pen = 300; 
        
        P_AC_prev = zeros(T_steps, 1); P_EV_prev = zeros(T_steps, 1);
        
        for iter = 1:Max_Iter
            net_params_safe = net_params;
            net_params_safe.ShedDist = zeros(N_bus, 1);

            [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast(...
                P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
                Physical_AC_Up, Physical_EV_Up, ...
                R_Gen_Max, R_Shed_Max, ...
                cost_params, risk_p, net_params_safe);
            
            start_row_net = 2 * N_scenarios; 
            for t = 1:T_steps
                if direction_signal(t) == 1 
                    rows_t = start_row_net + (t-1)*2*N_line + (1 : 2*N_line);
                    A(rows_t, info.idx_P_AC(t)) = -A(rows_t, info.idx_P_AC(t));
                    A(rows_t, info.idx_P_EV(t)) = -A(rows_t, info.idx_P_EV(t));
                end
            end
            
            idx_pow = [info.idx_P_AC, info.idx_P_EV, info.idx_P_Gen, info.idx_P_Shed];
            for idx = idx_pow, H(idx, idx) = H(idx, idx) * dt; end
            f(idx_pow) = f(idx_pow) * dt;
            
            if iter > 1
                f(info.idx_P_AC) = f(info.idx_P_AC) + (lambda_SDCI * P_EV_prev * dt);
                f(info.idx_P_EV) = f(info.idx_P_EV) + (lambda_SDCI * P_AC_prev * dt);
            end
            
            [x_opt, ~, exitflag] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
            
            if exitflag > 0
                P_AC_curr = x_opt(info.idx_P_AC);
                P_EV_curr = x_opt(info.idx_P_EV);
                P_Gen_curr = x_opt(info.idx_P_Gen);
                P_Shed_curr = x_opt(info.idx_P_Shed);
                eta_val = x_opt(info.idx_eta);
                z_val = x_opt(info.idx_z);
                
                P_AC_prev = P_AC_curr; P_EV_prev = P_EV_curr;
                
                strategies{i}.P_AC = P_AC_curr; strategies{i}.P_EV = P_EV_curr;
                strategies{i}.P_Shed = P_Shed_curr; strategies{i}.P_Gen = P_Gen_curr;
                
                % SDCI/Rho 记录
                n_dummy = ones(T_steps, 1);
                val_SDCI = calculate_SDCI_local(n_dummy, n_dummy, P_AC_curr, P_EV_curr);
                val_Rho  = calculate_Rho_local(n_dummy, P_AC_curr, n_dummy, P_EV_curr);
                if iter == 1
                    strategies{i}.SDCI_History = zeros(Max_Iter, 1);
                    strategies{i}.Rho_History = zeros(Max_Iter, 1);
                end
                strategies{i}.SDCI_History(iter) = val_SDCI;
                strategies{i}.Rho_History(iter) = val_Rho;

                cost_gen = sum((cost_params.c1_ac*P_AC_curr + cost_params.c2_ac*P_AC_curr.^2)*dt + ...
                               (cost_params.c1_ev*P_EV_curr + cost_params.c2_ev*P_EV_curr.^2)*dt + ...
                               (cost_params.c1_gen*P_Gen_curr + cost_params.c2_gen*P_Gen_curr.^2)*dt);
                cost_slack = sum(cost_params.c1_shed * P_Shed_curr * dt); 
                total_real_cost = cost_gen + cost_slack;
                
                cvar_val = eta_val + (1 / (N_scenarios * (1 - risk_p.confidence))) * sum(z_val);
                b_run_cost(i) = cost_gen;
                b_slack_sum(i) = sum(P_Shed_curr + P_Gen_curr) * dt;
                b_risk_val(i) = cvar_val;

                if iter == Max_Iter
                    fprintf('    发电成本: %.2f (元), 切负荷量: %.2f (MWh), 总成本: %.2f (元), CVaR风险: %.2f (MW), rho: %.4f, sdci: %.4f\n', ...
                        cost_gen, b_slack_sum(i), total_real_cost, cvar_val, val_Rho, val_SDCI);
                end
            else
                fprintf('    失败 (Exitflag %d)\n', exitflag); break;
            end
        end
    end

    % --- 绘图 B ---
    if any(~isnan(b_slack_sum))
        figure('Name', '场景B_风险灵敏度', 'Color', 'w', 'Position', [100, 100, 900, 400]);
        yyaxis left; 
        b = bar(1:3, b_slack_sum, 0.5, 'FaceColor', [0.8 0.3 0.3]); 
        ylabel('切负荷量 + 火电调度量 (MWh)'); 
        set(gca, 'XTick', 1:3, 'XTickLabel', beta_values);
        for i = 1:length(b_slack_sum)
            text(i, b_slack_sum(i), sprintf('%.2f', b_slack_sum(i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'FontSize', 12, 'Color', [0.6 0.1 0.1], 'FontWeight', 'bold');
        end
        
        yyaxis right; 
        plot(1:3, b_risk_val, 'b-o', 'LineWidth', 2, 'MarkerSize', 8); 
        ylabel('CVaR 潜在违约风险 (MW)'); 
        for i = 1:length(b_risk_val)
            text(i, b_risk_val(i), sprintf('%.2f', b_risk_val(i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'FontSize', 12, 'Color', 'b', 'FontWeight', 'bold');
        end
        xlabel('风险厌恶系数 \beta');
        legend('切负荷 + 火电 (安全性)', '潜在违约风险 (经济性)', 'Location', 'best');
        grid on;
        print(gcf, '风险偏好灵敏度分析.png', '-dpng', '-r300');
    end
end

