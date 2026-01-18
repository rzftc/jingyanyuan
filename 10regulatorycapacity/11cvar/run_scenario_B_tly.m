function strategies = run_scenario_B_tly(beta_values, Max_Iter, N_scenarios, N_bus, N_line, dt, ...
    P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, Physical_AC_Up, Physical_EV_Up, ...
    R_Gen_Max, R_Shed_Max, cost_params, net_params, direction_signal, lambda_SDCI, lambda_Rho, options)
    fprintf('\n>>> 场景 B: 风险偏好灵敏度分析 <<<\n');
    b_run_cost = nan(1, length(beta_values)); 
    b_slack_sum = nan(1, length(beta_values));
    b_agg_sum = nan(1, length(beta_values)); % [新增] 记录聚合体(AC+EV)总调度量
    b_risk_val = nan(1, length(beta_values));
    strategies = cell(1, length(beta_values));
    T_steps = length(P_grid_demand);
    for i = 1:length(beta_values)
        beta = beta_values(i);
        fprintf('  工况 %d (Beta=%d): \n', i, beta); 
        
        risk_p.beta = beta;
        risk_p.confidence = 0.95;
        risk_p.rho_pen = 5000; 
        risk_p.tight_factor = 0.9;
        
        P_AC_prev = zeros(T_steps, 1); P_EV_prev = zeros(T_steps, 1);
        
        for iter = 1:Max_Iter
            net_params_safe = net_params;
            net_params_safe.ShedDist = zeros(N_bus, 1);
            [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast_ramp_tly(...
                P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
                Physical_AC_Up, Physical_EV_Up, ...
                R_Gen_Max, R_Shed_Max, ...
                cost_params, risk_p, net_params_safe);
            
            % start_row_net = 2 * N_scenarios; 
            start_row_net = (T_steps + 1) * N_scenarios;
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
            
            % --- 关键修改：顺序线性化惩罚迭代 (包含 SDCI 和 Rho) ---
            if iter > 1
                % 1. 计算上一轮的均值 (用于 Rho 的协方差去均值处理)
                Mean_AC = mean(P_AC_prev);
                Mean_EV = mean(P_EV_prev);
                
                % 2. SDCI 惩罚项 (减少同向重叠)
                % 原理：最小化 P_AC * P_EV，线性化为 P_AC * P_EV_prev
                f(info.idx_P_AC) = f(info.idx_P_AC) + (lambda_SDCI * P_EV_prev * dt);
                f(info.idx_P_EV) = f(info.idx_P_EV) + (lambda_SDCI * P_AC_prev * dt);
                
                % 3. Rho 惩罚项 (减少同向趋势 / 相关性)
                % 原理：最小化 Cov(AC, EV，线性化为 P_AC * (P_EV_prev - Mean_EV)
                % 对应论文公式 (4-36)
                f(info.idx_P_AC) = f(info.idx_P_AC) + (lambda_Rho * (P_EV_prev - Mean_EV) * dt);
                f(info.idx_P_EV) = f(info.idx_P_EV) + (lambda_Rho * (P_AC_prev - Mean_AC) * dt);
            end
            
            
         % [x_opt, ~, exitflag] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
         [x_opt, ~, exitflag] = cplexqp(H, f, A, b, Aeq, beq, lb, ub, [], options);
           
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
                cost_gen = sum((cost_params.c1_ac.*P_AC_curr + cost_params.c2_ac*P_AC_curr.^2)*dt + ...
                               (cost_params.c1_ev.*P_EV_curr + cost_params.c2_ev*P_EV_curr.^2)*dt + ...
                               (cost_params.c1_gen.*P_Gen_curr + cost_params.c2_gen*P_Gen_curr.^2)*dt);
                cost_slack = sum(cost_params.c1_shed * P_Shed_curr * dt); 
                total_real_cost = cost_gen + cost_slack;
                
                cvar_val = eta_val + (1 / (N_scenarios * (1 - risk_p.confidence))) * sum(z_val);
                b_run_cost(i) = cost_gen;
                b_slack_sum(i) = sum(P_Shed_curr + P_Gen_curr) * dt;
                b_agg_sum(i)   = sum(P_AC_curr + P_EV_curr) * dt; % [新增] 聚合体总能量
                b_risk_val(i) = cvar_val;
                if iter == Max_Iter
                    fprintf('    发电成本: %.2f (元), 火电调节量: %.2f (MWh), 总成本: %.2f (元), CVaR风险: %.2f (MW), rho: %.4f, sdci: %.4f\n', ...
                        cost_gen, b_slack_sum(i), total_real_cost, cvar_val, val_Rho, val_SDCI);
                    
                    % --- [新增] 目标函数成分量级诊断 (用于验证强凸主导假设) ---
                    Term_Cost = cost_gen + cost_slack;
                    Term_Risk = beta * cvar_val; 
                    Term_SDCI = lambda_SDCI * val_SDCI;
                    Term_Rho  = lambda_Rho * val_Rho;
                    
                    Total_Obj_Approx = Term_Cost + Term_Risk + Term_SDCI + Term_Rho;
                    if Total_Obj_Approx == 0, Total_Obj_Approx = 1; end % 防止分母为0

                    Ratio_Cost = Term_Cost / Total_Obj_Approx * 100;
                    Ratio_Risk = Term_Risk / Total_Obj_Approx * 100;
                    Ratio_Penalty = (Term_SDCI + Term_Rho) / Total_Obj_Approx * 100;

                    fprintf('    [量级诊断] ------------------------------------------------------\n');
                    fprintf('    1. 基础成本 (Quadratic Cost):  %12.2f  (占比 %5.1f%%)\n', Term_Cost, Ratio_Cost);
                    fprintf('    2. 风险成本 (CVaR Risk):       %12.2f  (占比 %5.1f%%)\n', Term_Risk, Ratio_Risk);
                    fprintf('    3. 协同惩罚 (SDCI + Rho):      %12.2f  (占比 %5.1f%%)\n', Term_SDCI + Term_Rho, Ratio_Penalty);
                    fprintf('    -----------------------------------------------------------------\n');
                end
            else
                fprintf('    失败 (Exitflag %d)\n', exitflag); break;
            end
        end
    end
    % --- 绘图 B: 风险灵敏度 (修改版: 格式优化) ---
    if any(~isnan(b_slack_sum))
        figure('Name', '场景B_风险灵敏度', 'Color', 'w', 'Position', [100, 100, 900, 500]);
        yyaxis left; 
        
        % 准备数据：第一列是聚合体(AC+EV)，第二列是备用(Gen+Shed)
        data_to_plot = [b_agg_sum(:), b_slack_sum(:)];
        b = bar(1:length(beta_values), data_to_plot, 'grouped');
        
        % 设置颜色区分
        b(1).FaceColor = [0.2 0.6 0.8]; % 蓝色系：聚合体
        b(2).FaceColor = [0.8 0.3 0.3]; % 红色系：备用资源
        
        % [修改] 坐标轴标签大字号加粗
        ylabel('调度能量 (MWh)', 'FontSize', 20, 'FontWeight', 'bold'); 
        % [修改] 坐标轴刻度大字号，框线加粗
        set(gca, 'XTick', 1:length(beta_values), 'XTickLabel', beta_values, 'FontSize', 18, 'LineWidth', 1.5);
        
        % 为分组柱状图添加数值标签
        for k = 1:2
            xt = b(k).XEndPoints;
            yt = b(k).YEndPoints;
            for i = 1:length(xt)
                % [修改] 数值标签字号调整为 14
                text(xt(i), yt(i), sprintf('%.2f', yt(i)), ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                    'FontSize', 14, 'Color', b(k).FaceColor * 0.8, 'FontWeight', 'bold');
            end
        end
        
        yyaxis right; 
        plot(1:length(beta_values), b_risk_val, 'k-o', 'LineWidth', 2.0, 'MarkerSize', 8, 'MarkerFaceColor', 'y'); 
        % [修改] 坐标轴标签大字号加粗
        ylabel('CVaR 潜在违约风险/元', 'FontSize', 20, 'FontWeight', 'bold'); 
        for i = 1:length(b_risk_val)
            % [修改] 数值标签字号调整为 16
            text(i, b_risk_val(i), sprintf('%.2f', b_risk_val(i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'FontSize', 16, 'Color', 'k', 'FontWeight', 'bold');
        end
        
        % [修改] 横坐标标签大字号加粗
        xlabel('风险厌恶系数 \beta', 'FontSize', 20, 'FontWeight', 'bold');
        % [修改] 图例大字号
        legend({'聚合体 (AC+EV)', '备用', '潜在违约风险 (CVaR)'}, ...
       'Location', 'SouthOutside', 'Orientation', 'horizontal', 'FontSize', 16);
        
        % [修改] 去除网格线
        % grid on; % 已注释
        print(gcf, '风险偏好灵敏度分析.png', '-dpng', '-r300');
    end
end