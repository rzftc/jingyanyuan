function run_scenario_I_lambda_sensitivity(beta_fixed, lambda_values, Max_Iter, ...
    N_scenarios, N_bus, N_line, dt, ...
    P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
    Physical_AC_Up, Physical_EV_Up, ...
    R_Gen_Max, R_Shed_Max, cost_params, net_params, direction_signal, options)

    fprintf('\n>>> 场景 I: 协同惩罚权重灵敏度分析 (分图优化版) <<<\n');
    
    % --- 初始化结果存储 ---
    num_lam = length(lambda_values);
    res_cost_total = nan(1, num_lam);
    res_cvar_val   = nan(1, num_lam); 
    res_sdci       = nan(1, num_lam);
    res_rho        = nan(1, num_lam);
    
    T_steps = length(P_grid_demand);
    
    % --- 循环遍历 Lambda ---
    for i = 1:num_lam
        lam = lambda_values(i);
        fprintf('  权重工况 %d (Lambda=%.1f): \n', i, lam);
        
        % 构造风险参数 (固定 Beta)
        risk_p.beta = beta_fixed;
        risk_p.confidence = 0.95;
        risk_p.rho_pen = 5000; 
        risk_p.tight_factor = 0.9;
        
        % 上一轮结果初始化
        P_AC_prev = zeros(T_steps, 1); 
        P_EV_prev = zeros(T_steps, 1);
        
        % --- 线性化迭代 ---
        for iter = 1:Max_Iter
            % 1. 构造基础 QP
            net_params_safe = net_params;
            net_params_safe.ShedDist = zeros(N_bus, 1); 
            
            [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast_ramp_tly(...
                P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
                Physical_AC_Up, Physical_EV_Up, ...
                R_Gen_Max, R_Shed_Max, ...
                cost_params, risk_p, net_params_safe);
            
            % 2. 符号翻转逻辑
            start_row_net = (T_steps + 1) * N_scenarios;
            for t = 1:T_steps
                if direction_signal(t) == 1
                    rows_t = start_row_net + (t-1)*2*N_line + (1 : 2*N_line);
                    A(rows_t, info.idx_P_AC(t)) = -A(rows_t, info.idx_P_AC(t));
                    A(rows_t, info.idx_P_EV(t)) = -A(rows_t, info.idx_P_EV(t));
                    A(rows_t, info.idx_P_Gen(t)) = -A(rows_t, info.idx_P_Gen(t));
                end
            end
            
            % 3. dt 缩放
            idx_pow = [info.idx_P_AC, info.idx_P_EV, info.idx_P_Gen, info.idx_P_Shed];
            for idx = idx_pow, H(idx, idx) = H(idx, idx) * dt; end
            f(idx_pow) = f(idx_pow) * dt;
            
            % 4. 引入协同惩罚项
            if iter > 1
                Mean_AC = mean(P_AC_prev);
                Mean_EV = mean(P_EV_prev);
                
                f(info.idx_P_AC) = f(info.idx_P_AC) + (lam * P_EV_prev * dt);
                f(info.idx_P_EV) = f(info.idx_P_EV) + (lam * P_AC_prev * dt);
                
                f(info.idx_P_AC) = f(info.idx_P_AC) + (lam * (P_EV_prev - Mean_EV) * dt);
                f(info.idx_P_EV) = f(info.idx_P_EV) + (lam * (P_AC_prev - Mean_AC) * dt);
            end
            
            % 5. Cplex 求解
            cplex = Cplex('ScenarioI');
            cplex.Model.sense = 'minimize';
            cplex.Model.Q = H;
            cplex.Model.obj = f;
            cplex.Model.lb = lb;
            cplex.Model.ub = ub;
            
            if isempty(A)
                A_combined = Aeq; lhs = beq; rhs = beq;
            elseif isempty(Aeq)
                A_combined = A; lhs = -inf(size(b)); rhs = b;
            else
                A_combined = [A; Aeq];
                lhs = [-inf(size(b)); beq];
                rhs = [b; beq];
            end
            cplex.Model.A = A_combined;
            cplex.Model.lhs = lhs;
            cplex.Model.rhs = rhs;
            cplex.DisplayFunc = []; 
            
            cplex.solve();
            
            % 6. 结果提取
            if isfield(cplex.Solution, 'x') && ~isempty(cplex.Solution.x)
                x_opt = cplex.Solution.x;
                status = cplex.Solution.status;
                if status == 1 || status == 101 || status == 102
                    exitflag = 1;
                else
                    exitflag = -2;
                end
            else
                exitflag = -2;
            end
            
            if exitflag > 0
                P_AC_curr = x_opt(info.idx_P_AC);
                P_EV_curr = x_opt(info.idx_P_EV);
                P_Gen_curr = x_opt(info.idx_P_Gen);
                P_Shed_curr = x_opt(info.idx_P_Shed);
                eta_val = x_opt(info.idx_eta);
                z_val = x_opt(info.idx_z);
                
                P_AC_prev = P_AC_curr;
                P_EV_prev = P_EV_curr;
                
                cost_gen = sum((cost_params.c1_ac.*P_AC_curr + cost_params.c2_ac*P_AC_curr.^2)*dt + ...
                               (cost_params.c1_ev.*P_EV_curr + cost_params.c2_ev*P_EV_curr.^2)*dt + ...
                               (cost_params.c1_gen.*P_Gen_curr + cost_params.c2_gen*P_Gen_curr.^2)*dt);
                cost_slack = sum(cost_params.c1_shed * P_Shed_curr * dt);
                
                cvar_val = eta_val + (1 / (N_scenarios * (1 - risk_p.confidence))) * sum(z_val);
                risk_cost_val = beta_fixed * cvar_val;
                
                n_dummy = ones(T_steps, 1);
                val_SDCI = calculate_SDCI_local(n_dummy, n_dummy, P_AC_curr, P_EV_curr);
                val_Rho  = calculate_Rho_local(n_dummy, P_AC_curr, n_dummy, P_EV_curr);
                
                if iter == Max_Iter
                    res_cvar_val(i)   = cvar_val;
                    res_cost_total(i) = cost_gen + cost_slack + risk_cost_val;
                    res_sdci(i)       = val_SDCI;
                    res_rho(i)        = val_Rho;
                    
                    fprintf('    > 结果: 总成本=%.2f, CVaR=%.2f, SDCI=%.4f, Rho=%.4f\n', ...
                        res_cost_total(i), cvar_val, val_SDCI, val_Rho);
                end
            else
                warning('    求解失败 (Exitflag %d)', exitflag);
                break;
            end
        end
    end
    
    %% ================= 绘图 I-1: 成本与CVaR (上下分图版) =================
    % 解决量级差异过大导致趋势看不清的问题
    if any(~isnan(res_cost_total))
        % 增加画布高度，容纳两个子图
        figure('Name', '场景I_权重灵敏度_成本_分图', 'Color', 'w', 'Position', [100, 100, 800, 800]);
        
        % --- 上图：系统总成本 ---
        subplot(2, 1, 1); % 上半部分
        plot(1:num_lam, res_cost_total, '-o', 'LineWidth', 2.0, 'MarkerSize', 8, ...
             'MarkerFaceColor', [0.2 0.6 0.8], 'Color', [0.2 0.6 0.8]);
        
        ylabel('系统总成本 (元)', 'FontSize', 20, 'FontWeight', 'bold');
        grid on;
        set(gca, 'XTick', 1:num_lam, 'XTickLabel', lambda_values, 'FontSize', 18, 'LineWidth', 1.5);
        % axis tight; % 自动紧凑坐标轴，凸显变化
        
        % --- 下图：CVaR 值 ---
        subplot(2, 1, 2); % 下半部分
        plot(1:num_lam, res_cvar_val, '-^', 'LineWidth', 2.0, 'MarkerSize', 8, ...
             'MarkerFaceColor', [0.8 0.3 0.3], 'Color', [0.8 0.3 0.3]);
        
        ylabel('CVaR 风险值', 'FontSize', 20, 'FontWeight', 'bold');
        xlabel('协同惩罚权重 \lambda', 'FontSize', 20, 'FontWeight', 'bold');
        grid on;
        set(gca, 'XTick', 1:num_lam, 'XTickLabel', lambda_values, 'FontSize', 18, 'LineWidth', 1.5);
        % axis tight; % 关键：自动适应 CVaR 的微小变化范围
        
        % 保存
        print(gcf, '协同权重灵敏度_成本_分图.png', '-dpng', '-r300');
    end

    %% ================= 绘图 I-2: 指标灵敏度 (保持双Y轴) =================
    if any(~isnan(res_sdci))
        figure('Name', '场景I_权重灵敏度_指标', 'Color', 'w', 'Position', [100, 650, 900, 500]);
        
        yyaxis left
        plot(1:num_lam, res_sdci, '-s', 'LineWidth', 2.5, 'MarkerSize', 10, ...
            'MarkerFaceColor', [0 0.4470 0.7410], 'Color', [0 0.4470 0.7410]);
        ylabel('互补性指数 (SDCI)', 'FontSize', 20, 'FontWeight', 'bold');
        set(gca, 'YColor', [0 0.4470 0.7410]);
        
        yyaxis right
        plot(1:num_lam, res_rho, '-^', 'LineWidth', 2.5, 'MarkerSize', 10, ...
            'MarkerFaceColor', [0.8500 0.3250 0.0980], 'Color', [0.8500 0.3250 0.0980]);
        ylabel('相关系数 (Rho)', 'FontSize', 20, 'FontWeight', 'bold');
        set(gca, 'YColor', [0.8500 0.3250 0.0980]);
        
        set(gca, 'XTick', 1:num_lam, 'XTickLabel', lambda_values, 'FontSize', 18, 'LineWidth', 1.5);
        xlabel('协同惩罚权重 \lambda', 'FontSize', 20, 'FontWeight', 'bold');
        grid on;
        
        yyaxis left;  ylim_l = ylim; ylim([ylim_l(1), ylim_l(2)*1.2]);
        yyaxis right; ylim_r = ylim; ylim([ylim_r(1), ylim_r(2)*1.2]);
        
        legend({'SDCI', 'Rho'}, ...
            'Location', 'North', 'Orientation', 'horizontal', 'FontSize', 16);
            
        print(gcf, '协同权重灵敏度_指标.png', '-dpng', '-r300');
    end
end