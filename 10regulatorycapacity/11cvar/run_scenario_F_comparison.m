function run_scenario_F_comparison(beta_val, Max_Iter, N_scenarios, N_bus, N_line, dt, ...
    P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
    Reliable_AC_Up, Reliable_EV_Up, ... 
    R_Gen_Max, R_Shed_Max, cost_params, net_params, direction_signal, lambda_SDCI, lambda_Rho, options)

    fprintf('\n==========================================================\n');
    fprintf('>>> 场景 F: 协同约束效益对比分析 (含风险成本验证) <<<\n');
    fprintf('==========================================================\n');
    
    % 定义两种情况
    cases = {'无协同约束 (Baseline)', '有协同约束 (Optimized)'};
    lambdas_S = [0, lambda_SDCI]; % SDCI 权重
    lambdas_R = [0, lambda_Rho];  % Rho 权重
    
    results = struct();
    T_steps = length(P_grid_demand);
    
    % 统一风险参数
    risk_p.beta = beta_val;
    risk_p.confidence = 0.95;
    risk_p.rho_pen = 5000; 
    risk_p.tight_factor = 0.95; 
    
    % 循环运行两种情况
    for k = 1:2
        fprintf('  正在运行: %s ...\n', cases{k});
        
        l_sdci = lambdas_S(k);
        l_rho  = lambdas_R(k);
        
        P_AC_prev = zeros(T_steps, 1); 
        P_EV_prev = zeros(T_steps, 1);
        
        % 迭代求解
        for iter = 1:Max_Iter
            % 1. 构建 QP 问题
            [H, f, A, b, Aeq, beq, lb, ub, info] = construct_risk_constrained_qp_fast_ramp_tly(...
                P_grid_demand, Scenarios_AC_Up, Scenarios_EV_Up, ...
                Reliable_AC_Up, Reliable_EV_Up, ... 
                R_Gen_Max, R_Shed_Max, ...
                cost_params, risk_p, net_params);
            
            % 2. 处理方向信号
            % start_row_net = 2 * N_scenarios; 
            start_row_net = (T_steps + 1) * N_scenarios;
            for t = 1:T_steps
                if direction_signal(t) == 1 % Up Regulation
                    rows_t = start_row_net + (t-1)*2*N_line + (1 : 2*N_line);
                    A(rows_t, info.idx_P_AC(t)) = -A(rows_t, info.idx_P_AC(t));
                    A(rows_t, info.idx_P_EV(t)) = -A(rows_t, info.idx_P_EV(t));
                    A(rows_t, info.idx_P_Gen(t)) = -A(rows_t, info.idx_P_Gen(t));
                end
            end
            
            % 3. 时间步长修正
            idx_pow = [info.idx_P_AC, info.idx_P_EV, info.idx_P_Gen, info.idx_P_Shed];
            for idx = idx_pow, H(idx, idx) = H(idx, idx) * dt; end
            f(idx_pow) = f(idx_pow) * dt;
            
            % 4. 添加协同惩罚 (仅在迭代 > 1 时)
            if iter > 1 && (l_sdci > 0 || l_rho > 0)
                Mean_AC = mean(P_AC_prev);
                Mean_EV = mean(P_EV_prev);
                
                f(info.idx_P_AC) = f(info.idx_P_AC) + (l_sdci * P_EV_prev * dt);
                f(info.idx_P_EV) = f(info.idx_P_EV) + (l_sdci * P_AC_prev * dt);
                
                f(info.idx_P_AC) = f(info.idx_P_AC) + (l_rho * (P_EV_prev - Mean_EV) * dt);
                f(info.idx_P_EV) = f(info.idx_P_EV) + (l_rho * (P_AC_prev - Mean_AC) * dt);
            end
            
            % 5. 求解 (修改部分：使用 Cplex 类对象)
            % -----------------------------------------------------------
            cplex = Cplex('ScenarioF');
            cplex.Model.sense = 'minimize';
            
            cplex.Model.Q = H;
            cplex.Model.obj = f;
            cplex.Model.lb = lb;
            cplex.Model.ub = ub;
            
            % 处理约束 A*x <= b 和 Aeq*x = beq
            if isempty(A)
                A_combined = Aeq;
                lhs = beq;
                rhs = beq;
            elseif isempty(Aeq)
                A_combined = A;
                lhs = -inf(size(b));
                rhs = b;
            else
                A_combined = [A; Aeq];
                lhs = [-inf(size(b)); beq];
                rhs = [b; beq];
            end
            cplex.Model.A = A_combined;
            cplex.Model.lhs = lhs;
            cplex.Model.rhs = rhs;
            
            cplex.DisplayFunc = []; % 关闭输出
            
            cplex.solve();
            
            if isfield(cplex.Solution, 'x') && ~isempty(cplex.Solution.x)
                x_opt = cplex.Solution.x;
                status = cplex.Solution.status;
                if status == 1 || status == 101 || status == 102
                    exitflag = 1;
                else
                    exitflag = -2;
                end
            else
                x_opt = [];
                exitflag = -2;
            end
            % -----------------------------------------------------------

            if exitflag > 0
                P_AC_curr = x_opt(info.idx_P_AC);
                P_EV_curr = x_opt(info.idx_P_EV);
                P_AC_prev = P_AC_curr; 
                P_EV_prev = P_EV_curr;
                
                % --- 新增：计算 CVaR 风险值 ---
                eta_val = x_opt(info.idx_eta);
                z_val   = x_opt(info.idx_z);
                % CVaR = VaR + 尾部期望损失
                current_cvar = eta_val + (1 / (N_scenarios * (1 - risk_p.confidence))) * sum(z_val);
            else
                warning('场景 F 在 %s 迭代 %d 求解失败', cases{k}, iter);
                current_cvar = NaN;
            end
        end
        
        % 存储结果
        results(k).P_AC = P_AC_curr;
        results(k).P_EV = P_EV_curr;
        results(k).P_Total = P_AC_curr + P_EV_curr;
        results(k).CVaR = current_cvar; % 记录 CVaR
        
        results(k).Max_Peak = max(results(k).P_Total); 
        results(k).Std_Dev  = std(results(k).P_Total); 
        
        dummy = ones(T_steps, 1);
        results(k).Val_SDCI = calculate_SDCI_local(dummy, dummy, P_AC_curr, P_EV_curr);
    end
    
    %% ================= 结果对比输出 =================
    fprintf('\n----------------------------------------------------------\n');
    fprintf('指标对比 \t\t| 无约束 (Baseline) \t| 有约束 (Optimized) \t| 改善/变化率\n');
    fprintf('----------------------------------------------------------\n');
    
    imp_peak = (results(2).Max_Peak - results(1).Max_Peak) / results(1).Max_Peak * 100; % 变化率
    imp_std  = (results(1).Std_Dev - results(2).Std_Dev) / results(1).Std_Dev * 100;    % 改善率(下降为优)
    imp_cvar = (results(1).CVaR - results(2).CVaR) / results(1).CVaR * 100;             % 改善率(下降为优)
    
    fprintf('聚合峰值 (MW) \t| %.4f \t\t| %.4f \t\t| %+.2f%% (潜力释放)\n', ...
        results(1).Max_Peak, results(2).Max_Peak, imp_peak);
    fprintf('波动标准差 \t| %.4f \t\t| %.4f \t\t| -%.2f%% (平抑波动)\n', ...
        results(1).Std_Dev, results(2).Std_Dev, imp_std);
    fprintf('CVaR 风险 (MW)\t| %.4f \t\t| %.4f \t\t| -%.2f%% (风险降低)\n', ...
        results(1).CVaR, results(2).CVaR, imp_cvar);
    fprintf('----------------------------------------------------------\n');
    
    %% ================= 独立绘图与保存 =================
    % 通用绘图设置
    font_name = 'Microsoft YaHei'; 
    font_size = 14; % 修改：从12增大到14
    
    % --- 图 1: 聚合功率曲线对比 ---
    fig1 = figure('Name', 'F_Power_Curve', 'Color', 'w', 'Position', [100, 100, 600, 400]);
    hold on;
    x = 1:T_steps;
    
    % 绘制堆叠区域 (有约束)
    h_ac = area(x, results(2).P_AC, 'FaceColor', [0.6 0.8 1], 'EdgeColor', 'none', 'DisplayName', 'P_{AC} (协同)');
    h_tot = area(x, results(2).P_AC + results(2).P_EV, 'FaceColor', [0.6 1 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'P_{Total} (协同)');
    
    % 绘制线条 (无约束)
    plot(x, results(1).P_Total, 'r--', 'LineWidth', 1.5, 'DisplayName', 'P_{Total} (未协同)');
    
    % 修改：增加 FontWeight normal
    xlabel('时间步', 'FontName', font_name, 'FontSize', font_size, 'FontWeight', 'normal'); 
    ylabel('功率(MW)', 'FontName', font_name, 'FontSize', font_size, 'FontWeight', 'normal');
    legend([h_ac, h_tot], {'P_{AC} (协同)', 'P_{Total} (协同)'}, 'Location', 'northwest', 'FontName', font_name, 'FontSize', font_size);
    grid on; box on;
    set(gca, 'FontName', font_name, 'FontSize', font_size, 'FontWeight', 'normal');
    
    % 保存 图1
    print(fig1, '场景F峰值功率对比.png', '-dpng', '-r600');
    print(fig1, '场景F峰值功率对比.emf', '-dmeta');
    
    
    % --- 图 2: 物理指标对比 (峰值 & 波动率) ---
    fig2 = figure('Name', 'F_Physical_Metrics', 'Color', 'w', 'Position', [150, 150, 500, 400]);
    
    % 数据准备
    data_phy = [results(1).Max_Peak, results(2).Max_Peak; 
                results(1).Std_Dev*10, results(2).Std_Dev*10]; % 波动率x10以便同框
    
    b = bar(data_phy, 0.6);
    b(1).FaceColor = [0.8 0.3 0.3]; % 无约束 (红)
    b(2).FaceColor = [0.3 0.6 0.3]; % 有约束 (绿)
    
    % 修改：增加 FontWeight normal
    ylabel('功率值', 'FontName', font_name, 'FontSize', font_size, 'FontWeight', 'normal');
    set(gca, 'XTickLabel', {'峰值功率 (MW)', '波动率 (x10)'}, 'FontName', font_name, 'FontSize', font_size, 'FontWeight', 'normal');
    legend({'未协同', '协同'}, 'Location', 'best', 'FontName', font_name, 'FontSize', font_size);
    grid on; box on;
    
    % 数值标注
    for i = 1:2
        xt = b(i).XEndPoints;
        yt = b(i).YEndPoints;
        for j = 1:length(xt)
            % 修改：FontSize 10->12, 增加 FontWeight normal
            text(xt(j), yt(j), sprintf('%.4f', yt(j)), 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontName', font_name, 'FontWeight', 'normal');
        end
    end
    
    % 保存 图2
    print(fig2, '场景F物理指标对比.png', '-dpng', '-r600');
    print(fig2, '场景F物理指标对比.emf', '-dmeta');
    
    
    % --- 图 3: 风险成本 (CVaR) 对比 ---
    fig3 = figure('Name', 'F_Risk_Metrics', 'Color', 'w', 'Position', [200, 200, 400, 400]);
    
    data_risk = [results(1).CVaR, results(2).CVaR];
    
    b_risk = bar(data_risk, 0.5);
    b_risk.FaceColor = 'flat';
    b_risk.CData(1,:) = [0.8 0.3 0.3]; % 红
    b_risk.CData(2,:) = [0.3 0.6 0.3]; % 绿
    
    % 修改：增加 FontWeight normal
    ylabel('CVaR(元)', 'FontName', font_name, 'FontSize', font_size, 'FontWeight', 'normal');
    set(gca, 'XTickLabel', {'未协同', '协同'}, 'FontName', font_name, 'FontSize', font_size, 'FontWeight', 'normal');
    grid on; box on;
    
    % 数值标注
    xt = b_risk.XEndPoints;
    yt = b_risk.YEndPoints;
    % 修改：FontSize 12->14, 增加 FontWeight normal
    text(xt(1), yt(1), sprintf('%.4f', yt(1)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14, 'FontName', font_name, 'FontWeight', 'normal');
    text(xt(2), yt(2), sprintf('%.4f', yt(2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14, 'FontName', font_name, 'FontWeight', 'normal');
    
    % 增加下降箭头或文字说明
    mid_x = mean(xt);
    % 修改：移除 FontWeight bold，改为 normal
    text(mid_x, max(yt)*0.5, sprintf('风险减小\n-%.1f%%', imp_cvar), 'HorizontalAlignment', 'center', 'Color', 'b', 'FontWeight', 'normal', 'FontName', font_name, 'FontSize', 12);
    
    % 保存 图3
    print(fig3, '场景F风险对比.png', '-dpng', '-r600');
    print(fig3, '场景F风险对比.emf', '-dmeta');
end