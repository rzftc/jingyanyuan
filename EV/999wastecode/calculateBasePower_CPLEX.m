% function [P_base, S_agg] = calculateBasePower_CPLEX(EVs, current_time, H, dt, M1, M2, M3)
%     % 参数验证
%     validateattributes(M1, {'double'}, {'scalar', 'real'});
%     validateattributes(M2, {'double'}, {'scalar', 'real'});
%     validateattributes(M3, {'double'}, {'scalar', 'real'});
%     validateattributes(dt, {'double'}, {'positive', 'scalar'});
% 
%     % 获取功率边界
%     [P_min, P_max] = predictPowerLimits(EVs, current_time, H);
%     P_min = max(P_min, 0);
% 
%     %% 构建CPLEX模型
%     cplex = Cplex();
%     cplex.Model.sense = 'minimize';
% 
%     % ========== 关键修正1：变量定义 ==========
%     num_vars = H;
%     obj = ones(num_vars, 1);   % 必须为列向量 (原错误原因)
%     lb  = -ones(num_vars, 1);   % 列向量
%     ub  = ones(num_vars, 1);    % 列向量
%     cplex.addCols(obj, [], lb, ub);  % 正确参数格式
% 
%     %% 动态方程约束（论文式34）
%     A_eq = [];
%     b_eq = [];
%     S_agg_current = mean([EVs.S_original]);
% 
%     for t = 1:H
%         if t == 1
%             % 方程: M1*S1 + M2*S0 = P_min - M3
%             row = sparse(1, 1, M1, 1, num_vars);
%             rhs = P_min(t) - (M2*S_agg_current + M3);
%         else
%             % 方程: M1*St + M2*S(t-1) = P_min - M3
%             row = sparse(1, [t-1, t], [M2, M1], 1, num_vars);
%             rhs = P_min(t) - M3;
%         end
%         A_eq = [A_eq; row];
%         b_eq = [b_eq; rhs];
%     end
% 
%     %% 功率边界约束（论文式31）
%     A_ineq = [];
%     b_low = [];
%     b_high = [];
% 
%     for t = 1:H
%         if t == 1
%             row = sparse(1, 1, M1, 1, num_vars);
%             const_term = M2*S_agg_current + M3;
%         else
%             row = sparse(1, [t-1, t], [M2, M1], 1, num_vars);
%             const_term = M3;
%         end
%         A_ineq = [A_ineq; row];
%         b_low  = [b_low;  P_min(t) - const_term];
%         b_high = [b_high; P_max(t) - const_term];
%     end
% 
%     % ========== 关键修正2：正确添加约束 ==========
%     cplex.addRows(b_low, A_ineq, b_high);  % 不等式约束
%     cplex.addRows(b_eq,  A_eq,   b_eq);    % 等式约束
% 
%     %% 求解参数
%     cplex.Param.simplex.tolerances.feasibility.Cur = 1e-6;
%     cplex.solve();
% 
%     %% 结果处理
%     if cplex.Solution.status == 1
%         S_agg = cplex.Solution.x;
%         P_base = M1*[S_agg(2:end); 0] + M2*[S_agg_current; S_agg(1:end-1)] + M3;
%         P_base = P_base(1:H); % 截断最后一个无效值
%     else
%         warning('安全模式：使用边界中值');
%         P_base = (P_min + P_max)/2;
%         S_agg = linspace(S_agg_current, 0, H)';
%     end
% 
%     % ========== 关键修正3：SOC限幅 ==========
%     S_agg = max(min(S_agg, 1), -1);
% end
