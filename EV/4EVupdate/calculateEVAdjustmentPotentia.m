function [DeltaP_pu_plus_t, DeltaP_pu_minus_t] = calculateEVAdjustmentPotentia(C_EV, r, eta, E_tar, E_in, E_current, t_dep, t_in, Pmax, Pmin, P_base, S_curr, delta_t_adj)
    % 输入参数说明:
    % C_EV       - 电池容量向量 [Nx1]
    % eta        - 充放电效率向量 [Nx1]
    % E_tar      - 目标电量向量 [Nx1]
    % E_in       - 初始电量向量 [Nx1]
    % t_dep      - 离网时间向量 [Nx1]
    % t_in       - 入网时间向量 [Nx1]
    % Pmax       - 最大功率限制向量 [Nx1]
    % Pmin       - 最小功率限制向量 [Nx1]
    % P_base     - 基线功率向量 [Nx1]
    % S_curr     - 当前虚拟SOC向量 [Nx1]
    % delta_t_adj- 调节时长标量
    
    % 参数验证
    N = length(C_EV);
   
    DeltaP_pu_plus_t = zeros(N,1);
    DeltaP_pu_minus_t = zeros(N,1);
    
    for i = 1:N
        % 1. 计算功率约束潜力 (式35-36)
        DeltaP_P_plus = Pmax(i) - P_base(i);
        DeltaP_P_minus = Pmin(i) - P_base(i);
        DeltaP_P_plus=max(DeltaP_P_plus,0);
        DeltaP_P_minus=min(DeltaP_P_minus,0);
        % 2. 计算需求功率P_req,i (式8)
        if t_dep(i) == t_in(i)
            P_req = 0; % 避免除零
        else
            P_req = (E_tar(i) - E_in(i)) / (eta(i) * (t_dep(i) - t_in(i)));
        end
        
        % 3. 计算能量约束参数 (式39)
        A = -C_EV(i)*r(i)/eta(i);
        B = C_EV(i)*r(i)/eta(i);
        C = P_req;
        
        % 4. 计算能量约束潜力 (式40-41)
        term_plus = (-A/delta_t_adj) + (B*S_curr(i)/delta_t_adj) + C;
        DeltaP_E_plus = term_plus - P_base(i);
        
        term_minus = (A/delta_t_adj) + (B*S_curr(i)/delta_t_adj) + C;
        DeltaP_E_minus = term_minus - P_base(i);
        
        % 5. 边界修正 (式42-43)
        DeltaP_E_plus_adj = max(DeltaP_E_plus, 0);
        DeltaP_E_minus_adj = min(DeltaP_E_minus, 0);
        
        % 6. 综合约束 (式44)
        DeltaP_pu_plus_t(i) = min(DeltaP_P_plus, DeltaP_E_plus_adj);
        DeltaP_pu_minus_t(i) = max(DeltaP_P_minus, DeltaP_E_minus_adj);
        if E_current(i) >= E_tar(i)
            DeltaP_pu_plus_t(i)=0;
        end
    end
end



% % calculateEVAdjustmentPotentia 中增加容错
% function [DeltaP_pu_plus_t, DeltaP_pu_minus_t] = calculateEVAdjustmentPotentia(C_EV, eta,Pbase, Pmax, Pmin,m_3, S_curr, delta_t_adj, t)
%     if Pbase(t) == 0
%         DeltaP_pu_plus_t = 0;
%         DeltaP_pu_minus_t = 0;
%     else
%         % DeltaP_pu_plus = (Pmax - Pbase(t)) / Pbase(t);
%         % DeltaP_pu_minus = (Pbase(t) - Pmin) / Pbase(t);
%         DeltaP_pu_plus = (Pmax - Pbase(t)) ;
%         DeltaP_pu_minus = (Pmin-Pbase(t)) ;
%         % 修正公式：基于论文式(38)-(39)
%         A=-C_EV/eta;
%         B=C_EV/eta;
%         C=m_3;
%         term_1 = -A/delta_t_adj + B*S_curr/delta_t_adj+C;
%         term_2 = A/delta_t_adj + B*S_curr/delta_t_adj+C;
%         DeltaP_E_pu_plus = term_1 - Pbase(t);
%         DeltaP_E_pu_minus = term_2-Pbase(t);
%         % DeltaP_E_pu_plus = (term_1 - Pbase(t)) / Pbase(t);
%         % DeltaP_E_pu_minus = (Pbase(t) - term_2) / Pbase(t);
%         % 限制结果范围
%         if DeltaP_E_pu_plus<0
%             DeltaP_E_pu_plus=0;
%         end
%         if DeltaP_E_pu_minus>0
%             DeltaP_E_pu_minus=0;
%         end
%         DeltaP_pu_plus_t = min(DeltaP_pu_plus, DeltaP_E_pu_plus);  % 防止Inf
%         DeltaP_pu_minus_t = max(DeltaP_pu_minus, DeltaP_E_pu_minus);
%     end
% end
