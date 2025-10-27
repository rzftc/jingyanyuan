function Delta_E = calculateDeltaE(EV, P_real)
    P_r = P_real;
    
    if P_r >= EV.P_h_max
        % 条件1：P_r ≥ P_h_max
        Delta_E = EV.Delta_E_h_max;
        
    elseif P_r >= EV.P_0
        % 条件2：P_h_max ≥ P_r ≥ P_0 
        delta_p = P_r - EV.P_0;
        term1 = -(EV.Delta_E_h_max * delta_p^2) / (EV.P_h_max - EV.P_0)^2;
        term2 = (2 * EV.Delta_E_h_max * delta_p) / (EV.P_h_max - EV.P_0);
        Delta_E = term1 + term2;
        
    elseif P_r >= EV.P_l_min
        % 条件3：P_0 ≥ P_r ≥ P_l_min（关键修正）
        delta_p = P_r - EV.P_0;  % 修正Δp定义
        numerator = delta_p^2 * (3 * (EV.P_0 - EV.P_l_min) + 2 * delta_p);
        denominator = (EV.P_0 - EV.P_l_min)^3;
        Delta_E = EV.Delta_E_q_max * numerator / denominator;  % 变量名修正
        
    else
        % 条件4：P_r ≤ P_l_min
        Delta_E = EV.Delta_E_q_max;  % 变量名修正
    end
end
