function EV = calculateVirtualSOC(EV, t, dt)
    % 判断是否在并网时段
    if t < EV.t_in || t > EV.t_dep
        EV.S_original = 0; % 原始SOC
        EV.S_modified = 0; % 修正SOC
        return;
    end
    
    %% ----------------- 原始SOC计算（公式24）-----------------
    % 计算期望电量（公式22-23）
    P_req = (EV.E_tar - EV.E_ini) / (EV.eta * (EV.t_dep - EV.t_in) / 60);
    if P_req > EV.P_N
        warning('EV%d: P_req=%.2f > P_N=%.2f', EV.EV_ID, P_req, EV.P_N);
    end
    EV.E_exp = EV.E_exp + EV.eta * P_req * (dt / 60);
    
    % 实际电量变化（假设理想充电）
    EV.E_actual = EV.E_actual + EV.eta * EV.P_current * (dt / 60);
    
    % 原始SOC（公式24）
    EV.S_original = -(EV.E_actual - EV.E_exp) / (EV.C * EV.r);
    % if EV.S_original>1
    %     EV.S_original=1;
    % elseif EV.S_original<-1
    %      EV.S_original=-1;
    % end

    %% ---------- 修正SOC计算（基于 I 指标） ----------
    % 1) 额定功率一次性充满所需时间 t_ch  (单位: min)
    t_ch = 60 * (EV.E_tar - EV.E_ini) / EV.P_N;      % min
    
    % 2) 初始化剩余充电时间 tau_rem（只在首次并网调用时运行）
    if ~isfield(EV, 'tau_rem')
        EV.tau_rem = t_ch;                           % τ(t_in) = t_ch
    end
    
    % 3) 递推 tau_rem   τ_{k+1} = τ_k - (P/P_N)·Δt
    %    只在并网窗口内更新，且不允许为负
    if t >= EV.t_in && t <= EV.t_dep
        EV.tau_rem = max(EV.tau_rem - (EV.P_current/EV.P_N) * dt, 0);
    end
    
    % 4) 计算剩余时间比率 ρ(t) = τ / (t_dep - t)
    time_left = max(EV.t_dep - t, eps);              % eps 避免除 0
    rho       = EV.tau_rem / time_left;
    
    % 5) 指标 I(t) = tanh[ k (ρ - 0.5) ]    —— 取值范围 (-1,1)
    kappa   = 4;                                     % 陡峭度，可根据需要调整
    I_value = tanh( kappa * (rho - 0.5) );
    
    % 6) 综合得到修正 SOC
    alpha1 = 1;            % 对原始 SOC 的权重
    alpha2 = 0;            % 对指标   I 的权重
    EV.S_modified = alpha1 * EV.S_original + alpha2 * I_value;
end
