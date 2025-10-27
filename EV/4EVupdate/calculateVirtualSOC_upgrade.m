%% calculateVirtualSOC_upgrade.m (优化尝试版)
function EV = calculateVirtualSOC_upgrade(EV, t, dt)
    % 非并网时段，返回零值SOC
    if t < EV.t_in || t > EV.t_dep %
        EV.S_original = 0;
        EV.S_modified = 0;
        return; %
    end

    %% ----------------- 原始SOC计算 (S_original) -----------------
    % 期望平均功率 (避免重复计算时间差)
    time_diff_hours = (EV.t_dep - EV.t_in) / 60; % 单位: 小时
    if time_diff_hours <= 0 || EV.eta <= 0 % 防御性编程，避免除零或无效值
        P_req = 0;
         if time_diff_hours <= 0 && EV.E_tar > EV.E_ini
             % 如果时间窗口为0但仍需充电，可能是一个数据问题
             if isfield(EV, 'EV_ID')
                warning('calculateVirtualSOC_upgrade:ZeroTime', 'EV %d: 充电时间窗口为零或负，但仍需充电 (t_in=%.1f, t_dep=%.1f)', EV.EV_ID, EV.t_in, EV.t_dep);
             else
                warning('calculateVirtualSOC_upgrade:ZeroTime', '充电时间窗口为零或负，但仍需充电');
             end
         end
    else
        P_req = (EV.E_tar - EV.E_ini) / (EV.eta * time_diff_hours); %
    end
    % 更新期望与实际累计充电量
    dt_hours = dt / 60; % 预计算步长时间(小时)
    EV.E_exp = EV.E_exp + EV.eta * P_req * dt_hours; %
    EV.E_actual = EV.E_actual + EV.eta * EV.P_current * dt_hours; %

    % 计算原始虚拟SOC
    denominator_s_orig = EV.C * EV.r; % 预计算分母
    if denominator_s_orig ~= 0 % 避免除零
        EV.S_original = -(EV.E_actual - EV.E_exp) / denominator_s_orig; %
    else
        EV.S_original = 0; %
    end

    %% ---------- 修正SOC计算 (S_modified, 基于改进I指标) ----------
    % 1) 额定功率下充满所需总时间 t_ch (优化判断)
    energy_to_charge = EV.E_tar - EV.E_ini; % 直接用减法，下面判断非负
    if EV.P_N <= 0 || energy_to_charge <= 0 % 如果额定功率<=0 或 无需充电
        t_ch = 0; %
    else
        t_ch = 60 * energy_to_charge / EV.P_N; %
    end

    % 2) 初始化剩余充电时间 tau_rem (仅首次调用时)
    %   优化: 直接检查字段是否存在即可，无需 isempty
    if ~isfield(EV, 'tau_rem') %|| isempty(EV.tau_rem) % isempty 检查通常不是必需的
        EV.tau_rem = t_ch; %
    end

    % 3) 更新剩余充电时间 tau_rem (优化判断)
    % if t >= EV.t_in && t <= EV.t_dep % 此判断在函数开头已做，此处冗余
    if EV.P_N > 0 % 只需判断 P_N
        tau_reduction = (EV.P_current / EV.P_N) * dt; %
        EV.tau_rem = max(EV.tau_rem - tau_reduction, 0); %
    end
    % end

    % 4) 计算剩余时间比率 rho (优化变量名)
    time_remaining = EV.t_dep - t; %
    if time_remaining <= eps % 使用 eps 避免除零和负值
        rho = inf; % 如果剩余时间为0或负，认为时间压力无穷大
    else
        rho = EV.tau_rem / time_remaining; %
    end

    % 5) 计算改进后的紧急程度指标 I(t) (优化常数和计算)
    persistent kappa gamma % 使用 persistent 避免重复赋值
    if isempty(kappa)
        kappa = 4; %
        gamma = 0.05; %
    end

    t_elapsed = max(t - EV.t_in, 0); %

    if EV.P_N > 0
        power_ratio = EV.P_current / EV.P_N; %
    else
        power_ratio = 0; % 避免除零
    end
    argument = kappa * (rho - 0.5) - gamma * power_ratio * t_elapsed; %
    I_value = tanh(argument); %

    % 6) 计算修正虚拟 SOC (S_modified)
    persistent alpha1 alpha2 % 使用 persistent
    if isempty(alpha1)
        alpha1 = 0.8; %
        alpha2 = 0.2; %
    end
    EV.S_modified = alpha1 * EV.S_original + alpha2 * I_value; %
end