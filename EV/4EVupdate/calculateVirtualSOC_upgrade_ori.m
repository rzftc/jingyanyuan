function EV = calculateVirtualSOC_upgrade_ori(EV, t, dt)
%CALCULATEVIRTUALSOC_UPGRADE 计算电动汽车的原始及修正虚拟SOC。
%
%   输入参数:
%       EV: 结构体，包含电动汽车状态和参数:
%           t_in, t_dep     - 入网/离网时间 (min)
%           E_ini, E_tar    - 初始/目标电量 (kWh)
%           P_N             - 额定充电功率 (kW)
%           eta             - 充放电效率
%           C               - 电池容量 (kWh)
%           r               - S_original 计算因子
%           P_current       - 当前充电功率 (kW)
%           EV_ID           - 车辆ID
%           E_exp, E_actual - 期望/实际累计充电量 (内部更新)
%           tau_rem         - 剩余充电时间 (内部更新)
%       t:  当前时间 (min)
%       dt: 时间步长 (min)
%
%   输出参数:
%       EV: 更新后的电动汽车结构体

    % 非并网时段，返回零值SOC
    if t < EV.t_in || t > EV.t_dep
        EV.S_original = 0;
        EV.S_modified = 0;
        return;
    end

    %% ----------------- 原始SOC计算 (S_original) -----------------
    % 计算反映实际充电量与期望(匀速)充电量偏差的 SOC

    % 期望平均功率
    P_req = (EV.E_tar - EV.E_ini) / (EV.eta * (EV.t_dep - EV.t_in) / 60);
    if P_req > EV.P_N
        warning('EV%d: 期望平均功率 P_req=%.2f > 额定功率 P_N=%.2f', EV.EV_ID, P_req, EV.P_N);
    end

    % 更新期望与实际累计充电量
    EV.E_exp = EV.E_exp + EV.eta * P_req * (dt / 60);       % 单位: kWh
    EV.E_actual = EV.E_actual + EV.eta * EV.P_current * (dt / 60); % 单位: kWh

    % 计算原始虚拟SOC
    if EV.C > 0 && EV.r ~= 0 % 避免除零
        EV.S_original = -(EV.E_actual - EV.E_exp) / (EV.C * EV.r);
    else
        EV.S_original = 0;
    end
    % % 可选: 限制 S_original 在 [-1, 1]
    % EV.S_original = max(min(EV.S_original, 1), -1);


    %% ---------- 修正SOC计算 (S_modified, 基于改进I指标) ----------
    % 计算反映充电紧急程度的 SOC

    % 1) 额定功率下充满所需总时间 t_ch (假设从 E_ini 到 E_tar)
    energy_to_charge = max(EV.E_tar - EV.E_ini, 0);
    if EV.P_N <= 0
        t_ch = inf;
    else
        % 注意: 此处 t_ch 未考虑充电效率 eta (与原代码一致)
        t_ch = 60 * energy_to_charge / EV.P_N; % 单位: min
    end

    % 2) 初始化剩余充电时间 tau_rem (仅首次调用时)
    if ~isfield(EV, 'tau_rem')
        EV.tau_rem = t_ch; % τ(t_in) = t_ch
    end

    % 3) 更新剩余充电时间 tau_rem: τ_{k+1} = τ_k - (P_current / P_N) * Δt
    if t >= EV.t_in && t <= EV.t_dep % 仅在并网时段内
        if EV.P_N > 0 % 避免除零
            tau_reduction = (EV.P_current / EV.P_N) * dt;
            EV.tau_rem = max(EV.tau_rem - tau_reduction, 0);
        end
    end

    % 4) 计算剩余时间比率 rho: ρ(t) = τ_rem(t) / (t_dep - t)
    time_left = max(EV.t_dep - t, eps); % 剩余可用时间 (min), eps避免除零
    rho = EV.tau_rem / time_left;       % 反映充电时间压力

    % 5) 计算改进后的紧急程度指标 I(t)
    kappa = 4;    % tanh 函数陡峭度
    gamma = 0.05; % 时间衰减项强度 (单位: 1/min), 需调优

    t_elapsed = max(t - EV.t_in, 0); % 入网后经过时间 (min)

    % 指标 I 计算公式
    argument = kappa * (rho - 0.5) - gamma * (EV.P_current / EV.P_N) * t_elapsed;
    I_value = tanh(argument); % I_value 范围为 (-1, 1)

    % 6) 计算修正虚拟 SOC (S_modified)
    %    !! 注意：当前 alpha 设置使 S_modified = S_original (未启用指标I)。
    %    !! 若需使用指标I, 请调整 alpha1/alpha2 (如设为0/1)。
    alpha1 = 0.8;   % S_original 权重
    alpha2 = 0.2;   % 指标 I_value 权重
    EV.S_modified = alpha1 * EV.S_original + alpha2 * I_value;

    % % 可选: 限制 S_modified 在 [-1, 1]
    % EV.S_modified = max(min(EV.S_modified, 1), -1);

end