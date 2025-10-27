function EV = updateLockState(EV, t)
    % UPDATELOCKSTATE 根据电量需求及SOC边界更新EV闭锁状态和子状态
    % 输入参数:
    %   EV  - 电动汽车对象 (结构体)
    %   t   - 当前时间/分钟 (标量)
    % 输出参数:
    %   EV  - 更新状态后的电动汽车对象 (结构体)

    %% --------------------- 输入校验 ---------------------
    % 离网时间未定义或当前时间超限时强制转为LockOFF
    if isempty(EV.t_dep) || t > EV.t_dep
        EV.state    = 'LockOFF';
        EV.substate = 'OFF';
        EV.P_current = 0;
        return;
    end
    
    %% ----------------- 状态基础判断 -----------------
    % 根据时间窗口设置基础状态
    if t >= EV.t_in && t <= EV.t_dep
        EV.state = 'ON';
    elseif t < EV.t_in
        EV.state = 'LockOFF';
    elseif t>EV.t_dep
        EV.state='LockOFF';
    end 
    
    %% --------------------- 核心计算 ---------------------
    % 计算剩余可用时间（分钟）和待充电量（kWh）
    t_remain        = EV.t_dep - t; 
    required_energy = EV.E_tar - EV.E_actual; 

    %% ----------------- 闭锁状态决策 -----------------
    % 条件1: 电量已满足目标（考虑浮点计算容差）
    if EV.E_actual >= EV.E_tar - 1e-6
        EV.state    = 'LockOFF';
        EV.substate = 'OFF';

    % 条件2: 剩余时间不足完成最低充电需求
    elseif required_energy >= EV.eta * EV.P_N * (t_remain / 60)
        EV.state    = 'LockON';
        EV.substate = 'ON';

    % 条件3: 正常状态动态切换
    else
        % 根据充放电功率方向切换子状态
        if EV.P_current > 0
            EV.substate = 'ON';  % 充电状态
        else
            EV.substate = 'OFF'; % 空闲状态
        end
    end

    %% ----------------- SOC边界保护 -----------------
    % 当修正SOC触及边界时强制状态切换（防止数值振荡）
    if EV.S_modified >= 1 - 1e-6
        EV.substate = 'ON';  % SOC上限保护
    elseif EV.S_modified <= -1 + 1e-6
        EV.substate = 'OFF'; % SOC下限保护
    end
    % if EV.S_original >= 1 - 1e-6
    %     EV.substate = 'ON';  % SOC上限保护
    % elseif EV.S_original <= -1 + 1e-6
    %     EV.substate = 'OFF'; % SOC下限保护
    % end
end
