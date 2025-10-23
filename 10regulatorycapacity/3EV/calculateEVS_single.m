function [E_exp, E_current, P_current, SOC] = calculateEVS_single(m_3, E_exp, E_tar_ori, eta, E_current, P_current_in, C_EV, r, p_on, dt, t_dep, t_current)
    % calculateEVS_single: 根据控制信号和约束，计算单辆电动汽车（EV）在一个时间步长的状态。
    %
    % 输入:
    %   m_3          - 基础功率曲线或控制信号 (标量)
    %   E_exp        - 进入该时间步时的期望能量 (标量)
    %   E_tar_ori    - 原始设定的目标能量 (标量)
    %   eta          - 充电效率 (标量)
    %   E_current    - 进入该时间步时的当前能量 (标量)
    %   P_current_in - 进入该时间步时的功率 (此函数内会重新决策，可视为前一时刻状态)
    %   C_EV         - 电池容量 (标量)
    %   r            - SOC归一化因子 (标量)
    %   p_on         - EV最大充电功率 (标量)
    %   dt           - 时间步长 (小时)
    %   t_dep        - 计划离开时间 (绝对时间)
    %   t_current    - 当前时间 (绝对时间)
    %
    % 输出:
    %   E_exp        - 更新后的期望能量
    %   E_current    - 更新后的当前能量
    %   P_current    - 此时间步最终决策的充电功率
    %   SOC          - 计算出的虚拟电池状态（Virtual State of Charge）

    % --- 1. 初始化和计算剩余时间 ---
    t_rem = max(0, t_dep - t_current); % 剩余充电时间，避免负数
    P_current = P_current_in; % 默认使用上一时刻的功率

    % --- 2. 计算虚拟SOC ---
    % 这个SOC是基于期望能量和当前能量差异的控制变量，并非物理SOC（E_current/C_EV）
    SOC = -(E_current - E_exp) / (C_EV * r);

    % --- 3. 决策当前时间步的充电功率 P_current ---
    % 决策优先级: 紧急情况 > 常规调节
    
    % **紧急情况判断**
    % a) 如果当前电量已满足或超过目标电量，则必须停止充电。
    if E_current >= E_tar_ori
        P_current = 0;
    % b) 如果剩余时间内不以最大功率充电就无法达到目标，则必须以最大功率充电。
    elseif (E_tar_ori - E_current) >= t_rem * p_on * eta
        P_current = p_on;
    % **常规调节**
    % 在非紧急情况下，根据虚拟SOC进行调节
    else
        if SOC >= 1
            P_current = p_on; % 能量低于期望，增加功率
        elseif SOC <= -1
            P_current = 0;    % 能量高于期望，减小功率
        % else 在-1和1之间，保持上一时刻的功率 P_current_in
        end
    end
    
    % --- 4. 更新能量状态 ---
    % 根据最终确定的 P_current 来更新当前能量和期望能量
    E_exp = E_exp + eta * m_3 * dt;
    E_current = E_current + eta * P_current * dt;
    
    % --- 5. 确保能量不会超过电池物理容量 ---
    if E_current > C_EV
        E_current = C_EV;
    end
end