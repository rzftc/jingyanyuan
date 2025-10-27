function [P_base_sequence, S_sequence] = EVbaseP(EV, H, dt)
    % generatePbase_singleEV - 根据“尽快充满”策略为单个EV计算基线功率和SOC。
    % 该函数模拟EV从接入电网开始就以其额定功率充电，直到充满或离开。
    % 仿真时间范围固定为从第一天早上6点到第二天早上6点。
    %
    % 输入:
    %   EV            - 单个EV对象
    %   H             - 总仿真时长 (小时, 例如, 24)
    %   dt            - 时间步长 (小时, 例如, 0.5)
    %
    % 输出:
    %   P_base_sequence - 基线功率的轨迹 (1 x (H/dt))
    %   S_sequence      - 虚拟SOC的轨迹 (1 x (H/dt))

    % --- 参数提取 ---
    E_current = EV.E_ini;
    E_tar = EV.E_tar;
    t_in_minutes = EV.t_in;
    t_dep_minutes = EV.t_dep;
    p_on = EV.P_N;
    eta = EV.eta;
    
    % --- 初始化 ---
    num_steps = round(H / dt); % 计算总的仿真步数
    P_base_sequence = zeros(1, num_steps);
    S_sequence = zeros(1, num_steps);
    t_start_minutes = 6 * 60; % 仿真从早上6点开始

    % --- 遍历所有时间步 ---
    for t_idx = 1:num_steps
        % 计算当前时间点相对于第一天0点的绝对分钟数
        current_minute_abs = t_start_minutes + (t_idx - 1) * dt;

        % 检查车辆是否在线 (包含跨天逻辑)
        is_online = false;
        if t_dep_minutes > t_in_minutes
            % Case 1: 停车不跨天
            is_online = (current_minute_abs >= t_in_minutes) && (current_minute_abs < t_dep_minutes);
        else
            % Case 2: 停车跨天 (例如，晚上8点到第二天早上7点)
            % 在线时段为 [t_in, 1440) U [0, t_dep)
            is_online = (current_minute_abs >= t_in_minutes) || (current_minute_abs < t_dep_minutes);
        end

        % 计算功率
        if is_online && E_current < E_tar && p_on > 0
            % 在线、未满、能充电 -> 以额定功率充电
            P_base_sequence(t_idx) = p_on;
            E_current = E_current + p_on * eta * dt/60;
            if E_current > E_tar
                E_current = E_tar; % 防止过度充电
            end
        else
            % 离线或已充满 -> 功率为0
            P_base_sequence(t_idx) = 0;
        end
        
        % % 计算并存储当前步的虚拟SOC (S_original)
        % if (EV.E_tar - EV.E_exp) > 1e-6 % 避免除以零
        %     S_sequence(t_idx) = (E_current - EV.E_exp) / (EV.E_tar - EV.E_exp);
        % else
        %     S_sequence(t_idx) = (E_current >= EV.E_tar); % 如果E_exp和E_tar相等，充满则为1
        % end
    end
end
