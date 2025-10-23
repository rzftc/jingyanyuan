function P_base_sequence = EVbaseP_ChargeUntilFull(C_EV, eta, E_tar, E_in, ...
                                                 t_dep, t_in, dt, varargin)
% EVbaseP_ChargeUntilFull 计算电动汽车 (EV) 的基线功率序列。
% EV 从接入电网 (t_in) 开始，以其额定功率 (p_on) 充电，
% 直到电量达到目标电量 (E_tar) 或到达离开时间 (t_dep)。
% 在其他时间或充满电后，功率为零。

% 输入参数:
%   C_EV                - EV 电池容量 (kWh)
%   eta                 - 充电效率 (0 到 1)
%   E_tar               - 目标电量水平 (kWh)。当 E_current >= E_tar 时，EV停止充电。
%   E_in                - 接入电网时的初始电量水平 (kWh)
%   t_dep               - 离开时间 (小时, 从仿真开始计算，例如 17.5 代表下午5:30)
%   t_in                - 接入时间 (小时, 从仿真开始计算，例如 8.0 代表上午8:00)
%   dt                  - 仿真时间步长 (小时)
%   varargin            - 可变输入参数，以匹配原始调用签名。
%                         期望输入: r, p_on, SOC_initial, num_time_points
%                           - r (此函数未使用)
%                           - p_on: EV的额定充电功率 (kW)
%                           - SOC_initial (此函数未使用 E_in)
%                           - num_time_points: 仿真时间范围内的总时间点数
% 输出参数:
%   P_base_sequence     - 基线功率序列 (列向量, num_time_points x 1)，单位 kW。
    
    if length(varargin) < 4
        error('输入参数不足。期望通过varargin输入p_on和num_time_points。');
    end
    
    p_on = varargin{2}; % 额定充电功率
    num_time_points = varargin{4}; % 总时间点数
    
    % 初始化输出序列和当前电量状态
    P_base_sequence = zeros(num_time_points, 1);
    E_current = E_in;

    % 确保电量水平的实际限制
    if E_tar > C_EV
        E_tar = C_EV; % 目标电量不能超过电池容量
    end
    if E_current > C_EV
        E_current = C_EV; % 初始电量不能超过电池容量
    end
 
    % 遍历仿真时间范围内的每个时间点
    for t_idx = 1:num_time_points
        current_time = (t_idx - 1) * dt; % 当前时间区间的开始时间

        % 判断 EV 当前是否接入电网 (在线)
        % EV 在线条件: current_time 在 [t_in, t_dep) 区间内
        % 这意味着如果 t_dep 晚于当前时间，则可以在从 current_time 开始的区间内充电。
        is_online = (current_time >= t_in) && (current_time < t_dep);

        if is_online
            if E_current < E_tar && p_on > 0 
                % EV 需要充电且额定充电功率大于0
                
                P_base_sequence(t_idx) = p_on; % 以额定功率充电
                
                % 更新当前电量水平
                energy_charged_this_step = p_on * eta * dt;
                E_current = E_current + energy_charged_this_step;
                
                % 如果达到或超过 E_tar，将 E_current 限制为 E_tar。
                % EV 在下一个时间步将被视为充满。
                if E_current >= E_tar
                    E_current = E_tar;
                end
                
                % 确保 E_current 不会物理上超过 C_EV
                if E_current > C_EV
                    E_current = C_EV;
                end
                
            else
                % EV 在线但已充满 (E_current >= E_tar) 或 p_on 为 0。
                P_base_sequence(t_idx) = 0;
            end
        else
            % EV 不在线 (在 t_in 之前或在 t_dep 及之后)
            P_base_sequence(t_idx) = 0;
            % 注意: E_current 在离线时保持不变 (此处未模拟自放电)。
            % 如果 EV 稍后重新接入 (在此模型结构中不典型，但需注意)，
            % 它将以其最后的 E_current 继续。
        end
    end
end