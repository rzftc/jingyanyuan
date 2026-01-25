function P_base_sequence = EVbaseP_Hysteresis(C_EV, eta, E_tar, E_in, ...
                                              t_dep, t_in, dt, varargin)
% EVbaseP_Hysteresis 基于滞环控制(Bang-Bang)计算EV物理功率基线
% 
% 核心原理:
%   该函数不再假设车辆以恒定平均功率充电(物理不可达)，也不假设车辆即插即充(无调节能力)。
%   它构建一个围绕"理想平均能量轨迹"的约束通道(死区)，利用EV的 0/P_max 开关特性
%   在通道内进行滞环控制。
%
%   - 当电量落后于下限时 -> 全速充电 (P_max)
%   - 当电量超过上限时   -> 停止充电 (0)
%   - 当电量在通道内时   -> 保持上一时刻状态
%
%   这种策略在微观上是离散脉冲，宏观上紧贴调度目标，且保证离网时电量偏差在允许范围内。
%
% 输入参数 (与 EVbaseP_ChargeUntilFull 保持一致):
%   C_EV            - 电池容量 (kWh)
%   eta             - 充电效率 (0-1)
%   E_tar           - 目标电量 (kWh)
%   E_in            - 初始电量 (kWh)
%   t_dep           - 离网时间 (小时)
%   t_in            - 入网时间 (小时)
%   dt              - 时间步长 (小时)
%   varargin        - {1} r (偏差度), {2} p_on (额定功率), {3} SOC_init, {4} num_points, {5} time_axis
%
% 输出:
%   P_base_sequence - 功率基线序列 (kW)

    %% 1. 解析输入参数
    if length(varargin) < 5
        error('输入参数不足，请检查调用格式。');
    end
    
    r_flex = varargin{1};           % 用户设定的偏差度系数 (如 0.05 代表 5%)
    P_max = varargin{2};            % 额定充电功率 (kW)
    % SOC_initial = varargin{3};    % (本函数不直接使用，通过 E_in 体现)
    num_time_points = varargin{4};  % 总仿真点数
    time_points_abs = varargin{5};  % 绝对时间轴向量
    
    %% 2. 初始化
    P_base_sequence = zeros(num_time_points, 1);
    E_current = E_in;
    
    % 物理约束检查
    E_tar = min(E_tar, C_EV);
    E_current = min(E_current, C_EV);
    
    % 计算驻留时长
    T_dwell = t_dep - t_in;
    
    % 计算理想平均功率 P_avr (Schedule Reference)
    % 这是数学上的理想值，用于生成中心轨迹
    if T_dwell > 0 && E_tar > E_in
        P_avr = (E_tar - E_in) / (T_dwell * eta);
    else
        P_avr = 0;
    end
  
    %% 3. 定义滞环参数 (Deadband)
    % 设置允许波动的半带宽 delta
    % 建议取允许偏差的一半，以确保绝对安全，且防止最后时刻越界
    % 偏差度约束通常指: |E_final - E_tar| <= C * r
    safe_margin_factor = 0.1; 
    delta_width = C_EV * r_flex * safe_margin_factor;
    
    % 状态记忆变量 (0 或 P_max)
    last_power_state = 0; 
    
    %% 4. 时间步仿真循环
    for k = 1:num_time_points
        current_t = time_points_abs(k);
        
        % 判断是否在线
        is_online = (current_t >= t_in) && (current_t < t_dep);
        
        if is_online
            % --- A. 计算当前时刻的理想能量基准 (E_ideal) ---
            % 线性插值：从 t_in 到 t_dep，电量应从 E_in 均匀增加到 E_tar
            time_elapsed = current_t - t_in;
            
            % 理想累积充电量 (Ideal Expected Energy)
            E_ideal_now = E_in + P_avr * eta * time_elapsed;
            
            % 修正：理想轨迹不能超过物理目标(虽然P_avr已处理，但防万一)
            E_ideal_now = min(E_ideal_now, E_tar);
            
            % --- B. 构建滞环通道 (Corridor) ---
            E_upper_bound = E_ideal_now + delta_width;
            E_lower_bound = E_ideal_now - delta_width;
            
            % --- C. 滞环控制逻辑 (Bang-Bang Control) ---
            if E_current >= E_upper_bound
                % 触碰上限：充太快了，停机等待
                current_P = 0;
                
            elseif E_current <= E_lower_bound
                % 触碰下限：落后了，全速追赶
                current_P = P_max;
                
            else
                % 在通道中间：保持上一时刻状态 (惯性)
                % 这样可以避免在理想线附近高频抖动，形成较长的脉冲
                current_P = last_power_state;
            end
            
            % --- D. 物理截止约束 ---
            % 1. 电池容量限制
            if E_current >= C_EV
                current_P = 0;
            end
            
            % 2. 目标电量硬约束 (可选)
            % 如果接近离网时间且电量已达标，可以强制停止。
            % 但滞环逻辑本身会在 E_ideal 达到 E_tar 后，让 E_current 在 E_tar + delta 处停止。
            % 为了防止过充超过 E_tar 太多，这里加一个针对 E_tar 的软保护
            % (虽然 E_upper_bound 已经包含了 E_tar 的限制)
            if E_current >= E_tar + delta_width
                 current_P = 0;
            end

            % --- E. 更新状态 ---
            P_base_sequence(k) = current_P;
            
            % 能量积分
            E_charged = current_P * eta * dt;
            E_current = E_current + E_charged;
            
            % 记忆状态
            last_power_state = current_P;
            
        else
            % 离线状态
            P_base_sequence(k) = 0;
            last_power_state = 0;
            % 离线时不更新 E_current (假设无自放电)
        end
    end
end