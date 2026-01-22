function [DeltaP_pu_plus_t, DeltaP_pu_minus_t] = calculateEVAdjustmentPotentia(C_EV, r, eta, E_tar, E_in, E_current, t_dep, t_in, Pmax, Pmin, P_base, S_curr, delta_t_adj)
    % 输入参数说明:
    %   C_EV        - 电动汽车电池容量向量 [N×1] (kWh)
    %   r           - 电池充放电速率系数向量 [N×1]
    %   eta         - 充放电效率向量 [N×1] (0-1)
    %   E_tar       - 目标电量向量 [N×1] (kWh)
    %   E_in        - 初始电量向量 [N×1] (kWh)
    %   E_current   - 当前电量向量 [N×1] (kWh)
    %   t_dep       - 离网时间向量 [N×1] (分钟)
    %   t_in        - 入网时间向量 [N×1] (分钟)
    %   Pmax        - 最大允许功率向量 [N×1] (kW)
    %   Pmin        - 最小允许功率向量 [N×1] (kW)
    %   P_base      - 基线功率向量 [N×1] (kW)
    %   S_curr      - 当前虚拟SOC向量 [N×1] (0-1)
    %   delta_t_adj - 功率调节时间窗口 (分钟)
    
    % 参数维度一致性验证
    N = length(C_EV);
    
    % 初始化输出向量
    DeltaP_pu_plus_t = zeros(N, 1);
    DeltaP_pu_minus_t = zeros(N, 1);
    
    % 逐台电动汽车计算调节潜力
    for i = 1:N
        % ========== 功率约束潜力计算 ========== 
        % 式35-36: 计算功率上下限与基线的偏差
        DeltaP_P_plus = Pmax(i) - P_base(i);
        DeltaP_P_minus = Pmin(i) - P_base(i);
        
        % 确保功率潜力值在有效范围内
        DeltaP_P_plus = max(DeltaP_P_plus, 0);   % 正向潜力非负
        DeltaP_P_minus = min(DeltaP_P_minus, 0); % 负向潜力非正
        
        % ========== 需求功率计算 ========== 
        % 式8: 计算满足目标电量的需求功率
        if t_dep(i) == t_in(i)
            P_req = 0;  % 时间窗口为零时的保护处理
        else
            P_req = (E_tar(i) - E_in(i)) / (eta(i) * (t_dep(i) - t_in(i)));
        end
        
        % ========== 能量约束参数计算 ========== 
        % 式39: 定义能量约束参数
        A = -C_EV(i) * r(i) / eta(i);
        B = C_EV(i) * r(i) / eta(i);
        C = P_req;
        
        % ========== 能量约束潜力计算 ========== 
        % 式40-41: 计算考虑SOC的能量调节潜力
        term_plus = (-A/delta_t_adj) + (B*S_curr(i)/delta_t_adj) + C;
        DeltaP_E_plus = term_plus - P_base(i);
        
        term_minus = (A/delta_t_adj) + (B*S_curr(i)/delta_t_adj) + C;
        DeltaP_E_minus = term_minus - P_base(i);
        
        % ========== 潜力值边界处理 ========== 
        % 确保能量潜力值在有效范围内
        DeltaP_E_plus_adj = max(DeltaP_E_plus, 0);    % 正向潜力非负
        DeltaP_E_minus_adj = min(DeltaP_E_minus, 0);  % 负向潜力非正
        
        % ========== 综合约束处理 ========== 
        % 式44: 取功率约束和能量约束的交集
        DeltaP_pu_plus_t(i) = min(DeltaP_P_plus, DeltaP_E_plus_adj);
        DeltaP_pu_minus_t(i) = max(DeltaP_P_minus, DeltaP_E_minus_adj);
        
        % 当前电量已达目标时的特殊处理
        if E_current(i) >= E_tar(i)
            DeltaP_pu_plus_t(i) = 0;  % 禁止继续充电
            DeltaP_pu_minus_t(i) = 0;
        end
    end
end
