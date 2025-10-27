function [DeltaP_plus, DeltaP_minus] = calculateEVAdjustmentPotentia_new(E_reg_min, E_reg_max, E_current, t_dep, t_current, p_on, P_base, eta, dt)
% calculateEVAdjustmentPotentia: 根据能量灵活性窗口计算调节潜力
% (此版本已根据 calculateEVAdjustmentPotentia_new.m 的逻辑完全改造)
%
% 输入:
%   E_reg_min   - 可接受的能量下限 (kWh)
%   E_reg_max   - 可接受的能量上限 (kWh)
%   E_current   - 当前电量 (kWh)
%   t_dep       - 离网时间 (小时, 绝对时间)
%   t_current   - 当前时间 (小时, 绝对时间)
%   p_on        - 额定充电功率 (kW)
%   P_base      - 当前时间步的基准功率 (kW)
%   eta         - 充电效率
%   dt          - 时间步长 (小时)
%
% 输出:
%   DeltaP_plus  - 上调潜力 (kW, >= 0)
%   DeltaP_minus - 下调潜力 (kW, <= 0)

    % 计算剩余并网时间
    t_rem = t_dep - t_current;
    
    % 如果剩余时间不足一个步长，则无调节能力
    if t_rem < dt
        DeltaP_plus = 0;
        DeltaP_minus = 0;
        return;
    end

    % --- 上调潜力计算 (Up-regulation) ---
    % 逻辑：当前时间步能增加的最大功率，使得最终电量不超过 E_reg_max。
    %      我们做一个最保守的假设：此步之后，未来都以最低功率（0）充电。
    
    % 在当前步充电后，电量不能超过 E_reg_max
    % E_current + P_up_potential * eta * dt <= E_reg_max
    % P_up_potential <= (E_reg_max - E_current) / (eta * dt)
    max_power_from_energy_limit = (E_reg_max - E_current) / (eta * dt);
    
    % 实际可用的最大充电功率受物理额定功率和能量上限双重约束
    potential_up_power = min(p_on, max_power_from_energy_limit);
    
    % 上调潜力是相对于基准功率的增量，且必须为正
    DeltaP_plus = max(0, potential_up_power - P_base);

    % --- 下调潜力计算 (Down-regulation) ---
    % 逻辑：当前时间步能减少的最小功率，使得最终电量不低于 E_reg_min。
    %      我们做一个最乐观的假设：此步之后，未来都以最大功率（p_on）充电。

    % 未来所有剩余时间（不含当前步）能充入的最大电量
    max_energy_from_future_steps = p_on * eta * (t_rem - dt);
    
    % 为了最终达到 E_reg_min，当前时间步至少需要充入的电量
    energy_needed_now = E_reg_min - E_current - max_energy_from_future_steps;
    
    % 换算成当前时间步至少需要的功率
    min_power_from_energy_limit = energy_needed_now / (eta * dt);
    
    % 实际可用的最小充电功率受物理下限（0）和能量下限双重约束
    potential_down_power = max(0, min_power_from_energy_limit);
    
    % 下调潜力是相对于基准功率的减量，且必须为负
    DeltaP_minus = min(0, potential_down_power - P_base);

end