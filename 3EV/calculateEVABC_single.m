function [m_1, m_2, m_3] = calculateEVABC_single(C_EV, eta, E_tar, E_in, t_dep, t_in, time_step, r)
    % 计算 EV 充放电参数
    % 修正说明：明确区分时间步长 (time_step) 与调节时长 (delta_t_adj)
    
    % 计算式 (20) 的参数
    m_1 = -C_EV * r / (eta * time_step);  % 计算充放电速率的系数
    m_2 = C_EV * r / (eta * time_step);   % 计算充放电功率的系数
    m_3 = (E_tar - E_in) / (eta * (t_dep - t_in)); % 计算目标电量变化率
end
