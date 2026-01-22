function SOC_next = updateACSOC_single(SOC_current, Delta_Pj, alpha, beta, gamma)
    % updateACSOC_single: 根据功率指令更新单个AC的SOC状态
    %
    % 该函数实现了理论模型中的空调状态转移方程 (3-10):
    %   soc_j(t+1) = alpha_j * soc_j(t) + beta_j * Delta_P_j(t) + gamma_j
    %
    % 输入:
    %   SOC_current - 当前时间步的SOC (soc_j(t))
    %   Delta_Pj    - 应用于当前时间步的调节功率指令 (Delta_P_j(t)) 
    %   alpha       - 系数 alpha_j (来自 calculateACABC_single) [cite: 32]
    %   beta        - 系数 beta_j (来自 calculateACABC_single) [cite: 32]
    %   gamma       - 系数 gamma_j (来自 calculateACABC_single) [cite: 32]
    %
    % 输出:
    %   SOC_next    - 下一个时间步的SOC (soc_j(t+1))，已约束在 [0, 1] 范围

    % 1. 应用状态转移方程 (3-10)
    SOC_next_raw = alpha * SOC_current + beta * Delta_Pj + gamma; 
    
    % 2. 应用SOC边界约束 (3-9)
    % 确保SOC值始终保持在 [0, 1] 范围内
    SOC_next = max(0, min(1, SOC_next_raw)); 

end