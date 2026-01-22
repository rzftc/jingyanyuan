function [DeltaP_pu_plus_t, DeltaP_pu_minus_t] = calculateACAdjustmentPotentia(Pbase, Pmax, Pmin, alpha, beta, gamma, S_curr, delta_t_adj)
    % 输入参数说明:
    % Pbase       - 基线功率向量 [Nx1]
    % Pmax, Pmin  - 功率上下限向量 [Nx1]
    % alpha, beta, gamma - 系数向量 [Nx1]
    % S_curr      - 当前SOC向量 [Nx1]
    % delta_t_adj - 调节时长标量
    
    % 获取AC总数
    N = length(alpha);
    DeltaP_pu_plus_t = zeros(N, 1);
    DeltaP_pu_minus_t = zeros(N, 1);
    
    for i = 1:N
        % 1. 计算功率约束潜力 (式26-27)
        DeltaP_P_plus = Pmax(i) - Pbase(i);
        DeltaP_P_minus = Pmin(i) - Pbase(i);
        DeltaP_P_plus = max(DeltaP_P_plus, 0);
        DeltaP_P_minus = min(DeltaP_P_minus, 0);
        
        % 2. 计算能量约束潜力 (式30-31)
        numerator_plus = 1 - (1 - alpha(i) * delta_t_adj) * S_curr(i) - gamma(i) * delta_t_adj;
        denominator_plus = beta(i) * delta_t_adj;
        DeltaP_E_plus = numerator_plus / denominator_plus;
        
        numerator_minus = -(1 - alpha(i) * delta_t_adj) * S_curr(i) - gamma(i) * delta_t_adj;
        denominator_minus = beta(i) * delta_t_adj;
        DeltaP_E_minus = numerator_minus / denominator_minus;
        
        % 3. 修正能量约束边界 (式32-33)
        DeltaP_E_plus_adj = max(DeltaP_E_plus, 0);
        DeltaP_E_minus_adj = min(DeltaP_E_minus, 0);
        
        % 4. 综合功率与能量约束 (式34)
        DeltaP_pu_plus_t(i) = min(DeltaP_P_plus, DeltaP_E_plus_adj);
        DeltaP_pu_minus_t(i) = max(DeltaP_P_minus, DeltaP_E_minus_adj);
    end
end
