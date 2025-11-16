function Delta_P_j_commands = dispatchACOrderByStateConsistent(Delta_P_S, ACs_participating, AggParams)
    % dispatchACOrderByStateConsistent: 根据“状态一致”原则分解总功率指令
    %
    % 实现了论文中的 (式 2-35)
    %
    % 输入:
    %   Delta_P_S         - 聚合体接收到的总调节功率指令 (kW)
    %   ACs_participating - 参与聚合的空调结构体数组 (N x 1)
    %                       (必须包含 .alpha, .beta, .gamma)
    %   AggParams         - 从 calculateAggregatedACParams 获取的聚合参数
    %                       (必须包含 .sum_term_1, .sum_term_2)
    %
    % 输出:
    %   Delta_P_j_commands - 分配给每个空调的单独指令 (N x 1 向量)

    num_ac = length(ACs_participating);
    Delta_P_j_commands = zeros(num_ac, 1);

    % 检查聚合参数是否有效
    if abs(AggParams.sum_term_1) < 1e-9
        % 如果 sum((1-alpha)/beta) 为 0，无法分解，返回平均值
        warning('指令分解失败：sum_term_1 为零，采用平均分配作为后备。');
        Delta_P_j_commands(:) = Delta_P_S / num_ac;
        return;
    end

    % --- 1. 计算 (式 2-35) 中的公共因子 ---
    % common_factor = ( Delta_P_S + Σ(γ_j/β_j) ) / ( Σ((1-α_j)/β_j) )
    common_factor = (Delta_P_S + AggParams.sum_term_2) / AggParams.sum_term_1;

    % --- 2. 计算每个单体的指令 (式 2-35) ---
    for j = 1:num_ac
        beta_j = ACs_participating(j).beta;
        
        if abs(beta_j) < 1e-9
            % 该空调无法调节，指令为 0
            Delta_P_j_commands(j) = 0;
        else
            alpha_j = ACs_participating(j).alpha;
            gamma_j = ACs_participating(j).gamma;
            
            % 计算该空调独有的项
            term_1_j = (1 - alpha_j) / beta_j;
            term_2_j = gamma_j / beta_j;
            
            % 应用 (式 2-35)
            Delta_P_j_commands(j) = term_1_j * common_factor - term_2_j;
        end
    end
end