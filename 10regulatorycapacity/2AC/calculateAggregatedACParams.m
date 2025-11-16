function AggParams = calculateAggregatedACParams(ACs_participating)
    % calculateAggregatedACParams: 计算空调聚合体的状态转移参数 A, B, C
    %
    % 实现了论文中的 (式 2-37)
    %
    % 输入:
    %   ACs_participating - 参与聚合的空调结构体数组。
    %                       要求每个元素 ACs(i) 必须包含
    %                       .alpha, .beta, .gamma (由 calculateACABC_single.m 预计算)
    %
    % 输出:
    %   AggParams         - 包含聚合参数的结构体:
    %                       .A, .B, .C (用于聚合状态转移)
    %                       .sum_term_1, .sum_term_2 (用于指令分解)

    if isempty(ACs_participating)
        warning('输入的空调列表为空，无法计算聚合参数。');
        AggParams = struct('A', 0, 'B', 0, 'C', 0, 'sum_term_1', 0, 'sum_term_2', 0);
        return;
    end

    num_ac = length(ACs_participating);
    
    % 预分配向量
    term_1_j = zeros(num_ac, 1); % (1 - alpha_j) / beta_j
    term_2_j = zeros(num_ac, 1); % gamma_j / beta_j
    term_3_j = zeros(num_ac, 1); % alpha_j / beta_j
    term_4_j = zeros(num_ac, 1); % 1 / beta_j

    % --- 1. 计算所有空调的中间项 ---
    for j = 1:num_ac
        beta_j = ACs_participating(j).beta;
        
        % 增加对 beta_j 接近零的保护
        if abs(beta_j) < 1e-9
            % 如果 beta 为 0 (例如 Tmax=Tmin)，该空调无法调节，不应参与聚合
            % 为避免除零，我们将这些项设为0或NaN，求和时忽略
            term_1_j(j) = NaN;
            term_2_j(j) = NaN;
            term_3_j(j) = NaN;
            term_4_j(j) = NaN;
        else
            alpha_j = ACs_participating(j).alpha;
            gamma_j = ACs_participating(j).gamma;
            
            term_1_j(j) = (1 - alpha_j) / beta_j; % 用于指令分解 (式 2-35)
            term_2_j(j) = gamma_j / beta_j;       % 用于指令分解 (式 2-35)
            term_3_j(j) = alpha_j / beta_j;       % 用于 A (式 2-37)
            term_4_j(j) = 1 / beta_j;             % 用于 A, B (式 2-37)
        end
    end

    % --- 2. 计算各项总和 (忽略NaN) ---
    sum_term_1 = sum(term_1_j, 'omitnan'); % sum( (1-alpha_j)/beta_j )
    sum_term_2 = sum(term_2_j, 'omitnan'); % sum( gamma_j/beta_j )
    sum_term_3 = sum(term_3_j, 'omitnan'); % sum( alpha_j/beta_j )
    sum_term_4 = sum(term_4_j, 'omitnan'); % sum( 1/beta_j )

    % --- 3. 计算聚合参数 A, B, C (式 2-37) ---
    if abs(sum_term_4) < 1e-9
        warning('聚合参数计算失败：sum(1/beta_j) 总和为零。');
        AggParams.A = 0;
        AggParams.B = 0;
        AggParams.C = 0;
    else
        AggParams.A = sum_term_3 / sum_term_4;
        AggParams.B = 1 / sum_term_4;
        AggParams.C = AggParams.B * sum_term_2;
    end
    
    % 存储用于指令分解的中间项
    AggParams.sum_term_1 = sum_term_1;
    AggParams.sum_term_2 = sum_term_2;
end