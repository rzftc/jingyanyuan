%% calculateDeltaE_vec.m (修正版)
function Delta_E_vec = calculateDeltaE_vec(P_real_vec, P_h_max_vec, P_0_vec, P_l_min_vec, Delta_E_h_max_vec, Delta_E_q_max_vec)
% CALCULATEDELTAE_VEC 向量化计算 Delta_E. (修正数组维度不兼容问题)
%   输入均为列向量.

    num_evs = length(P_real_vec);
    Delta_E_vec = zeros(num_evs, 1);

    % 条件掩码 (Logical Masks)
    mask1 = (P_real_vec >= P_h_max_vec);
    mask4 = (P_real_vec < P_l_min_vec); % 注意这里用 <
    mask2 = (P_real_vec >= P_0_vec) & ~mask1; % P_0 <= P_real < P_h_max
    mask3 = (P_real_vec >= P_l_min_vec) & ~mask2 & ~mask1; % P_l_min <= P_real < P_0

    % --- 条件 1 ---
    Delta_E_vec(mask1) = Delta_E_h_max_vec(mask1);

    % --- 条件 4 ---
    Delta_E_vec(mask4) = Delta_E_q_max_vec(mask4); % 注意是 Delta_E_q_max

    % --- 条件 2 ---
    if any(mask2)
        % 提取 mask2 对应的数据子集
        delta_p2 = P_real_vec(mask2) - P_0_vec(mask2);
        P_h_diff2 = P_h_max_vec(mask2) - P_0_vec(mask2);
        Delta_E_h_max_vec_m2 = Delta_E_h_max_vec(mask2); % 先提取子集

        % 避免除零的掩码 (作用于子集上)
        valid_h_diff = P_h_diff2 ~= 0;

        term1 = zeros(size(delta_p2));
        term2 = zeros(size(delta_p2));

        if any(valid_h_diff)
           % --- !!! 修正索引 !!! ---
           term1(valid_h_diff) = -(Delta_E_h_max_vec_m2(valid_h_diff) .* delta_p2(valid_h_diff).^2) ./ (P_h_diff2(valid_h_diff).^2);
           term2(valid_h_diff) = (2 * Delta_E_h_max_vec_m2(valid_h_diff) .* delta_p2(valid_h_diff)) ./ P_h_diff2(valid_h_diff);
           % --- 修正结束 ---
        end
        Delta_E_vec(mask2) = term1 + term2; % 将子集结果写回
    end

    % --- 条件 3 ---
    if any(mask3)
        % 提取 mask3 对应的数据子集
        delta_p3 = P_real_vec(mask3) - P_0_vec(mask3); % 注意 delta_p 定义
        P_l_diff3 = P_0_vec(mask3) - P_l_min_vec(mask3);
        Delta_E_q_max_vec_m3 = Delta_E_q_max_vec(mask3); % 先提取子集

        % 避免除零的掩码 (作用于子集上)
        valid_l_diff = P_l_diff3 ~= 0;

        numerator = zeros(size(delta_p3));
        denominator = ones(size(delta_p3)); % 分母设为1避免下面除零警告

        if any(valid_l_diff)
            numerator(valid_l_diff) = delta_p3(valid_l_diff).^2 .* (3 * P_l_diff3(valid_l_diff) + 2 * delta_p3(valid_l_diff));
            denominator(valid_l_diff) = P_l_diff3(valid_l_diff).^3;

            % 再次检查分母是否为零（如果 P_0 == P_l_min）
            zero_denom_mask = denominator == 0;
            if any(zero_denom_mask)
                warning('calculateDeltaE_vec:ZeroDenominator', '条件3计算中出现零分母 (P_0 == P_l_min)，Delta_E 将设为0');
                numerator(zero_denom_mask) = 0; % 避免 NaN
                denominator(zero_denom_mask) = 1; % 避免 NaN
            end
        else
             warning('calculateDeltaE_vec:ZeroDenominator', '条件3计算中 P_0 - P_l_min 全为零');
             % 如果所有 P_0 - P_l_min 都为零，确保分子也为零
             numerator(:) = 0;
             denominator(:) = 1;
        end

        % --- !!! 修正索引 !!! ---
        % Delta_E_q_max_vec_m3 已经是对应 mask3 的子集，直接使用
        Delta_E_vec(mask3) = Delta_E_q_max_vec_m3 .* numerator ./ denominator;
        % --- 修正结束 ---
    end

end