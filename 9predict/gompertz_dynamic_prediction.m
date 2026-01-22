function N_pred = gompertz_dynamic_prediction(N_hist, predict_years, K_dynamic)
    % GOMPERTZ_DYNAMIC_PREDICTION 动态饱和水平 Gompertz 模型
    % 模型公式: N(t) = K(t) * exp(-b * exp(-r * t))
    % 输入:
    %   N_hist: 历史保有量
    %   predict_years: 预测年数
    %   K_dynamic: 动态饱和水平向量。
    %              如果是标量，则为固定饱和水平。
    %              如果是向量，长度应为 (length(N_hist) + predict_years)，覆盖历史和未来。
    % 输出:
    %   N_pred: 预测值向量

    t_hist = (1:length(N_hist))';
    y_hist = N_hist(:);
    
    % 处理动态 K
    len_total = length(N_hist) + predict_years;
    if isscalar(K_dynamic)
        K_vec = repmat(K_dynamic, len_total, 1);
    elseif length(K_dynamic) >= len_total
        K_vec = K_dynamic(1:len_total);
        K_vec = K_vec(:);
    else
        error('K_dynamic 向量长度不足，至少需要覆盖历史+预测年份 (%d 年)', len_total);
    end
    
    K_hist = K_vec(1:length(N_hist));
    
    % 参数估计 (线性化方法)
    % ln(N(t)/K(t)) = -b * exp(-r * t)
    % ln(-ln(N(t)/K(t))) = ln(b) - r * t
    % 令 Y = ln(-ln(N(t)/K(t))), A = ln(b), B = -r
    % Y = A + B * t
    
    % 确保历史数据小于当前的 K，否则对数无意义
    valid_idx = y_hist < K_hist & y_hist > 0;
    if sum(valid_idx) < 3
         warning('Gompertz 模型有效历史数据点不足（需 N(t) < K(t)），无法拟合。');
         N_pred = nan(1, predict_years);
         return;
    end
    
    try
        Y_trans = log(-log(y_hist(valid_idx) ./ K_hist(valid_idx)));
        t_valid = t_hist(valid_idx);
        
        % 线性回归拟合 Y = A + B*t
        X_reg = [ones(length(t_valid), 1), t_valid];
        coeffs = (X_reg' * X_reg) \ (X_reg' * Y_trans);
        A = coeffs(1);
        B = coeffs(2);
        
        b_est = exp(A);
        r_est = -B;
        
        % 预测未来
        t_pred = (length(N_hist) + 1 : length(N_hist) + predict_years)';
        K_future = K_vec(length(N_hist)+1 : end);
        
        % N(t) = K(t) * exp(-b * exp(-r * t))
        N_pred_col = K_future .* exp(-b_est .* exp(-r_est .* t_pred));
        N_pred = N_pred_col';
        
    catch ME
        warning('Gompertz 模型拟合失败: %s', ME.message);
        N_pred = nan(1, predict_years);
    end
end