function N_pred = svr_autoregressive_prediction(N_hist, X_hist, X_pred)
    % SVR_AUTOREGRESSIVE_PREDICTION (Autoregressive Version)
    % 使用支持向量回归 (SVR) 预测增长量 Delta_N(t) = f(N(t-1), X(t))
    %
    % 输入:
    %   N_hist: 历史保有量数据向量 (行或列向量)
    %   X_hist: 历史年份对应的外部影响因子矩阵 (行数必须等于 length(N_hist))
    %   X_pred: 未来年份对应的外部影响因子矩阵
    %
    % 输出:
    %   N_pred: 预测的未来保有量向量 (列向量)

    % 检查 Statistics and Machine Learning Toolbox
    if ~license('test', 'Statistics_Toolbox') || isempty(which('fitrsvm'))
        warning('未检测到 Statistics and Machine Learning Toolbox 或 fitrsvm 函数。SVR 预测无法执行。将返回 NaN。');
        N_pred = nan(size(X_pred, 1), 1); 
        return;
    end
    
    N_hist = double(N_hist(:));
    X_hist = double(X_hist);
    X_pred = double(X_pred);
    n_hist_samples = length(N_hist);
    n_pred_steps = size(X_pred, 1);

    % --- 1. 准备训练数据 (自回归) ---
    % 目标 Y = Delta_N(t) = N(t) - N(t-1)
    Y_train = diff(N_hist); % (n-1) x 1
    
    % 特征 X = [N(t-1), X_i(t)]
    if size(X_hist, 1) ~= n_hist_samples
        error('AR-SVR: X_hist 行数 (%d) 与 N_hist 长度 (%d) 必须匹配。', size(X_hist, 1), n_hist_samples);
    end
    % X_train 组合: [N(t-1), X_i(t)]
    % 即 N_hist(1:end-1) 对应 X_hist(2:end, :)
    X_train = [N_hist(1:end-1), X_hist(2:end, :)];
    
    if isempty(X_train) || isempty(Y_train) || length(Y_train) < 2
         warning('AR-SVR: 历史数据不足 (%d) 无法创建有效的SVR训练集。将返回 NaN。', n_hist_samples);
         N_pred = nan(n_pred_steps, 1); 
         return;
    end
    
    try
        % --- 2. 训练 SVR 模型 ---
        % 'KernelFunction', 'gaussian' (即RBF) 适用于非线性
        % 'Standardize', true 是一个好的实践
        fprintf('  [AR-SVR] 正在训练 SVR (高斯核) 模型...\n');
        Mdl = fitrsvm(X_train, Y_train, ...
            'KernelFunction', 'gaussian', ...
            'Standardize', true, ...
            'KernelScale', 'auto');
        
        % --- 3. 逐年迭代预测 ---
        N_pred = zeros(n_pred_steps, 1);
        N_current = N_hist(end); % 启动迭代的最后一个历史值
        
        for t = 1:n_pred_steps
            % 构造当前步的输入: [N(t-1), X_i(t)]
            X_in_t = [N_current, X_pred(t, :)];
            
            % 预测 *增长量*
            Delta_N_t = predict(Mdl, X_in_t);
            
            % (安全检查：防止因拟合不佳导致负增长)
            if N_current + Delta_N_t < N_current
                 Delta_N_t = 0; % 假设保有量至少不减少
            end
            
            % 计算 *保有量*
            N_current = N_current + Delta_N_t;
            N_pred(t) = N_current;
        end
        
    catch ME
        warning('自回归 SVR 模型训练或预测失败: %s', ME.message);
        N_pred = nan(n_pred_steps, 1);
    end
end