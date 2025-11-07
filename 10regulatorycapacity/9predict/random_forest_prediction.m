%% --- [修正] 随机森林预测函数 (自回归版本) ---
function N_pred = random_forest_prediction(N_hist, X_hist, X_pred, n_trees)
    % RANDOM_FOREST_PREDICTION (Autoregressive Version)
    % 预测增长量 Delta_N(t) = f(N(t-1), X(t))
    
    if nargin < 4 || isempty(n_trees), n_trees = 50; end
    
    % 检查 Statistics and Machine Learning Toolbox
    if ~license('test', 'Statistics_Toolbox') || isempty(which('TreeBagger'))
        warning('未检测到 Statistics and Machine Learning Toolbox 或 TreeBagger 函数。随机森林预测无法执行。将返回 NaN。');
        N_pred = nan(size(X_pred, 1), 1); return;
    end
    
    N_hist = double(N_hist(:));
    X_hist = double(X_hist);
    X_pred = double(X_pred);
    n_hist_samples = length(N_hist);
    n_pred_steps = size(X_pred, 1);

    % --- 1. 准备训练数据 ---
    % 目标 Y = Delta_N(t) = N(t) - N(t-1)
    Y_train = diff(N_hist); % (n-1) x 1
    
    % 特征 X = [N(t-1), X_i(t)]
    % 检查 X_hist 是否包含所有年份的数据 (与 N_hist 长度相同)
    if size(X_hist, 1) ~= n_hist_samples
        error('RF-AR: X_hist 行数 (%d) 与 N_hist 长度 (%d) 必须匹配。', size(X_hist, 1), n_hist_samples);
    end
    % X_train 组合: [N(t-1), X_i(t)]
    % 即 N_hist(1:end-1) 对应 X_hist(2:end, :)
    X_train = [N_hist(1:end-1), X_hist(2:end, :)];
    
    if isempty(X_train) || isempty(Y_train)
         warning('RF-AR: 历史数据不足 (%d) 无法创建训练集。', n_hist_samples);
         N_pred = nan(n_pred_steps, 1); return;
    end
    
    try
        % 2. 训练模型
        fprintf('  [RF-AR] 正在训练 %d 棵树的自回归随机森林模型...\n', n_trees);
        Mdl = TreeBagger(n_trees, X_train, Y_train, 'Method', 'regression', 'OOBPrediction', 'off');
        
        % 3. 逐年迭代预测
        N_pred = zeros(n_pred_steps, 1);
        N_current = N_hist(end); % 启动迭代的最后一个历史值
        
        for t = 1:n_pred_steps
            % 构造当前步的输入: [N(t-1), X_i(t)]
            X_in_t = [N_current, X_pred(t, :)];
            
            % 预测 *增长量*
            Delta_N_t = predict(Mdl, X_in_t);
            
            % 计算 *保有量*
            N_current = N_current + Delta_N_t;
            N_pred(t) = N_current;
        end
        
    catch ME
        warning('自回归随机森林模型训练或预测失败: %s');
        N_pred = nan(n_pred_steps, 1);
    end
end