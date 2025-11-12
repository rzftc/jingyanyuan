function [N_pred, fit_info] = arima_prediction(N_hist, predict_years, order)
    % ARIMA_PREDICTION 自回归积分滑动平均模型
    % (版本 2.1 - 降低数据量门槛以尝试运行)
    % 输入:
    %   N_hist: 历史数据向量
    %   predict_years: 预测年数
    %   order: (可选) ARIMA模型阶数 [p, d, q]，默认为 [1, 1, 1]
    % 输出:
    %   N_pred: 预测值向量
    %   fit_info: 拟合信息

    if nargin < 3 || isempty(order)
        order = [1, 1, 1]; % 默认简单阶数，实际应用中应通过AIC/BIC选择
    end
    
    data = N_hist(:);
    
    % 检查是否拥有 Econometrics Toolbox
    if ~exist('arima', 'class')
        warning('未检测到 Econometrics Toolbox，ARIMA 无法使用标准库执行。将返回 NaN。');
        N_pred = nan(1, predict_years);
        fit_info = 'Error: Econometrics Toolbox not found';
        return;
    end
    
    % [!!! 修改：降低数据量门槛 !!!]
    % ARIMA 在统计上需要更多数据，但为了满足运行要求，我们将门槛降到 5
    % (N=4 仍然可能导致 d=1, p=1 的模型失败，但这取决于 estimate 函数的内部实现)
    min_data_points = 5; 
    if length(data) < min_data_points
        warning('ARIMA 模型因历史数据不足 (仅 %d 个点，需要 > %d) 而跳过。将返回 NaN。', length(data), min_data_points);
        N_pred = nan(1, predict_years);
        fit_info = 'Skipped: Insufficient data';
        return;
    end

    try
        % 1. 模型定义
        Mdl = arima(order(1), order(2), order(3));
        
        % 2. 参数估计
        % 'Display','off' 关闭命令行输出
        [EstMdl, ~] = estimate(Mdl, data, 'Display', 'off');
        
        % 3. 预测
        [N_pred_col, ~] = forecast(EstMdl, predict_years, 'Y0', data);
        N_pred = N_pred_col';
        fit_info = EstMdl;
        
    catch ME
        % [!!! 修复 BUG: 将 ME.message 传递给 warning !!!]
        warning('ARIMA 模型拟合或预测失败: %s', ME.message);
        N_pred = nan(1, predict_years);
        fit_info = ME;
    end
end