function [N_pred, fit_info] = arima_prediction(N_hist, predict_years, order)
    % ARIMA_PREDICTION 自回归积分滑动平均模型
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
        fit_info = [];
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
        warning('ARIMA 模型拟合或预测失败: %s');
        N_pred = nan(1, predict_years);
        fit_info = ME;
    end
end