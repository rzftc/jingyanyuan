function [N_pred, models_params] = grey_prediction_gm11(N_hist, predict_years)
    % GREY_PREDICTION_GM11 灰色预测模型 GM(1,1)
    % 输入:
    %   N_hist: 历史保有量数据向量 (行或列向量)
    %   predict_years: 需要预测的未来年数
    % 输出:
    %   N_pred: 预测的未来保有量向量 (1 x predict_years)
    %   models_params: 结构体，包含发展系数 a 和灰作用量 b

    x0 = N_hist(:); % 确保为列向量
    n = length(x0);
    
    % 1. 级比检验 (可选，此处略过直接建模)
    
    % 2. 一次累加生成 (1-AGO)
    x1 = cumsum(x0);
    
    % 3. 构造数据矩阵 B 和数据向量 Y
    % 紧邻均值生成序列
    z1 = 0.5 * (x1(1:end-1) + x1(2:end));
    B = [-z1, ones(n-1, 1)];
    Y = x0(2:end);
    
    % 4. 最小二乘法求解参数 a, b
    % u = [a; b] = inv(B'*B) * B' * Y
    try
        u = (B'*B) \ (B'*Y);
    catch
        warning('GM(1,1) 矩阵求逆失败，可能数据共线性太高。');
        N_pred = nan(1, predict_years);
        models_params = struct('a', NaN, 'b', NaN);
        return;
    end
    a = u(1);
    b = u(2);
    models_params.a = a;
    models_params.b = b;
    
    % 5. 建立时间响应函数并预测
    % x1_pred(k+1) = (x0(1) - b/a) * exp(-a*k) + b/a
    N_pred_cum = zeros(n + predict_years, 1);
    N_pred_cum(1) = x0(1);
    for k = 1:(n + predict_years - 1)
        N_pred_cum(k+1) = (x0(1) - b/a) * exp(-a * k) + b/a;
    end
    
    % 6. 累减还原得到原始数据预测值
    N_pred_all = [N_pred_cum(1); diff(N_pred_cum)];
    
    % 提取未来预测部分
    N_pred = N_pred_all(n+1 : end)';
end