function [cov_plus, cov_minus] = calculateSpearmanCorrelation(AC_plus, AC_minus, EV_plus, EV_minus)
    % 输入向量化处理
    x_plus = AC_plus(:); y_plus = EV_plus(:);
    x_minus = AC_minus(:); y_minus = EV_minus(:);
    
    % 动态计算均值（式0-21, 0-22）
    mu_x_plus = mean(x_plus); mu_y_plus = mean(y_plus);
    mu_x_minus = mean(x_minus); mu_y_minus = mean(y_minus);
    
    % 计算中心化协方差和
    cov_plus = sum( (x_plus - mu_x_plus) .* (y_plus - mu_y_plus) );
    cov_minus = sum( (x_minus - mu_x_minus) .* (y_minus - mu_y_minus) );
end
