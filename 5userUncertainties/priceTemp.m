function T_set = priceTemp(p_real, p_min, p_max, p_0, T_max)
    % 修正说明：按论文式(22)-(23)重构为线性模型（原错误：p < p_min时使用三次函数）
    
    T_set = zeros(size(p_real));
    
    % Case 1: p >= p_max -> T_max
    T_set(p_real >= p_max) = T_max;
    
    % Case 2: p_min <= p < p_max -> 线性响应
    idx = (p_real >= p_min) & (p_real < p_max);
    T_set(idx) = (T_max / (p_max - p_min)) * (p_real(idx) - p_min);
    
    % Case 3: p < p_min -> 无响应（设为0）
    T_set(p_real < p_min) = 0;
end
