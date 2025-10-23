function [alpha, beta, gamma] = calculateACABC_single(R, C, eta, Tmax, Tmin, Tset, deltaT)
    %CALCULATEACABC_SINGLE 计算单个AC的 alpha, beta, gamma 系数
    %
    % 输入参数：
    %   R       - 电阻（标量）
    %   C       - 电容（标量）
    %   eta     - 效率（标量）
    %   Tmax    - 最大温度（标量）
    %   Tmin    - 最小温度（标量）
    %   Tset    - 设定温度（标量）
    %   deltaT  - 时间间隔（标量）
    %
    % 输出参数：
    %   alpha   - 衰减系数 (exp(-deltaT/(R*C)))
    %   beta    - 功率转换系数
    %   gamma   - 温度偏差系数

    % 计算衰减系数
    alpha = exp(-deltaT / (R * C));
    
    % 计算功率转换系数
    beta = (eta * R * (1 - exp(-deltaT / (R * C)))) / (Tmax - Tmin);
    
    % 计算温度偏差系数（强制按元素操作）
    gamma = ((Tmax - Tset) .* (1 - exp(-deltaT / (R * C)))) ./ (Tmax - Tmin);
end
