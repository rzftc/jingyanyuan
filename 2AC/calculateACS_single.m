function SOC = calculateACS_single(T_ja, Tmax, Tmin)
    % calculateSOC_single: 计算单个AC在下一时间步的 SOC
    % 输入：
    %   alpha      - 系数 alpha_j（标量）
    %   beta       - 系数 beta_j（标量）
    %   gamma      - 系数 gamma_j（标量）
    %   SOC_current - 当前时间步的 SOC（标量）
    %   delta_P    - 功率变化量 Delta P_j(t)（标量）
    %
    % 输出：
    %   SOC_next   - 下一时间步的 SOC（标量）
    
    % 计算下一个时间步的 SOC
    SOC = (Tmax-T_ja)/(Tmax-Tmin);
end
