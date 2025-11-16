function SOC = calculateACS_single_corrected(T_j_current, Tmax, Tmin)
    % calculateACS_single_corrected: 根据【室内】温度计算单个AC的 SOC
    %
    % 实现了论文中的 (式 2-8): Soc_j(t) = (T_jmax - T_j(t)) / (T_jmax - T_jmin)
    %
    % 输入:
    %   T_j_current - 空调【当前】的室内温度 (T_j(t))
    %   Tmax        - 用户设定的舒适温度上限
    %   Tmin        - 用户设定的舒适温度下限
    %
    % 输出:
    %   SOC         - 计算得到的当前SOC值 [0, 1]

    % 应用 (式 2-8)
    SOC_raw = (Tmax - T_j_current) / (Tmax - Tmin);
    
    % 确保SOC在 [0, 1] 范围内
    SOC = max(0, min(1, SOC_raw));
end