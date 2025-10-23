function P_baseline = ACbaseP_single(T_ja, T_jset, R, eta)
    %ACBASEP_SINGLE 计算单个空调系统的基准功率
    %
    % 输入参数：
    %   T_ja    : 当前环境温度（标量，单位：℃）
    %   T_jset  : 目标设定温度（标量，单位：℃）
    %   R       : 系统热阻值（标量，单位：℃/kW）
    %   eta     : 系统热效率（标量，无单位）
    %
    % 输出参数：
    %   P_baseline : 单体系统基准功率（标量，单位：kW）
    %
    % 计算公式：
    %   P_baseline = (T_ja - T_jset) / (R * eta)
    
    % 执行功率计算
    P_baseline = (T_ja - T_jset) / (R * eta);
end
