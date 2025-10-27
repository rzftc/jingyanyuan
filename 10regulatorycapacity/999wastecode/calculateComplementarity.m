function [SDCI_plus, SDCI_minus, ICI] = calculateComplementarity(AC_plus, AC_minus, EV_plus, EV_minus)
    % 计算空调与电动汽车调节能力的互补性指标
    % 输入:
    %   AC_plus  - 空调上调节能力时间序列 [T×1]
    %   AC_minus - 空调下调节能力时间序列 [T×1]
    %   EV_plus  - 电动汽车上调节能力时间序列 [T×1]
    %   EV_minus - 电动汽车下调节能力时间序列 [T×1]
    % 输出:
    %   SDCI_plus  - 上调节互补性指数
    %   SDCI_minus - 下调节互补性指数
    %   ICI        - 综合互补性指数
    AC_plus = AC_plus(:); EV_plus = EV_plus(:);
    AC_minus = AC_minus(:); EV_minus = EV_minus(:);
    % 计算上调节互补性指数 (式4-29)
    numerator_plus = sum(min(AC_plus, EV_plus));
    denominator_plus = sum(max(AC_plus, EV_plus));
    SDCI_plus = numerator_plus / denominator_plus;
    
    % 计算下调节互补性指数 (式4-30)
    numerator_minus = sum(min(abs(AC_minus), abs(EV_minus)));
    denominator_minus = sum(max(abs(AC_minus), abs(EV_minus)));
    SDCI_minus = numerator_minus / denominator_minus;
    
    % 计算综合互补性指数 (式4-31)
    ICI = (SDCI_plus + SDCI_minus) / 2;
end