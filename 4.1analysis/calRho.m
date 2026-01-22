function [rho] = calRho(n_AC, deltaP_AC, n_EV, deltaP_EV)
% 计算斯皮尔曼秩相关系数（含设备数量参数）
% 输入参数：
%   n_AC      - [T×1] 各时段空调参与数量
%   deltaP_AC - [T×1] 单台空调调节能力
%   n_EV      - [T×1] 各时段EV参与数量
%   deltaP_EV - [T×1] 单台EV调节能力
% 输出：
%   rho       - 斯皮尔曼秩相关系数[-1,1]

    % ===== 步骤1：生成调节能力序列 =====
    AC_series = n_AC .* deltaP_AC; % 空调总调节能力
    EV_series = n_EV .* deltaP_EV; % EV总调节能力
    
    % ===== 步骤2：数据校验 =====
    assert(length(AC_series)==length(EV_series), '时间序列长度不一致');
    T = length(AC_series);
    if T < 2
        rho = NaN;
        warning('样本量不足，无法计算相关系数');
        return;
    end
    
    % ===== 步骤3：计算独立排名（处理并列值）=====
    rank_AC = tiedrank(AC_series); % MATLAB内置函数处理并列排名
    rank_EV = tiedrank(EV_series);
    
    % ===== 步骤4：计算排名差平方和 =====
    d_sq = (rank_AC - rank_EV).^2;
    sum_d_sq = sum(d_sq);
    
    % ===== 步骤5：计算相关系数 =====
    rho = 1 - (6 * sum_d_sq) / (T * (T^2 - 1));
    
    % ===== 数值边界处理 =====
    rho = max(min(rho, 1), -1); % 强制限制在[-1,1]
end