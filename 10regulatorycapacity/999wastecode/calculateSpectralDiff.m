function [spec_plus, spec_minus, Constraints] = calculateSpectralDiff(n_AC, n_EV, AC_plus, AC_minus, EV_plus, EV_minus)
    L = 3; % 取前3个主频分量
    T = size(AC_plus,1);
    
    % ===== 上调节频谱差异 =====
    % 预计算角度矩阵
    theta = 2*pi*(0:L-1)'*(0:T-1)/T;
    cos_mat = cos(theta);
    sin_mat = -sin(theta); % 注意负号
    
    % 定义频域分量
    a_AC_plus = sdpvar(L,1);
    b_AC_plus = sdpvar(L,1);
    a_EV_plus = sdpvar(L,1);
    b_EV_plus = sdpvar(L,1);
    
    % 空调分量约束
    AC_plus_flat = sum(n_AC.*AC_plus,2); % 按时间维度求和
    EV_plus_flat = sum(n_EV.*EV_plus,2);
    Constraints = [];
    for k = 1:L
        % 实部计算
        Constraints = [Constraints;
            a_AC_plus(k) == sum(AC_plus_flat' .* cos_mat(k,:))];
        Constraints = [Constraints;
            a_EV_plus(k) == sum(EV_plus_flat' .* cos_mat(k,:))];
        
        % 虚部计算
        Constraints = [Constraints;
            b_AC_plus(k) == sum(AC_plus_flat' .* sin_mat(k,:))];
        Constraints = [Constraints;
            b_EV_plus(k) == sum(EV_plus_flat' .* sin_mat(k,:))];
    end
    
    % 频谱差异计算
    spec_plus = 0;
    for k = 1:L
        % 定义幅值差辅助变量
        d_plus = sdpvar(1);
        
        % 二阶锥约束 |(a_AC - a_EV, b_AC - b_EV)| <= d_plus
        Constraints = [Constraints;
            cone([(a_AC_plus(k)-a_EV_plus(k)); (b_AC_plus(k)-b_EV_plus(k))], d_plus)];
        
        spec_plus = spec_plus + d_plus^2;
    end
    
    % ===== 下调节频谱差异 =====
    AC_minus_flat = sum(n_AC.*abs(AC_minus),2);
    EV_minus_flat = sum(n_EV.*abs(EV_minus),2);
    
    a_AC_minus = sdpvar(L,1);
    b_AC_minus = sdpvar(L,1);
    a_EV_minus = sdpvar(L,1);
    b_EV_minus = sdpvar(L,1);
    
    for k = 1:L
        Constraints = [Constraints;
            a_AC_minus(k) == sum(AC_minus_flat' .* cos_mat(k,:))];
        Constraints = [Constraints;
            a_EV_minus(k) == sum(EV_minus_flat' .* cos_mat(k,:))];
        Constraints = [Constraints;
            b_AC_minus(k) == sum(AC_minus_flat' .* sin_mat(k,:))];
        Constraints = [Constraints;
            b_EV_minus(k) == sum(EV_minus_flat' .* sin_mat(k,:))];
    end
    
    spec_minus = 0;
    for k = 1:L
        d_minus = sdpvar(1);
        Constraints = [Constraints;
            cone([(a_AC_minus(k)-a_EV_minus(k)); (b_AC_minus(k)-b_EV_minus(k))], d_minus)];
        spec_minus = spec_minus + d_minus^2;
    end
end
