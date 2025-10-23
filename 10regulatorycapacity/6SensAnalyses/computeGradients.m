function [dP_up, dP_down] = computeGradients(p, up, down)
    % 输入维度验证
    assert(isrow(p), 'p必须是行向量');
    [T, P] = size(up);
    dP_up = zeros(T, P);
    dP_down = zeros(T, P);
    
    % 中心差分法
    for i = 1:P
        if i == 1
            h = p(2) - p(1);
            dP_up(:,i) = (up(:,2) - up(:,1)) / h;
            dP_down(:,i) = (down(:,2) - down(:,1)) / h;
        elseif i == P
            h = p(end) - p(end-1);
            dP_up(:,i) = (up(:,end) - up(:,end-1)) / h;
            dP_down(:,i) = (down(:,end) - down(:,end-1)) / h;
        else
            h_forward = p(i+1) - p(i);
            h_backward = p(i) - p(i-1);
            dP_up(:,i) = (up(:,i+1)*h_backward^2 - up(:,i-1)*h_forward^2 - ...
                up(:,i)*(h_backward^2 - h_forward^2)) / (h_backward*h_forward*(h_backward + h_forward));
            dP_down(:,i) = (down(:,i+1)*h_backward^2 - down(:,i-1)*h_forward^2 - ...
                down(:,i)*(h_backward^2 - h_forward^2)) / (h_backward*h_forward*(h_backward + h_forward));
        end
    end
end
