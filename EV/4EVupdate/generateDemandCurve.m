%% generateDemandCurve.m %%
function EV = generateDemandCurve(EV)
    % 需求曲线定义（严格按原文数学表达式实现）：
    %   D(s) = P_ON 当 -1 ≤ s < s_prime 
    %   D(s) = 0    当 s_prime ≤ s ≤ 1
    switch EV.state
        case 'LockON'
            EV.demandCurve = @(s) EV.P_N;  % 强制输出P_N
            
        case 'LockOFF'
            EV.demandCurve = @(s) 0;       % 强制输出0
        case 'OFF'   
            EV.demandCurve = @(s) 0; 
        otherwise
            if strcmp(EV.substate, 'ON')
                % ON状态: s_prime = (S_modified + 1)/2 ∈ [0,1]
                s_prime = (EV.S_modified + 1) / 2;
                EV.demandCurve = @(s) EV.P_N .* (-1 <= s & s < s_prime);
                
            elseif strcmp(EV.substate, 'OFF')
                % OFF状态: s_prime = (S_modified - 1)/2 ∈ [-1,0]
                s_prime = (EV.S_modified - 1) / 2;
                EV.demandCurve = @(s) EV.P_N .* (-1 <= s & s < s_prime);
                
            else
                error('未知substate状态: %s', EV.substate);
            end
    end
    
    % 防御性处理s超出[-1,1]的情况
    % EV.demandCurve = @(s) EV.demandCurve(s) .* (s >= -1 & s <= 1);
end
