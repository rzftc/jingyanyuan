function S_prime = getSPrime(ev)
    % 根据substate状态计算sprime
    if strcmp(ev.substate, 'ON')
        S_prime = (ev.S_modified + 1) / 2;  % 开状态 (值域[0,1])
    elseif strcmp(ev.substate, 'OFF')
        S_prime = (ev.S_modified - 1) / 2;  % 关状态 (值域[-1,0])
    else
        error('未定义substate状态: %s', ev.substate);
    end
end
