function [P_agg, S_agg_next] = calculateVirtualSOC_agg(EVs, dt)
    % 筛选可用EV（ON/LockON状态）
    active_mask = strcmp({EVs.state}, 'ON') | strcmp({EVs.state}, 'LockON');
    active_EVs = EVs(active_mask);
    
    P_agg = 0;
    P_req = 0;
    E_exp = 0;
    E_actual = 0;
    C = 0;
    r = 0;
    
    % 计算所有可用EV的当前功率总和
    for i = 1:length(active_EVs) 
        P_agg = P_agg + active_EVs(i).P_current;
        
        % 计算期望电量（公式22-23）
        P_req = (active_EVs(i).E_tar - active_EVs(i).E_ini) / (active_EVs(i).eta * (active_EVs(i).t_dep - active_EVs(i).t_in) / 60);
        E_exp = E_exp + active_EVs(i).E_exp + active_EVs(i).eta * P_req * (dt / 60);
        
        % 实际电量变化（假设理想充电）
        E_actual = E_actual + active_EVs(i).E_actual + active_EVs(i).eta * active_EVs(i).P_current * (dt / 60);
        
        % 累加电池容量和平均电阻
        C = C + active_EVs(i).C;
        r = mean(active_EVs(i).r);
    end
    
    % 根据聚合系数计算下一个时间步的虚拟SOC
    S_agg_next = -(E_actual - E_exp) / (C * r);
end
