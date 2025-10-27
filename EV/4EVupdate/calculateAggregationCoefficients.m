%% 新版参数聚合函数（集成功率边界计算）
function [M1, M2, M3, P_agg_max] = calculateAggregationCoefficients(EVs, dt)
    % 筛选可用EV（ON/LockON状态）
    active_mask = strcmp({EVs.state}, 'ON') | strcmp({EVs.state}, 'LockON');
    active_EVs = EVs(active_mask);
    
    % 初始化聚合系数
    M1 = 0;
    M2 = 0;
    M3 = 0;
    
    % 遍历所有激活的电动汽车，计算聚合系数
    for i = 1:length(active_EVs) % 修改了循环变量的范围
        M1 = M1 + (-active_EVs(i).C) / (active_EVs(i).eta * dt);
        M2 = -M1; % M2 是 M1 的负值
        M3 = M3 + (active_EVs(i).E_tar - active_EVs(i).E_ini) / (active_EVs(i).eta * (active_EVs(i).t_dep - active_EVs(i).t_in)/60);
    end
    
    % 关键修正：计算实时最大功率
    P_agg_max = sum([EVs.P_N]); % 所有激活EV的额定功率之和
end
