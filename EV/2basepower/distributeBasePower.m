function EVs = distributeBasePower(EVs, ~, P_base_opt)
    % 将P_base_opt序列按照当前时刻参与聚合的EV的额定功率加权分配
    % 输入：
    %   EVs - 电动汽车结构体数组
    %   current_time - 当前时间(分钟)
    %   P_base_opt - 待分配的基准功率序列 (N x 1)
    %   dt - 时间步长(分钟)
    %
    % 输出：
    %   EVs - 更新后的电动汽车结构体数组，包含分配的功率
    
    %% 1. 筛选当前时刻参与聚合的EV
    % 条件：在充电时间范围内且状态为ON或LockON
    active_EV_indices = [];
    total_P_N = 0;
    
    for i = 1:length(EVs)
        ev = EVs(i);
        if (strcmp(ev.state, 'ON') || strcmp(ev.state, 'LockON'))
            active_EV_indices = [active_EV_indices, i];
            total_P_N = total_P_N + ev.P_N;
        end
    end
    
    %% 2. 如果没有活跃EV，将所有功率设为0
    if isempty(active_EV_indices)
        for i = 1:length(EVs)
            EVs(i).P_base = 0;
            EVs(i).BaseSchedule = zeros(size(P_base_opt));
        end
        return;
    end
    
    %% 3. 计算每个活跃EV的权重(基于额定功率)
    weights = zeros(length(active_EV_indices), 1);
    for i = 1:length(active_EV_indices)
        idx = active_EV_indices(i);
        weights(i) = EVs(idx).P_N / total_P_N;
    end
    
    %% 4. 分配功率给每个活跃EV
    for i = 1:length(active_EV_indices)
        idx = active_EV_indices(i);
        
        % 计算分配的功率序列
        allocated_P = P_base_opt * weights(i);
        
        % 确保不超过EV的额定功率
        allocated_P = min(allocated_P, EVs(idx).P_N);
        
        % 确保功率非负
        allocated_P = max(allocated_P, 0);
        
        % 更新EV的功率信息
        EVs(idx).P_base = allocated_P(1);  % 当前时刻功率
        EVs(idx).BaseSchedule = allocated_P; % 完整功率序列
    end
    
    %% 5. 将非活跃EV的功率设为0
    for i = 1:length(EVs)
        if ~ismember(i, active_EV_indices)
            EVs(i).P_base = 0;
            EVs(i).BaseSchedule = zeros(size(P_base_opt));
        end
    end
end
