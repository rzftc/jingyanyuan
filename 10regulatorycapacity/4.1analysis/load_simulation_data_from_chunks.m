% load_simulation_data_from_chunks.m
function SimData = load_simulation_data_from_chunks(chunk_dir, selected_mat_files, direction)
    % (v3) 从 chunk_results 加载数据，并验证 AC 和 EV 数据是否 *同时* 存在
    
    if isempty(selected_mat_files)
        error('您没有在主脚本的 "selected_mat_files" 变量中指定任何文件。');
    end
    
    all_p_AC = {};
    all_p_EV = {};
    T_ac = 0;
    T_ev = 0;
    
    % 确定字段名
    if strcmpi(direction, 'Up')
        p_field_ac = 'AC_Up_Individual';
        p_field_ev = 'EV_Up_Individual';
    elseif strcmpi(direction, 'Down')
        p_field_ac = 'AC_Down_Individual';
        p_field_ev = 'EV_Down_Individual';
    else
        error('无效的调节方向: %s. 必须是 "Up" 或 "Down".', direction);
    end

    % 遍历指定的文件列表
    for i = 1:length(selected_mat_files)
        file_name = selected_mat_files{i};
        file_path = fullfile(chunk_dir, file_name);
        
        if ~exist(file_path, 'file')
            warning('指定的文件未找到: %s. 跳过此文件。', file_path);
            continue;
        end
        
        fprintf('  正在加载指定文件: %s\n', file_name);
        try
            data = load(file_path);
            
            % 检查 AC 数据
            if isfield(data, 'results') && isfield(data.results, p_field_ac)
                ac_data = data.results.(p_field_ac);
                if strcmpi(direction, 'Down'), ac_data = abs(ac_data); end
                all_p_AC{end+1} = ac_data;
                if T_ac == 0, T_ac = size(ac_data, 2); end
            end
            
            % 检查 EV 数据
            if isfield(data, 'results') && isfield(data.results, p_field_ev)
                ev_data = data.results.(p_field_ev);
                 if strcmpi(direction, 'Down'), ev_data = abs(ev_data); end
                all_p_EV{end+1} = ev_data;
                if T_ev == 0, T_ev = size(ev_data, 2); end
            end
            
        catch ME
            warning('加载文件 %s 时出错: %s. 跳过此文件。', file_name, ME.message);
        end
    end
    
    % --- 数据验证 ---
    if isempty(all_p_AC) && isempty(all_p_EV)
        error('未能从指定的 .mat 文件中加载任何 AC 或 EV 数据。');
    end
    
    % *** 关键检查：确保两种数据都存在 ***
    if isempty(all_p_AC)
        error(['加载失败: 未能在 .mat 文件中找到 AC 数据 (字段: %s)。\n', ...
               '请确保 ac_ev_simulation_block.m 在 runAC=true 的设置下运行。'], p_field_ac);
    end
    if isempty(all_p_EV)
        error(['加载失败: 未能在 .mat 文件中找到 EV 数据 (字段: %s)。\n', ...
               '请确保 ac_ev_simulation_block.m 在 runEV=true 的设置下运行。'], p_field_ev);
    end

    if T_ac ~= T_ev
        warning('AC 和 EV 数据的时间步长不匹配 (%d vs %d)，可能导致错误。', T_ac, T_ev);
    end
    
    % 沿设备维度 (维度1) 合并所有分块
    SimData.p_AC = cat(1, all_p_AC{:}); % [nAC x T]
    SimData.p_EV = cat(1, all_p_EV{:}); % [nEV x T]
    
    SimData.nAC = size(SimData.p_AC, 1);
    SimData.nEV = size(SimData.p_EV, 1);
    SimData.N_total = SimData.nAC + SimData.nEV;
    SimData.T = T_ac; % 假设 AC 的时间步长是权威的
end