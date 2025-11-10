function SimData = load_simulation_data_from_chunks(chunk_dir, chunk_load_config, direction)
    % load_simulation_data_from_chunks (v9 - 支持按行数限制加载个体数据)
    
    if isempty(chunk_load_config)
        error('未指定要加载的配置 (chunk_load_config 为空)。');
    end
    
    all_p_AC = {};
    all_p_EV = {};
    
    % [!!! 移除：不再累加文件中预先聚合的数据 !!!]
    % P_AC_agg_accumulator = [];
    % P_EV_agg_accumulator = [];

    T_ac = 0;
    T_ev = 0;
    
    % 根据调节方向确定要读取的字段名
    if strcmpi(direction, 'Up')
        p_field_ac_ind = 'AC_Up_Individual';
        p_field_ev_ind = 'EV_Up_Individual';
    elseif strcmpi(direction, 'Down')
        p_field_ac_ind = 'AC_Down_Individual';
        p_field_ev_ind = 'EV_Down_Individual';
    else
        error('无效的调节方向: %s. 必须是 "Up" 或 "Down".', direction);
    end

    % [!!! 修改点：遍历 chunk_load_config 结构体 !!!]
    for i = 1:length(chunk_load_config)
        config = chunk_load_config(i); % 获取当前配置
        file_name = config.file;
        file_path = fullfile(chunk_dir, file_name);
        
        if ~exist(file_path, 'file')
            warning('文件未找到: %s，已跳过。', file_path);
            continue;
        end
        
        fprintf('  正在加载文件: %s ...\n', file_name);
        try
            data = load(file_path);
            if ~isfield(data, 'results')
                warning('文件 %s 中无 results 字段，跳过。', file_name);
                continue;
            end
            res = data.results;

            % --- 1. 加载个体数据 (!!! 应用行数限制 !!!) ---
            if isfield(res, p_field_ac_ind) && ~isempty(res.(p_field_ac_ind))
                ac_data_chunk = res.(p_field_ac_ind);
                % 确定要加载的行数
                rows_to_take_ac = min(config.ac_rows, size(ac_data_chunk, 1));
                % 提取指定行数的数据
                all_p_AC{end+1} = ac_data_chunk(1:rows_to_take_ac, :);
                
                if T_ac == 0, T_ac = size(all_p_AC{end}, 2); end
                fprintf('    -> %s: 加载 %d / %d 行 (限制: %s)\n', ...
                    p_field_ac_ind, rows_to_take_ac, size(ac_data_chunk, 1), num2str(config.ac_rows));
            else
                 fprintf('    -> %s: 未找到或为空。\n', p_field_ac_ind);
            end
            
            if isfield(res, p_field_ev_ind) && ~isempty(res.(p_field_ev_ind))
                ev_data_chunk = res.(p_field_ev_ind);
                % 确定要加载的行数
                rows_to_take_ev = min(config.ev_rows, size(ev_data_chunk, 1));
                % 提取指定行数的数据
                all_p_EV{end+1} = ev_data_chunk(1:rows_to_take_ev, :);
                
                if T_ev == 0, T_ev = size(all_p_EV{end}, 2); end
                 fprintf('    -> %s: 加载 %d / %d 行 (限制: %s)\n', ...
                    p_field_ev_ind, rows_to_take_ev, size(ev_data_chunk, 1), num2str(config.ev_rows));
            else
                 fprintf('    -> %s: 未找到或为空。\n', p_field_ev_ind);
            end

            % --- [!!! 移除：不再加载聚合数据 !!!] ---
            % (原始加载 P_AC_agg_accumulator 和 P_EV_agg_accumulator 的代码已删除)

        catch ME
            warning('加载文件 %s 时出错: %s', file_name, ME.message);
        end
    end
    
    % 数据验证
    if T_ac == 0 && T_ev == 0 && (isempty(all_p_AC) || isempty(all_p_EV))
        warning('未能加载到任何有效的 AC 或 EV 数据。');
        % 即使为空也继续，以便脚本可以处理 0 设备的情况
    end
    
    if T_ac ~= T_ev && T_ac > 0 && T_ev > 0
        warning('AC 和 EV 数据的时间步长不一致 (%d vs %d)。', T_ac, T_ev);
    end
    
    % 确定最终 T (优先使用AC，否则使用EV)
    if T_ac > 0
        SimData.T = T_ac;
    else
        SimData.T = T_ev;
    end
    
    % 合并个体数据
    if ~isempty(all_p_AC)
        SimData.p_AC = cat(1, all_p_AC{:});
    else
        SimData.p_AC = zeros(0, SimData.T); % 创建一个 0 x T 的空矩阵
    end
    
    if ~isempty(all_p_EV)
        SimData.p_EV = cat(1, all_p_EV{:});
    else
        SimData.p_EV = zeros(0, SimData.T); % 创建一个 0 x T 的空矩阵
    end

    SimData.nAC = size(SimData.p_AC, 1);
    SimData.nEV = size(SimData.p_EV, 1);
    SimData.N_total = SimData.nAC + SimData.nEV;
    
    
    % --- [!!! 修改点：从加载的个体数据重新计算聚合值 !!!] ---
    fprintf('  正在根据加载的 %d 台 AC 和 %d 台 EV 重新计算聚合潜力...\n', SimData.nAC, SimData.nEV);
    
    if SimData.nAC > 0
        SimData.P_AC_agg_loaded = sum(SimData.p_AC, 1);
    else
        SimData.P_AC_agg_loaded = zeros(1, SimData.T);
    end
    
    if SimData.nEV > 0
        SimData.P_EV_agg_loaded = sum(SimData.p_EV, 1);
    else
        SimData.P_EV_agg_loaded = zeros(1, SimData.T);
    end
    
end