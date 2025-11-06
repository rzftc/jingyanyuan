function SimData = load_simulation_data_from_chunks(chunk_dir, selected_mat_files, direction)
    % load_simulation_data_from_chunks (v8 - 直接读取聚合数据 & 保留原始符号)
    
    if isempty(selected_mat_files)
        error('未指定要加载的文件 (selected_mat_files 为空)。');
    end
    
    all_p_AC = {};
    all_p_EV = {};
    
    % [新增] 初始化聚合功率累加器
    P_AC_agg_accumulator = [];
    P_EV_agg_accumulator = [];

    T_ac = 0;
    T_ev = 0;
    
    % 根据调节方向确定要读取的字段名
    if strcmpi(direction, 'Up')
        p_field_ac_ind = 'AC_Up_Individual';
        p_field_ev_ind = 'EV_Up_Individual';
        p_field_ac_agg = 'AC_Up';    % 聚合数据字段
        p_field_ev_agg = 'EV_Up';
    elseif strcmpi(direction, 'Down')
        p_field_ac_ind = 'AC_Down_Individual';
        p_field_ev_ind = 'EV_Down_Individual';
        p_field_ac_agg = 'AC_Down';  % 聚合数据字段
        p_field_ev_agg = 'EV_Down';
    else
        error('无效的调节方向: %s. 必须是 "Up" 或 "Down".', direction);
    end

    for i = 1:length(selected_mat_files)
        file_name = selected_mat_files{i};
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

            % --- 1. 加载个体数据 (保留原始符号) ---
            if isfield(res, p_field_ac_ind)
                all_p_AC{end+1} = res.(p_field_ac_ind);
                if T_ac == 0, T_ac = size(all_p_AC{end}, 2); end
            end
            if isfield(res, p_field_ev_ind)
                all_p_EV{end+1} = res.(p_field_ev_ind);
                if T_ev == 0, T_ev = size(all_p_EV{end}, 2); end
            end

            % --- [新增] 2. 直接读取并累加聚合数据 ---
            % 确保数据为行向量 (1 x T) 以便统一累加
            if isfield(res, p_field_ac_agg)
                agg_data = res.(p_field_ac_agg);
                if iscolumn(agg_data), agg_data = agg_data'; end 
                if isempty(P_AC_agg_accumulator)
                    P_AC_agg_accumulator = agg_data;
                else
                    P_AC_agg_accumulator = P_AC_agg_accumulator + agg_data;
                end
            end
             if isfield(res, p_field_ev_agg)
                agg_data = res.(p_field_ev_agg);
                if iscolumn(agg_data), agg_data = agg_data'; end
                if isempty(P_EV_agg_accumulator)
                    P_EV_agg_accumulator = agg_data;
                else
                    P_EV_agg_accumulator = P_EV_agg_accumulator + agg_data;
                end
            end

        catch ME
            warning('加载文件 %s 时出错: %s', file_name, ME.message);
        end
    end
    
    % 数据验证
    if isempty(all_p_AC) || isempty(all_p_EV)
        error('未能加载到有效的 AC 或 EV 数据。');
    end
    if T_ac ~= T_ev
        warning('AC 和 EV 数据的时间步长不一致 (%d vs %d)。', T_ac, T_ev);
    end
    
    % 合并个体数据
    SimData.p_AC = cat(1, all_p_AC{:});
    SimData.p_EV = cat(1, all_p_EV{:});
    SimData.nAC = size(SimData.p_AC, 1);
    SimData.nEV = size(SimData.p_EV, 1);
    SimData.N_total = SimData.nAC + SimData.nEV;
    SimData.T = T_ac;
    
    % --- [新增] 存储直接读取的聚合数据 ---
    % 如果文件中没有聚合数据，作为回退，使用现场求和
    if isempty(P_AC_agg_accumulator)
        SimData.P_AC_agg_loaded = sum(SimData.p_AC, 1);
    else
        SimData.P_AC_agg_loaded = P_AC_agg_accumulator;
    end
    if isempty(P_EV_agg_accumulator)
        SimData.P_EV_agg_loaded = sum(SimData.p_EV, 1);
    else
        SimData.P_EV_agg_loaded = P_EV_agg_accumulator;
    end
end