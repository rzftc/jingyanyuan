% load_simulation_data_from_chunks.m
function SimData = load_simulation_data_from_chunks(chunk_dir, selected_mat_files, direction)
    % 从 chunk_results 目录加载并聚合 *指定* 的分块数据
    
    if isempty(selected_mat_files)
        error('您没有在 main_vpp_optimizer.m 的 "selected_mat_files" 变量中指定任何文件。');
    end
    
    % 初始化空的 cell 数组用于聚合
    all_p_AC = {};
    all_p_EV = {};
    
    T = 0; % 时间步长
    
    % *** 修改：遍历指定的文件列表，而不是扫描目录 ***
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
            
            % 确定调节方向
            if strcmpi(direction, 'Up')
                p_field_ac = 'AC_Up_Individual';
                p_field_ev = 'EV_Up_Individual';
            elseif strcmpi(direction, 'Down')
                % 注意：下调潜力是负值，优化时通常使用其绝对值
                % 为保持一致性，我们取绝对值，假设成本 c_i 也是正的
                p_field_ac = 'AC_Down_Individual';
                p_field_ev = 'EV_Down_Individual';
            else
                error('无效的调节方向: %s. 必须是 "Up" 或 "Down".', direction);
            end

            if isfield(data, 'results') && isfield(data.results, p_field_ac)
                ac_data = data.results.(p_field_ac);
                if strcmpi(direction, 'Down')
                    ac_data = abs(ac_data); % 取绝对值
                end
                all_p_AC{end+1} = ac_data;
                if T == 0, T = size(ac_data, 2); end
            end
            
            if isfield(data, 'results') && isfield(data.results, p_field_ev)
                ev_data = data.results.(p_field_ev);
                 if strcmpi(direction, 'Down')
                    ev_data = abs(ev_data); % 取绝对值
                end
                all_p_EV{end+1} = ev_data;
                if T == 0, T = size(ev_data, 2); end
            end
            
        catch ME
            warning('加载文件 %s 时出错: %s. 跳过此文件。', file_name, ME.message);
        end
    end
    
    % 沿设备维度 (维度1) 合并所有分块
    SimData.p_AC = cat(1, all_p_AC{:}); % [nAC x T]
    SimData.p_EV = cat(1, all_p_EV{:}); % [nEV x T]
    
    SimData.nAC = size(SimData.p_AC, 1);
    SimData.nEV = size(SimData.p_EV, 1);
    SimData.N_total = SimData.nAC + SimData.nEV;
    SimData.T = T;
    
    if SimData.N_total == 0
        error('未能从指定的 .mat 文件中加载任何 AC 或 EV 数据。');
    end
end