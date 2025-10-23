%% generateEVData_48.m
function generateEV_48(num_devices_ev, area_type, selected_models)
    % 参数验证
    if nargin < 3
        selected_models = 'all';  % 默认选择全部车型
    end
    if nargin < 2
        area_type = 'residential';  % 默认生成居民区数据
    end

    %% ===== 车型数据库 =====
    evModels = [
        struct('name','比亚迪秦PLUS',  'p_slow',6.6, 'p_fast',60,  'C_EV',[18.32,26.86], 'prob',0.35),...
        struct('name','比亚迪海鸥',    'p_slow',6.6, 'p_fast',40,  'C_EV',[30.08,38.88], 'prob',0.25),...
        struct('name','特斯拉Model Y','p_slow',11,  'p_fast',250, 'C_EV',[60,78.4],    'prob',0.15),...
        struct('name','比亚迪元PLUS', 'p_slow',7,   'p_fast',90,  'C_EV',[50.12,60.48],'prob',0.15),...
        struct('name','比亚迪宋PLUS DM-i', 'p_slow',6.6, 'p_fast',70,  'C_EV',[18.32,26.86],'prob',0.1),...
    ];

    % 车型筛选
    if ~strcmpi(selected_models, 'all')
        valid_models_idx = ismember({evModels.name}, selected_models);
        if ~any(valid_models_idx)
            error('未找到指定车型，可用车型: %s', strjoin({evModels.name}, ', '));
        end
        evModels = evModels(valid_models_idx);
        if ~isempty(evModels)
            total_prob = sum([evModels.prob]);
            if total_prob > 0
                probs_vector = [evModels.prob] / total_prob;
                for k_model = 1:length(evModels)
                    evModels(k_model).prob = probs_vector(k_model);
                end
            elseif ~isempty(evModels) % 如果总概率为0但仍有模型，则平均分配
                default_prob = 1/length(evModels);
                for k_model = 1:length(evModels)
                    evModels(k_model).prob = default_prob;
                end
            end
        else
             error('筛选后没有可用的EV模型。');
        end
    end

    %% ===== 区域参数配置 =====
    timeSlots = {}; 
    area_name_str = ''; 
    charge_type_handle = [];
    efficiency_range = [];
    mean_arrival_h_config = 0; 
    std_dev_h_config = 0;
    % min_arrival_h_boundary 和 max_arrival_h_boundary 不再用于严格的重采样截断，
    % 但可以保留用于指导参数选择或进行非常宽松的最终检查（如果需要）。
    % 当前逻辑将主要依赖正态分布和mod运算。

    switch lower(area_type)
        case 'residential'
            timeSlots = { 
                struct('range',0:4,  'dist',[70,20,10], 'max_dur',[24,8,3]),    
                struct('range',5:8,  'dist',[20,70,10], 'max_dur',[19,8,3]),   
                struct('range',9:15, 'dist',[10,20,70], 'max_dur',[15,8,3]),  
                struct('range',16:17,'dist',[20,70,10], 'max_dur',[8,8,3]),    
                struct('range',18:23,'dist',[70,20,10], 'max_dur',[24,8,3])    
            };
            area_name_str = '居民区';
            charge_type_handle = @(m) m.p_slow;  
            efficiency_range = [0.85, 0.95];
            % --- 居民区并网时间正态分布参数 ---
            mean_arrival_h_config = 19.0;  % 平均到达时间 19:00 (傍晚)
            std_dev_h_config = 2.0;        % 标准差2.0小时，以允许一定的分散性
            
        case 'workplace'
            timeSlots = { 
                struct('range',0:4,  'dist',[10,20,70], 'max_dur',[3,8,3]),    
                struct('range',5:8,  'dist',[20,10,70], 'max_dur',[3,8,3]),    
                struct('range',9:15, 'dist',[20,70,10], 'max_dur',[6,8,3]),   
                struct('range',16:17,'dist',[10,20,70], 'max_dur',[3,8,3]),   
                struct('range',18:23,'dist',[10,20,70], 'max_dur',[3,8,3])    
            };
            area_name_str = '工作区';
            charge_type_handle = @(m) m.p_fast;  
            efficiency_range = [0.90, 0.95];
            % --- 工作区并网时间正态分布参数 ---
            mean_arrival_h_config = 8.5;   % 平均到达时间 08:30 (早上7-9点之间)
            std_dev_h_config = 1.0;        % 标准差1.0小时，确保主要在7-9点但允许少量超出
            
        otherwise
            error('无效区域类型，可选: residential, workplace');
    end

    %% ===== 数据生成 =====
    [C_EV,eta,p_on,t_in,t_dep,E_in,E_tar] = deal(zeros(num_devices_ev,1));
    areaType_cell_output = repmat({area_name_str}, num_devices_ev, 1); 

    for i = 1:num_devices_ev
        cumProb = cumsum([evModels.prob]);
        selectedModel = evModels(find(rand() <= cumProb,1,'first'));

        C_EV(i) = selectedModel.C_EV(randi(length(selectedModel.C_EV)));
        p_on(i) = charge_type_handle(selectedModel);
        eta(i) = efficiency_range(1) + diff(efficiency_range)*rand();

        valid = false;
        while ~valid
            % --- 修改：根据区域类型生成并网时间 (正态分布，允许尾部环绕) ---
            generated_hour_float = mean_arrival_h_config + std_dev_h_config * randn();
            
            % 将生成的小时通过取模操作映射到 [0, 24) 区间
            % MATLAB的mod(a,m)当m为正时，结果与m同号或为0。
            % 例如 mod(25, 24) = 1; mod(-1, 24) = 23. 结果总是在 [0, 24)
            t_in_h_float_daily = mod(generated_hour_float, 24);
            
            t_in_h_for_slot_lookup = floor(t_in_h_float_daily); 
            % 确保 t_in_h_for_slot_lookup 在 [0,23] 范围内，尽管 mod(X,24) 后 floor 应该已经是了
            t_in_h_for_slot_lookup = max(0, min(23, t_in_h_for_slot_lookup));

            t_in_m = randi([0,59]); 
            t_in(i) = t_in_h_float_daily + t_in_m/60;
            % t_in(i) 是车辆在其典型一天中的接入时间（0-23.99小时制）
            % 后续 t_dep(i) = t_in(i) + duration; 可以使 t_dep(i) > 24
            % --- 结束并网时间生成修改 ---
            
            currentSlot = timeSlots{find(cellfun(@(x) ismember(t_in_h_for_slot_lookup, x.range), timeSlots),1)};
            
            chargeType_idx = 0; 
            minDur = 0; maxDur_from_slot = 0; 
            randVal = rand()*100;
            if randVal <= currentSlot.dist(1)
                chargeType_idx = 1; 
                minDur = (strcmp(area_name_str,'居民区')*7) + (strcmp(area_name_str,'工作区')*2); 
                maxDur_from_slot = currentSlot.max_dur(1);
            elseif randVal <= sum(currentSlot.dist(1:2))
                chargeType_idx = 2; 
                minDur = (strcmp(area_name_str,'居民区')*3) + (strcmp(area_name_str,'工作区')*1); 
                maxDur_from_slot = currentSlot.max_dur(2);
            else
                chargeType_idx = 3; 
                minDur = 0;
                maxDur_from_slot = currentSlot.max_dur(3);
            end

            effectiveMaxDur = maxDur_from_slot; 
            if strcmp(area_name_str,'工作区')
                time_until_workday_end = 0;
                workday_end_hour = 18.0; 
                if t_in(i) < workday_end_hour 
                    time_until_workday_end = workday_end_hour - t_in(i);
                else 
                    time_until_workday_end = 0.25; 
                end
                effectiveMaxDur = min(effectiveMaxDur, time_until_workday_end);
            end
            
            if effectiveMaxDur < minDur
                 if strcmp(area_name_str,'工作区') && chargeType_idx == 3 && t_in(i) >= (18.0 - effectiveMaxDur - 0.1) 
                     duration = max(0, effectiveMaxDur); 
                 else
                    continue; 
                 end
            end
            
            if effectiveMaxDur <= 1e-3 && minDur <= 1e-3 
                duration = 0;
            elseif effectiveMaxDur >= minDur 
                 duration = minDur + (effectiveMaxDur - minDur)*rand();
            else 
                 duration = 0; 
            end
            
            t_dep(i) = t_in(i) + duration;

            if strcmp(area_name_str, '居民区')
                E_in(i) = (0.1 + 0.2*rand()) * C_EV(i); 
            else 
                E_in(i) = (0.3 + 0.4*rand()) * C_EV(i); 
            end

            max_charge_possible = p_on(i) * eta(i) * duration;
            
            target_soc_upper_limit = 0.95;
            if strcmp(area_name_str, '工作区')
                target_soc_upper_limit = 0.90; 
            end
            E_tar_max_limit = min(E_in(i) + max_charge_possible, target_soc_upper_limit * C_EV(i));
            
            if E_tar_max_limit <= E_in(i) + 1e-3 
                 if duration > 1e-3 
                    continue; 
                 else 
                    E_tar(i) = E_in(i);
                 end
            else
                 E_tar(i) = E_in(i) + (E_tar_max_limit - E_in(i)) * (0.2 + 0.8*rand()^2); 
            end

            % --- 修改：增加离网时间约束 ---
            required_power_check = 0;
            if duration > 1e-6 
                required_power_check = (E_tar(i) - E_in(i)) / (eta(i) * duration);
            end
            
            if t_dep(i) < 32 && (required_power_check <= p_on(i) + 1e-3) && (E_tar(i) >= E_in(i) - 1e-3)
                 if duration < 1e-6 && abs(E_tar(i) - E_in(i)) < 1e-3
                     valid = true;
                 elseif duration >= 1e-6 && E_tar(i) > E_in(i) - 1e-3
                     valid = true;
                 elseif duration >= 1e-6 && abs(E_tar(i) - E_in(i)) < 1e-3 && p_on(i) == 0
                     valid = true; 
                 end
            end
            % --- 结束修改 ---
        end

        C_EV(i) = round(C_EV(i), 2);
        eta(i) = round(eta(i), 2);
        t_in(i) = round(t_in(i), 2);
        t_dep(i) = round(t_dep(i), 2);
        E_in(i) = round(E_in(i), 2); 
        E_tar(i) = round(E_tar(i), 2);
    end

    %% ===== 数据保存 =====
    evData = struct(...
        'ID', num2cell((1:num_devices_ev)'),...
        'C_EV', num2cell(C_EV),...
        'eta', num2cell(eta),...
        'p_on', num2cell(p_on),...
        'E_in', num2cell(E_in),...
        'E_tar', num2cell(E_tar),...
        't_in', num2cell(t_in),...
        't_dep', num2cell(t_dep),...
        'AreaType', areaType_cell_output,... 
        'r', num2cell(repmat(0.025, num_devices_ev,1)),...
        'SOC', 0,... 
        'p_incentive', num2cell(round(60*rand(num_devices_ev,1), 1))...
    );

    currentDir = fileparts(mfilename('fullpath'));
    projectRoot = fullfile(currentDir, '..'); 
    saveDir = fullfile(projectRoot, '0input_data');
    if ~exist(saveDir, 'dir'), mkdir(saveDir); end
    
    filename_to_save = sprintf('2EV_%s.xlsx', area_type); 
    writetable(struct2table(evData), fullfile(saveDir, filename_to_save));
    disp(['EV数据文件已生成：' filename_to_save]);
end