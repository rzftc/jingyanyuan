%% generateEVData.m
function generateEVData(num_devices_ev, area_type, selected_models)
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
        valid_models = ismember({evModels.name}, selected_models);
        if ~any(valid_models)
            error('未找到指定车型，可用车型: %s', strjoin({evModels.name}, ', '));
        end
        evModels = evModels(valid_models);
        % 重新归一化概率
        total_prob = sum([evModels.prob]);
        [evModels.prob] = arrayfun(@(x) x/total_prob, [evModels.prob]); % 修正概率分配方式
    end

    %% ===== 区域参数配置 =====
    timeSlots = {}; % 初始化 timeSlots
    area_name = '';
    charge_type_handle = [];
    efficiency_range = [];
    mean_arrival_hour = 0;
    std_dev_arrival_hour = 0;
    min_arrival_hour = 0; % 新增：允许的最早到达小时
    max_arrival_hour = 0; % 新增：允许的最晚到达小时


    switch lower(area_type)
        case 'residential'
            timeSlots = { % 居民区充电行为分布
                struct('range',0:4,  'dist',[70,20,10], 'max_dur',[24,8,3]),    
                struct('range',5:8,  'dist',[20,70,10], 'max_dur',[19,8,3]),   
                struct('range',9:15, 'dist',[10,20,70], 'max_dur',[15,8,3]),  
                struct('range',16:17,'dist',[20,70,10], 'max_dur',[8,8,3]),    
                struct('range',18:23,'dist',[70,20,10], 'max_dur',[6,8,3])     
            };
            area_name = '居民区';
            charge_type_handle = @(m) m.p_slow;  % 居民区使用慢充
            efficiency_range = [0.85, 0.95];
            % --- 修改：居民区并网时间参数 ---
            mean_arrival_hour = 19.0;  % 傍晚平均到达时间 (例如 19:00)
            std_dev_arrival_hour = 1.5;   % 标准差 (例如 1.5 小时)
            min_arrival_hour = 16.0;   % 例如，允许最早16:00
            max_arrival_hour = 23.99;  % 例如，允许最晚接近午夜
            % --- 结束修改 ---
            
        case 'workplace'
            timeSlots = { % 工作区充电行为分布
                struct('range',0:4,  'dist',[10,20,70], 'max_dur',[3,8,3]),    
                struct('range',5:8,  'dist',[20,10,70], 'max_dur',[3,8,3]),    
                struct('range',9:15, 'dist',[20,70,10], 'max_dur',[6,8,3]),   
                struct('range',16:17,'dist',[10,20,70], 'max_dur',[3,8,3]),   
                struct('range',18:23,'dist',[10,20,70], 'max_dur',[3,8,3])    
            };
            area_name = '工作区';
            charge_type_handle = @(m) m.p_fast;  % 工作区使用快充
            efficiency_range = [0.90, 0.95];
            % --- 修改：工作区并网时间参数 ---
            mean_arrival_hour = 8.0;   % 早上平均到达时间 (例如 08:00)
            std_dev_arrival_hour = 0.75;  % 标准差 (例如 0.75 小时)
            min_arrival_hour = 7.0;    % 例如，允许最早07:00
            max_arrival_hour = 9.99;   % 例如，允许最晚接近10:00 (确保集中在7-9点)
            % --- 结束修改 ---
            
        otherwise
            error('无效区域类型，可选: residential, workplace');
    end

    %% ===== 数据生成 =====
    [C_EV,eta,p_on,t_in,t_dep,E_in,E_tar] = deal(zeros(num_devices_ev,1));
    areaType_cell = repmat({area_name}, num_devices_ev, 1); % 使用新变量名以避免与area_type函数输入冲突

    for i = 1:num_devices_ev
        % 车型选择
        cumProb = cumsum([evModels.prob]);
        selectedModel = evModels(find(rand() <= cumProb,1));

        % 充电参数设置
        C_EV(i) = selectedModel.C_EV(randi(length(selectedModel.C_EV)));
        p_on(i) = charge_type_handle(selectedModel); % 使用句柄
        eta(i) = efficiency_range(1) + diff(efficiency_range)*rand();

        % 时间参数生成
        valid = false;
        t_in_h_integer = 0; % 初始化用于 timeslots 的整数小时

        while ~valid
            % --- 修改：根据区域类型生成并网时间 ---
            generated_hour_float = mean_arrival_hour + std_dev_arrival_hour * randn();
            
            % 将生成的小时限制在定义的范围内并进行循环处理或截断
            % 采用截断法，如果超出范围则重新生成，以更好地保证分布集中性
            if generated_hour_float < min_arrival_hour || generated_hour_float > max_arrival_hour
                continue; % 如果超出期望范围，则重新生成，避免mod导致的不自然分布
            end
            
            t_in_h_float = generated_hour_float;
            t_in_h_integer = floor(t_in_h_float); % 用于查找 timeSlots 的整数小时
                                                % 确保 t_in_h_integer 在 [0,23] 范围内
            t_in_h_integer = mod(t_in_h_integer, 24); 


            t_in_m = randi([0,59]); % 分钟数仍然随机
            t_in(i) = t_in_h_float + t_in_m/60;
            
            % 确保 t_in 在 [0, 24) 范围内 (如果上面的截断逻辑不够严格)
            % 或者如果允许跨天生成，这里的逻辑需要调整。当前假设在一天内。
            t_in(i) = mod(t_in(i), 24); 
            if t_in(i) < 0 
                t_in(i) = t_in(i) + 24;
            end
            % 更新 t_in_h_integer 以确保与最终的 t_in(i) 对应的小时部分一致
            t_in_h_integer = floor(t_in(i));
            % --- 结束修改 ---
            
            currentSlot = timeSlots{find(cellfun(@(x) ismember(t_in_h_integer, x.range), timeSlots),1)};
            
            % 充电类型选择
            randVal = rand()*100;
            chargeType = 0; % 初始化
            minDur = 0; maxDur = 0; % 初始化
            if randVal <= currentSlot.dist(1)
                chargeType = 1; 
                minDur = strcmp(area_name,'居民区')*7 + 2; 
                maxDur = currentSlot.max_dur(1);
            elseif randVal <= sum(currentSlot.dist(1:2))
                chargeType = 2; 
                minDur = strcmp(area_name,'居民区')*3 + 1;
                maxDur = currentSlot.max_dur(2);
            else
                chargeType = 3; 
                minDur = 0;
                maxDur = currentSlot.max_dur(3);
            end

            % 持续时间计算
            maxPossibleDuration = 0;
            if strcmp(area_name,'工作区')
                 % 工作区车辆通常在当天工作时间结束前离开
                 % 假设工作区最晚充电到18:00，如果t_in已经很晚，则可用时长很短
                if t_in(i) < 18.0 
                    maxPossibleDuration = 18.0 - t_in(i);
                else % 如果接入时间晚于或等于18点，则认为当天不再充电或充电时间极短
                    maxPossibleDuration = 0.5; % 给予一个很小的窗口，或者直接continue
                end
            else % 居民区
                maxPossibleDuration = 24.0 - t_in(i); % 当天剩余时间
            end
            
            % 有效最大持续时间不能超过 time slot 定义的 maxDur 和当天剩余可用时间
            effectiveMaxDur = min(maxDur, maxPossibleDuration);

            if effectiveMaxDur < minDur 
                if strcmp(area_name,'工作区') && t_in(i) >= 17.5 && chargeType == 3 % 工作区短时充电，在17:30后接入可以允许短时
                     duration = effectiveMaxDur; % 使用可能的短时
                else
                    continue; % 如果有效最大时长小于最小时长，重新生成接入时间等
                end
            end
            
            if effectiveMaxDur <= 0 && minDur <=0 % 如果只能进行0时长充电
                duration = 0;
            elseif effectiveMaxDur > minDur
                 duration = minDur + (effectiveMaxDur - minDur)*rand();
            else % effectiveMaxDur == minDur
                 duration = minDur;
            end

            t_dep(i) = t_in(i) + duration;

            % 电量参数生成
            % 根据区域类型调整初始电量范围
            if strcmp(area_name, '居民区')
                E_in(i) = (0.1 + 0.2*rand()) * C_EV(i); % 居民区初始电量较低: 10%-30%
            else % 工作区
                E_in(i) = (0.3 + 0.4*rand()) * C_EV(i); % 工作区初始电量可能稍高: 30%-70%
            end

            max_charge_possible = p_on(i) * eta(i) * duration;
            
            % 目标充电量上限调整，工作区可能不会充满到95%
            target_soc_upper_limit = 0.95;
            if strcmp(area_name, '工作区')
                target_soc_upper_limit = 0.9; % 工作区目标SOC上限稍低
            end
            E_tar_max_limit = min(E_in(i) + max_charge_possible, target_soc_upper_limit * C_EV(i));
            
            if E_tar_max_limit <= E_in(i) % 如果无法充入更多电量
                 if duration > 1e-3 % 只有在有充电时长时才认为无效
                    continue; % 重新生成
                 else % 0时长充电，目标等于初始
                    E_tar(i) = E_in(i);
                 end
            else
                 E_tar(i) = E_in(i) + (E_tar_max_limit - E_in(i)) * (0.2 + 0.8*rand()^2); % 非线性增加目标电量
            end


            % 验证条件
            required_power_check = 0;
            if duration > 1e-6 % 避免除以零
                required_power_check = (E_tar(i) - E_in(i)) / (eta(i) * duration);
            end
            
            % 确保离开时间不超过24点（对于单日模型）或满足特定逻辑
            % 并且实际需求功率小于等于额定功率，且目标电量大于初始电量（除非duration为0）
            if t_dep(i) <= 24 && (required_power_check <= p_on(i) + 1e-3) && (E_tar(i) >= E_in(i) - 1e-3)
                 if duration < 1e-6 && abs(E_tar(i) - E_in(i)) < 1e-3 % 0时长，目标应等于初始
                     valid = true;
                 elseif duration >= 1e-6 && E_tar(i) > E_in(i) - 1e-3 % 有时长，目标应不小于初始
                     valid = true;
                 elseif duration >= 1e-6 && abs(E_tar(i) - E_in(i)) < 1e-3 && p_on(i) == 0 % 功率为0，无法充电
                     valid = true; 
                 end
            end
        end

        % 数据精度处理
        C_EV(i) = round(C_EV(i), 2);
        eta(i) = round(eta(i), 2);
        t_in(i) = round(t_in(i), 2);
        t_dep(i) = round(t_dep(i), 2);
        E_in(i) = round(E_in(i), 2); % E_in, E_tar 精度调整
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
        'AreaType', areaType_cell,... % 使用新的变量名
        'r', num2cell(repmat(0.025, num_devices_ev,1)),...
        'SOC', 0,... % SOC通常在主程序中根据E_in, C_EV等计算或更新
        'p_incentive', num2cell(round(60*rand(num_devices_ev,1), 1))...
    );

    currentDir = fileparts(mfilename('fullpath'));
    projectRoot = fullfile(currentDir, '..'); 
    saveDir = fullfile(projectRoot, '0input_data');
    if ~exist(saveDir, 'dir'), mkdir(saveDir); end
    
    % 根据 area_type 的原始输入值来命名文件
    filename = sprintf('EV_%s.xlsx', area_type); % 使用原始 area_type
    writetable(struct2table(evData), fullfile(saveDir, filename));
    disp(['EV数据文件已生成：' filename]);
end

%% 应用实例
% generateEVData(100, 'residential')  % 生成100辆居民区EV
% generateEVData(50, 'workplace', {'特斯拉Model Y'}) % 生成50辆工作区特斯拉