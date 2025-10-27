function generateEVParameters_real_choose(filePath, numEV, areaRatio, varargin)
    % 增强版EV参数生成函数（支持车型筛选与区域类型指定）
    % 输入参数：
    %   filePath    - 输出文件路径（必选）
    %   numEV       - EV总数（默认1000）
    %   areaRatio   - 居民区比例（默认0.5，当areaType='混合'时生效）
    %   varargin    - 键值对参数：
    %       'ModelNames'  : 指定车型名称（默认全部车型）
    %       'AreaType'    : 区域类型['居民区','工作区','混合']（默认'混合'）
    %% 使用范例
    % 生成200辆工作区特斯拉
    % generateEVParameters_real_choose('tesla_work.xlsx', 200, 0,...
    %     'ModelNames', {'特斯拉Model Y'},...
    %     'AreaType', '工作区');
    % 
    % % 生成500辆居民区比亚迪
    % generateEVParameters_real_choose('byd_res.xlsx', 500, 1,...
    %     'ModelNames', {'比亚迪秦PLUS DM-i','比亚迪海鸥'},...
    %     'AreaType', '居民区');
    % % 生成1000辆混合区域车辆（30%居民区+70%工作区）
    % generateEVParameters_real_choose('mixed_area.xlsx', 1000, 0.3,...
    %     'AreaType', '混合');
   %  generateEVParameters_real_choose('resi_inc_1000.xlsx', 1000, 0,...
   % 'AreaType', '居民区');

    %% 参数解析系统
    p = inputParser;
    addRequired(p, 'filePath', @ischar);
    addOptional(p, 'numEV', 1000, @(x) x>0 && mod(x,1)==0);
    addOptional(p, 'areaRatio', 0.5, @(x) x>=0 && x<=1);
    addParameter(p, 'ModelNames', {}, @iscellstr);
    addParameter(p, 'AreaType', '混合', @(x) any(strcmp(x,{'居民区','工作区','混合'})));
    parse(p, filePath, numEV, areaRatio, varargin{:});
    
    %% 精确车型参数库（增加快充功率参数）
    models = [
        struct('Name','比亚迪秦PLUS DM-i', 'C',[18.32,26.86], 'SlowCharge',6.6, 'FastCharge',60,  'Ratio',0.28),...
        struct('Name','比亚迪海鸥',        'C',[30.08,38.88], 'SlowCharge',6.6, 'FastCharge',40,  'Ratio',0.25),...
        struct('Name','比亚迪宋PLUS DM-i','C',[18.32,26.86], 'SlowCharge',6.6, 'FastCharge',70,  'Ratio',0.20),...
        struct('Name','特斯拉Model Y',   'C',[60.0,78.4],   'SlowCharge',11,  'FastCharge',250, 'Ratio',0.15),...
        struct('Name','比亚迪元PLUS',     'C',[50.12,60.48],'SlowCharge',7,   'FastCharge',90,  'Ratio',0.12)...
    ];
    
    %% 车型筛选逻辑
    if ~isempty(p.Results.ModelNames)
        % 验证输入车型名称
        validModels = ismember({models.Name}, p.Results.ModelNames);
        if ~any(validModels)
            error('无效车型名称，可用车型: %s', strjoin({models.Name}, ', '));
        end
        models = models(validModels);
        % 重新归一化比例
        totalRatio = sum([models.Ratio]);
        for i = 1:length(models)
            models(i).Ratio = models(i).Ratio / totalRatio;
        end
    end

    %% 区域分配增强系统
    switch p.Results.AreaType
        case '居民区'
            isResidential = true(numEV,1);
        case '工作区'
            isResidential = false(numEV,1);
        case '混合'
            isResidential = rand(numEV,1) < p.Results.areaRatio;
        otherwise
            error('未知区域类型: %s', p.Results.AreaType);
    end
    
    %% 车型分配系统
    model_probs = cumsum([models.Ratio]);
    model_idx = arrayfun(@(x) find(x <= model_probs,1), rand(numEV,1));

    %% 初始化数据表（完整字段）
    data = table();
    data.EV_ID = (1:numEV)';
    data.Area = repmat("", numEV,1);
    
    %% 严格参数绑定
    for i = 1:numEV
        m = models(model_idx(i));
        data.C(i) = m.C(randi(numel(m.C))); % 电池容量二选一
        
        % 根据区域选择充电功率
        if isResidential(i)
            data.P_N(i) = m.SlowCharge;
            data.Area(i) = "居民区";
        else
            data.P_N(i) = m.FastCharge;
            data.Area(i) = "工作区";
        end
    end
    
    %% 电价参数系统（保持不变）
    basePrice = 0.5; % 基准电价元/kWh
    priceRange = 0.3; % 浮动范围
    data.p_incentive = round(60*rand(numEV,1), 1); % This is the line for p_incentive
    data.P_0 = basePrice + priceRange*randn(numEV,1)*0.2;
    data.P_h_max = data.P_0 + priceRange*rand(numEV,1);
    data.P_l_min = data.P_0 - priceRange*rand(numEV,1);
    
    data.P_l_min = max(data.P_l_min, 0.2);           % 最低电价保护
    data.P_h_max = max(data.P_h_max, data.P_0 + 0.1); % 合理价差保障
    
    data.Delta_E_h_max = 0.1 * data.C;  % 高功率区能量变化
    data.Delta_E_q_max = 0.05 * data.C; % 低功率区能量变化
    
    price_std = 0.1;
    data.p_real = data.P_0 + price_std * randn(numEV,1);

    %% 其他参数系统（保持不变）
    data.eta = 0.85 + 0.1*rand(numEV,1);  % 充电效率
    data.r = 0.05 * ones(numEV,1);        % SOC调节系数
    data.state = repmat("LockOFF", numEV,1);  % 初始状态

    %% 电量参数生成
    data.E_ini = data.C .* (0.1 + 0.3*rand(numEV,1));       % 初始电量10%-40%
    data.E_tar_set = data.C .* (0.8 + 0.1999*rand(numEV,1)); % 目标电量80%-99.99%

    %% 时间生成系统（增强约束验证）
    hour_edges = cumsum([0,5,4,7,2,6])/24; % [0,5,9,16,18,24]小时
    [~, hour_groups] = histc(rand(numEV,1), hour_edges);
    time_labels = {'0-4','5-8','9-15','16-17','18-23'};
    
    residential_weights = [
        0.70, 0.20, 0.10;   % 0-4时
        0.20, 0.70, 0.10;   % 5-8时
        0.10, 0.20, 0.70;   % 9-15时
        0.20, 0.70, 0.10;   % 16-17时
        0.70, 0.20, 0.10    % 18-23时
    ];
    
    workplace_weights = [
        0.10, 0.20, 0.70;   % 0-4时
        0.20, 0.10, 0.70;   % 5-8时
        0.20, 0.70, 0.10;   % 9-15时
        0.10, 0.20, 0.70;   % 16-17时
        0.10, 0.20, 0.70    % 18-23时
    ];

    data.t_in = zeros(numEV,1);
    data.t_dep = zeros(numEV,1);
    
    for i = 1:numEV
        group = hour_groups(i);
        time_range = sscanf(time_labels{group},'%d-%d')*60;
        
        % 入网时间生成
        data.t_in(i) = randi([time_range(1), time_range(2)]);
        
        % 区域权重选择
        if strcmp(data.Area(i), "居民区")
            time_weights = residential_weights(group,:);
        else
            time_weights = workplace_weights(group,:);
        end
        
        % 参数约束验证循环
        valid = false;
        while ~valid
            % 停车时长选择
            dur_probs = cumsum(time_weights);
            dur_type = find(rand() <= dur_probs, 1);
            
            % 持续时间生成
            switch dur_type
                case 1 % 长时停车
                    dep_add = randi([9*60, 24*60]);
                case 2 % 中时停车
                    dep_add = randi([4*60, 8*60]);
                case 3 % 短时停车
                    dep_add = randi([0, 3*60]);
            end
            
            % 计算最小充电时间
            delta_E = data.E_tar_set(i) - data.E_ini(i);
            t_min = (delta_E/(data.eta(i)*data.P_N(i)))*60;
            dep_add = max(dep_add, ceil(t_min) + 1); % 增加1分钟缓冲
            
            % 验证功率约束
            m3 = delta_E/(data.eta(i)*(dep_add/60));
            valid = (data.P_N(i) > m3);
            
            % 异常处理
            if dep_add > 24*60*3 % 超过3天
                error('EV%d生成异常: dep_add=%d', i, dep_add);
            end
        end
        
        % 处理跨天逻辑
        data.t_dep(i) = data.t_in(i) + dep_add;
        data.t_dep(i) = data.t_dep(i) + 1440 * floor(data.t_in(i)/1440);
        
        % 时间顺序验证
        assert(data.t_dep(i) > data.t_in(i), '时间异常 EV%d: t_in=%d, t_dep=%d', i, data.t_in(i), data.t_dep(i));
    end

    %% 数据完整性验证
    valid_C = [18.32,26.86,30.08,38.88,50.12,60.0,60.48,78.4];
    valid_PN = [6.6,7,11,40,60,70,90,250];
    assert(all(ismember(data.C,valid_C)),'电池容量校验失败');
    assert(all(ismember(data.P_N,valid_PN)),'充电功率校验失败');
    assert(all(data.t_dep > data.t_in),'时间顺序校验失败');
    
    % 功率约束最终验证
    m3_all = (data.E_tar_set - data.E_ini)./(data.eta.*(data.t_dep - data.t_in)/60);
    assert(all(data.P_N > m3_all),'存在不满足P_N>m3的情况');

    %% 数据存储
    writetable(data, filePath);
    fprintf('EV数据生成完成：%s\n居民区：%d辆（%.1f%%）\n工作区：%d辆（%.1f%%）\n',...
            filePath,...
            sum(isResidential), 100*mean(isResidential),...
            sum(~isResidential), 100*mean(~isResidential));
end