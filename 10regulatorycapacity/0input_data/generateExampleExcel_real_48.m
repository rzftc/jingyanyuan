%% generateExampleExcel.m
function generateExampleExcel_real_48(num_devices_ac, num_devices_ev, residential_ratio)
    % 参数初始化
    if nargin < 3
        residential_ratio = 0.5; % 默认居民区占比50%
    end
    if nargin < 2
        num_devices_ac = 200;
        num_devices_ev = 100; 
    end

    %% ===== 空调数据生成 =====
    acParams = struct(...
        'Tmin', 18,...
        'Tmax', 26,...
        'R_range', [1.5, 3.0],...
        'C_range', [0.3, 0.7]...
    );
    
    acData = struct(...
        'ID', num2cell(1:num_devices_ac)',...
        'R', num2cell(round(acParams.R_range(1) + diff(acParams.R_range)*rand(num_devices_ac,1), 2)),...
        'C', num2cell(round(acParams.C_range(1) + diff(acParams.C_range)*rand(num_devices_ac,1), 2)),...
        'eta', num2cell(round(0.85 + 0.1*rand(num_devices_ac,1), 2)),...
        'Tset', num2cell(round(acParams.Tmin + (acParams.Tmax-acParams.Tmin)*rand(num_devices_ac,1), 1)),...
        'Tmax', num2cell(repmat(acParams.Tmax, num_devices_ac,1)),...
        'Tmin', num2cell(repmat(acParams.Tmin, num_devices_ac,1)),...
        'SOC', num2cell(round(0.5*rand(num_devices_ac,1), 2)),...
        'p_incentive', num2cell(round(60*rand(num_devices_ac,1), 1))...
    );

    %% ===== 电动汽车数据生成（重大更新）=====
    % 精确车型数据库（2024年）
    evModels = [
        % 名称                 慢充功率  快充功率  电池容量选项
        struct('name','比亚迪秦PLUS',  'p_slow',6.6, 'p_fast',60,  'C_EV',[18.32,26.86], 'prob',0.35),...
        struct('name','比亚迪海鸥',    'p_slow',6.6, 'p_fast',40,  'C_EV',[30.08,38.88], 'prob',0.25),...
        struct('name','特斯拉Model Y','p_slow',11,  'p_fast',250, 'C_EV',[60,78.4],    'prob',0.15),...
        struct('name','比亚迪元PLUS', 'p_slow',7,   'p_fast',90,  'C_EV',[50.12,60.48],'prob',0.25)
    ];

    % 分区域时间分布（居民区 vs 工作区）
    timeSlots = {
        % 0-4点
        struct('range',0:4, 'residential',[70,20,10], 'workplace',[10,20,70]), 
        % 5-8点
        struct('range',5:8, 'residential',[20,70,10], 'workplace',[20,10,70]),
        % 9-15点
        struct('range',9:15,'residential',[10,20,70], 'workplace',[20,70,10]),
        % 16-17点
        struct('range',16:17,'residential',[20,70,10],'workplace',[10,20,70]),
        % 18-23点
        struct('range',18:23,'residential',[70,20,10],'workplace',[10,20,70])
    };

    % 初始化存储数组
    [C_EV,eta,p_on,t_in,t_dep,E_in,E_tar] = deal(zeros(num_devices_ev,1));
    areaType = cell(num_devices_ev,1); % 改为元胞数组存储字符串
    
    %% 主生成循环
    for i = 1:num_devices_ev
        % 确定区域类型（直接存储中文）
        if rand() < residential_ratio
            areaType{i} = '居民区';
            isResidential = true;
        else
            areaType{i} = '工作区';
            isResidential = false;
        end
        
        % 车型选择
        cumProb = cumsum([evModels.prob]);
        selectedModel = evModels(find(rand() <= cumProb,1));
        
        % 设置充电参数
        C_EV(i) = selectedModel.C_EV(randi(length(selectedModel.C_EV))); % 精确电池容量
        p_on(i) = isResidential*selectedModel.p_slow + (~isResidential)*selectedModel.p_fast;
        eta(i) = 0.85 + 0.1*rand(); % 充电效率85%-95%
        
        % 时间参数生成
        valid = false;
        while ~valid
            % 生成接入时间（允许0-47小时）
            t_in_h = randi([0,47]);  % 修改范围
            t_in_m = randi([0,59]);
            t_in(i) = t_in_h + t_in_m/60;
            
            % 获取时间分布（自动模24处理）
            current_hour = mod(t_in_h,24);
            currentSlot = timeSlots{find(cellfun(@(x) ismember(current_hour,x.range), timeSlots),1)};
           % 修正区域分布判断语法
            if isResidential
                dist = currentSlot.residential;
            else
                dist = currentSlot.workplace;
            end
            
            % 选择充电类型
            randVal = rand()*100;
            if randVal <= dist(1)
                chargeType = 1; % 长时充电
            elseif randVal <= sum(dist(1:2))
                chargeType = 2; % 中时充电
            else
                chargeType = 3; % 短时充电
            end
            
            % 设置充电时长限制
            switch chargeType
                case 1 % 长时：9-24小时
                    minDur = 9; maxDur = 24;
                case 2 % 中时：4-8小时
                    minDur = 4; maxDur = 8;
                case 3 % 短时：0-3小时
                    minDur = 0; maxDur = 3;
            end
            
             % 生成持续时间（允许超过24小时）
            duration = minDur + (maxDur - minDur)*rand();
            t_dep(i) = t_in(i) + duration;  % 直接累加不取模
            
            % 电量参数生成
            E_in(i) = 0.1*C_EV(i) + 0.3*C_EV(i)*rand(); % 初始电量10%-40%
            max_charge = p_on(i) * eta(i) * duration;
            E_tar_max = min(E_in(i) + max_charge, 0.95*C_EV(i));
            
            % 非线性充电曲线
            E_tar(i) = E_in(i) + (E_tar_max - E_in(i)) * (0.2 + 0.8*rand()^2);
            
            % 验证充电功率
           % 验证条件简化
            if t_dep(i) > t_in(i) && E_tar(i) > E_in(i)
                valid = true;
            end
        end
        
        % 数据精度处理
        C_EV(i) = round(C_EV(i), 2);
        eta(i) = round(eta(i), 2);
        t_in(i) = round(t_in(i), 2);
        t_dep(i) = round(t_dep(i), 2);
        E_in(i) = round(E_in(i), 1);
        E_tar(i) = round(E_tar(i), 1);
    end

    %% 构建数据结构（添加区域类型）
    evData = struct(...
        'ID', num2cell(1:num_devices_ev)',...
        'C_EV', num2cell(C_EV),...
        'eta', num2cell(eta),...
        'p_on', num2cell(p_on),...
        'E_in', num2cell(E_in),...
        'E_tar', num2cell(E_tar),...
        't_in', num2cell(t_in),...
        't_dep', num2cell(t_dep),...
        'AreaType', areaType,... % 直接使用元胞数组
        'r', num2cell(repmat(0.025, num_devices_ev,1)),...
        'SOC', num2cell(round(rand(num_devices_ev,1), 2)),...
        'p_incentive', num2cell(round(60*rand(num_devices_ev,1), 1))...
    );

     % 获取当前函数所在路径的上一级目录（项目根目录）
    currentDir = fileparts(mfilename('fullpath'));
    projectRoot = fullfile(currentDir, '..');
    % 构建目标保存路径
    saveDir = fullfile(projectRoot, '0input_data');
    % 确保目录存在
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    ac_file = fullfile(saveDir, 'AC_template.xlsx');
    ev_file = fullfile(saveDir, 'EV_template.xlsx');
    % ===== 写入文件 =====
    writetable(struct2table(acData), ac_file);
    writetable(struct2table(evData), ev_file);
    disp('数据文件已生成: AC_template.xlsx 和 EV_template.xlsx');
end
