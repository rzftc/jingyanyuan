%% generateResidentialEV.m
function generateResidentialEV(num_devices_ev)
    % 参数初始化
    if nargin < 1
        num_devices_ev = 100;  % 默认EV数量
    end

    %% ===== 电动汽车数据生成 =====
    % 完整车型数据库（保留全部原始车型）
    evModels = [
        struct('name','比亚迪秦PLUS',  'p_slow',6.6, 'p_fast',60,  'C_EV',[18.32,26.86], 'prob',0.35),...
        struct('name','比亚迪海鸥',    'p_slow',6.6, 'p_fast',40,  'C_EV',[30.08,38.88], 'prob',0.25),...
        struct('name','特斯拉Model Y','p_slow',11,  'p_fast',250, 'C_EV',[60,78.4],    'prob',0.15),...
        struct('name','比亚迪元PLUS', 'p_slow',7,   'p_fast',90,  'C_EV',[50.12,60.48],'prob',0.15),...
        struct('name','比亚迪宋PLUS DM-i', 'p_slow',6.6, 'p_fast',70,  'C_EV',[18.32,26.86],'prob',0.1),...
    ];

    % 居民区专用时间分布参数
    timeSlots = {
        struct('range',0:4,  'dist',[70,20,10], 'max_dur',[24,8,3]),    % 0-4点
        struct('range',5:8,  'dist',[20,70,10], 'max_dur',[19,8,3]),   % 5-8点
        struct('range',9:15, 'dist',[10,20,70], 'max_dur',[15,8,3]),  % 9-15点
        struct('range',16:17,'dist',[20,70,10], 'max_dur',[8,8,3]),    % 16-17点
        struct('range',18:23,'dist',[70,20,10], 'max_dur',[6,8,3])     % 18-23点
    };

    % 初始化存储数组
    [C_EV,eta,p_on,t_in,t_dep,E_in,E_tar] = deal(zeros(num_devices_ev,1));
    areaType = repmat({'居民区'}, num_devices_ev, 1);  % 固定区域类型

    %% 主生成循环
    for i = 1:num_devices_ev
        % 车型选择（概率归一化处理）
        cumProb = cumsum([evModels.prob]/sum([evModels.prob]));
        selectedModel = evModels(find(rand() <= cumProb,1));

        % 设置充电参数（固定使用慢充功率）
        C_EV(i) = selectedModel.C_EV(randi(length(selectedModel.C_EV)));
        p_on(i) = selectedModel.p_slow;  % 居民区仅使用慢充
        eta(i) = 0.85 + 0.1*rand();

        % 时间参数生成
        valid = false;
        while ~valid
            % 生成接入时间
            t_in_h = randi([0,23]);
            t_in_m = randi([0,59]);
            t_in(i) = t_in_h + t_in_m/60;
            
            % 匹配时间槽
            currentSlot = timeSlots{find(cellfun(@(x) ismember(t_in_h,x.range), timeSlots),1)};
            
            % 选择充电类型
            randVal = rand()*100;
            if randVal <= currentSlot.dist(1)
                chargeType = 1; minDur = 9; maxDur = currentSlot.max_dur(1);
            elseif randVal <= sum(currentSlot.dist(1:2))
                chargeType = 2; minDur = 4; maxDur = currentSlot.max_dur(2);
            else
                chargeType = 3; minDur = 0; maxDur = currentSlot.max_dur(3);
            end

            % 持续时间计算
            duration = min(maxDur, 24 - t_in(i));
            if duration < minDur, continue; end
            t_dep(i) = t_in(i) + duration;

            % 电量参数生成
            E_in(i) = 0.1*C_EV(i) + 0.3*C_EV(i)*rand();
            max_charge = p_on(i) * eta(i) * duration;
            E_tar_max = min(E_in(i) + max_charge, 0.95*C_EV(i));
            E_tar(i) = E_in(i) + (E_tar_max - E_in(i)) * (0.2 + 0.8*rand()^2);

            % 验证条件
            required_power = (E_tar(i) - E_in(i)) / (eta(i) * duration);
            valid = (t_dep(i) <= 24) && (required_power <= p_on(i)) && (E_tar(i) > E_in(i));
        end

        % 数据精度处理
        variables = {'C_EV','eta','t_in','t_dep','E_in','E_tar'};
        for var = variables
            eval(sprintf('%s(i) = round(%s(i), 2);', var{1}, var{1}));
        end
    end

    %% 构建数据结构
    evData = struct(...
        'ID', num2cell(1:num_devices_ev)',...
        'C_EV', num2cell(C_EV),...
        'eta', num2cell(eta),...
        'p_on', num2cell(p_on),...
        'E_in', num2cell(E_in),...
        'E_tar', num2cell(E_tar),...
        't_in', num2cell(t_in),...
        't_dep', num2cell(t_dep),...
        'AreaType', areaType,...
        'r', num2cell(repmat(0.025, num_devices_ev,1)),...
        'SOC', num2cell(round(rand(num_devices_ev,1), 2)),...
        'p_incentive', num2cell(round(60*rand(num_devices_ev,1), 1))...
    );

    %% 文件保存
    currentDir = fileparts(mfilename('fullpath'));
    saveDir = fullfile(currentDir, '..', '0input_data');
    if ~exist(saveDir, 'dir'), mkdir(saveDir); end
    
    ev_file = fullfile(saveDir, 'EV_Residential.xlsx');
    writetable(struct2table(evData), ev_file);
    disp(['居民区EV数据文件已生成：' ev_file]);
end
