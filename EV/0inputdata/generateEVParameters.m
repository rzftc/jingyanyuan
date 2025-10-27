function generateEVParameters(filePath, numEV)
    % 生成电价相关参数模板
    if nargin < 2
        numEV = 10;
    end

    % 电价参数范围（元/kWh）
    basePrice = 0.5; % 基准电价
    priceRange = 0.3; % 电价浮动范围
    
    % 创建数据表
    data = table();
    
    %% 基础参数
    data.EV_ID = (1:numEV)';
    data.P_N = 6 + 4*rand(numEV,1);       % 额定功率 [6,10] kW

    %% 电池参数（保持不变）
    data.C = 50 + 20*rand(numEV,1);       % 电池容量
    data.eta = 0.85 + 0.1*rand(numEV,1);  % 效率
    data.r = 0.05*ones(numEV,1);          % SOC调节系数
    
    %% 电价参数
    data.P_0 = basePrice + priceRange*randn(numEV,1)*0.2; % 基准电价
    data.P_h_max = data.P_0 + priceRange*rand(numEV,1);   % 最高电价阈值
    data.P_l_min = data.P_0 - priceRange*rand(numEV,1);   % 最低电价阈值
    
    % 确保电价关系 P_l_min < P_0 < P_h_max
    data.P_l_min = max(data.P_l_min, 0.2); % 最低电价不低于0.2元
    data.P_h_max = max(data.P_h_max, data.P_0 + 0.1); % 保持合理价差
    data.Delta_E_h_max = 0.1 * data.C;  % 高功率区能量变化
    data.Delta_E_q_max = 0.05 * data.C; % 低功率区能量变化
    price_std = 0.1; % 电价标准差
    data.p_real = data.P_0 + price_std * randn(numEV,1); % 实时电价初始化

    %% 时间参数（保持不变）
    data.t_in = randi([1, 6*60], numEV,1);
    data.t_dep = data.t_in + 8*60 + randi([-60,60],numEV,1);
    
    %% 电量参数（保持不变）
    data.E_ini = data.C .* (0.1 + 0.2*rand(numEV,1));
    data.E_tar_set = data.E_ini + 0.3*data.C + 0.55*data.C.*rand(numEV,1);
    
    %% 状态参数
    data.state = repmat("OFF", numEV,1);
    
    % 保存文件
    writetable(data, filePath);
    
end


