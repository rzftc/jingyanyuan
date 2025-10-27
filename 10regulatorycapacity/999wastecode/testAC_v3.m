%% 基于全流程的AC集群测试（结构体数组版本）
clear; close all; clc;

%% 1. 参数初始化
num_AC = 10;     % AC数量
T_total = 24;    % 总时间（小时）
dt = 0.5;        % 时间步长（小时）
rng(2023);       % 随机种子

% 生成AC结构体数组
ACs(num_AC) = struct(); % 预分配结构体数组
for i = 1:num_AC
    ACs(i).R = 1.5 + 1.5*rand();      % 热阻[1.5-3Ω]
    ACs(i).C = 0.3 + 0.4*rand();      % 电容[0.3-0.7F]
    ACs(i).eta = 0.85 + 0.1*rand();   % 效率[0.85-0.95]
    ACs(i).Tset = 22 + 3*randn();     % 初始设定温度
    ACs(i).SOC = 0.5*rand();          % 初始SOC [0-1]
    ACs(i).P_base = 0;                % 基准功率
    ACs(i).Tmax = 26;                 % 温度上限（初始值）
    ACs(i).Tmin = 18;                 % 温度下限（初始值）
    ACs(i).p_incentive = 60 * rand(1);
end

% 电价参数
base_price = 30;                    

% 存储历史数据
history = struct();
history.time = 0:dt:T_total;
history.SOC = zeros(num_AC, length(history.time));
history.P_adj_plus = zeros(num_AC, length(history.time));
history.P_adj_minus = zeros(num_AC, length(history.time));

%% 阶段2: 更新温度设定
p_min = 10; p_max = 50; 
p_min_prime = 15; p_max_prime = 45;
T_set_max = 3; % 最大温度偏移

% 为每个AC单独计算温度变化
for i = 1:num_AC
    participation = calculateParticipation(ACs(i).p_incentive, base_price);
    [~, ~, deltaT] = incentiveTempAC(... 
        ACs(i).p_incentive, p_min, p_max, p_min_prime, p_max_prime, T_set_max);
    
    % 80%概率响应温度变化
    if rand() < participation
        ACs(i).Tmax = ACs(i).Tset + deltaT;
        ACs(i).Tmin = ACs(i).Tset - deltaT;
    end
end

%% 2. 主时间循环
for t_idx = 1:length(history.time)
    t = history.time(t_idx);
    fprintf('\n=== 时间 %.1f小时 ===\n', t);
    
    %% 阶段3: 计算基准功率
    T_ja = [ACs.Tset] + 2*randn(num_AC, 1); % 环境温度
    for i = 1:num_AC
        ACs(i).P_base = ACbaseP_single(... 
            T_ja(i),...
            ACs(i).Tset,...
            ACs(i).R,...
            ACs(i).eta);
    end
    
    %% 阶段4: 更新状态参数
    for i = 1:num_AC
        % 计算ABC系数
        [alpha, beta, gamma] = calculateACABC_single(... 
            ACs(i).R,...
            ACs(i).C,...
            ACs(i).eta,...
            ACs(i).Tmax,...
            ACs(i).Tmin,...
            ACs(i).Tset,...
            dt);
        
        % 保存系数到结构体
        ACs(i).alpha = alpha;
        ACs(i).beta = beta;
        ACs(i).gamma = gamma;
    end
    
    %% 阶段5: SOC更新与归一化
    for i = 1:num_AC
        ACs(i).SOC = ACs(i).alpha * ACs(i).SOC + ... 
                    ACs(i).beta * 0 + ...  % delta_P=0
                    ACs(i).gamma;
        ACs(i).SOC = max(min(ACs(i).SOC, 1), 0); % 归一化
    end
    
    %% 阶段6: 计算调节潜力
    for i = 1:num_AC
        [DeltaP_plus, DeltaP_minus] = calculateACAdjustmentPotentia(... 
            ACs(i).P_base,...
            ACs(i).P_base*2, ...  % Pmax=2*Pbase
            0, ...                % Pmin=0
            ACs(i).alpha,...
            ACs(i).beta,...
            ACs(i).gamma,...
            ACs(i).SOC,...
            dt);
        
        % 存储调节潜力
        history.P_adj_plus(i, t_idx) = DeltaP_plus;
        history.P_adj_minus(i, t_idx) = DeltaP_minus;
    end
    
    %% 存储SOC数据
    history.SOC(:, t_idx) = [ACs.SOC];
end

%% 3. 可视化分析
figure;
subplot(3, 1, 1);
plot(history.time, history.SOC');
title('归一化SOC演化');
xlabel('时间(小时)'); ylabel('SOC');
ylim([0 1]);

subplot(3, 1, 2);
plot(history.time, history.P_adj_plus');
title('上调潜力');
xlabel('时间(小时)'); ylabel('ΔP+ (kW)');

subplot(3, 1, 3);
plot(history.time, history.P_adj_minus');
title('下调潜力');
xlabel('时间(小时)'); ylabel('ΔP- (kW)');

exportgraphics(gcf, 'example_image.png', 'Resolution', 600);
