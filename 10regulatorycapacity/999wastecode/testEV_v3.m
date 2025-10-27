%% 基于全流程的EV集群测试（时间循环版，结构体数组版本）
clear; close all; clc;

%% 1. 初始化参数（使用结构体数组 EVs(i)）
num_EV = 10;    % EV数量
T_total = 24;   % 总时间（小时）
dt = 0.5;       % 时间步长（小时）
rng(2023);      % 固定随机种子

% 生成EV参数（结构体数组初始化）
EVs = struct('C_EV', {}, 'eta', {}, 'p_on', {}, 'E_in', {}, 'E_tar', {}, ...
             't_in', {}, 't_dep', {}, 'r', {}, 'E_exp', {}, 'P_current', {}, ...
             'SOC', {}, 'P_base', {}, 'E_current', {});

for i = 1:num_EV
    EVs(i).C_EV = 40 + 60 * rand();          % 电池容量[40-100kWh]
    EVs(i).eta = 0.8 + 0.15 * rand();        % 充放电效率[0.8-0.95]
    EVs(i).p_on = 19 + 3 * rand();           % 额定功率[3-22kW]
    EVs(i).E_in = 10 + 30 * rand();          % 初始电量[10-40kWh]
    EVs(i).E_tar = EVs(i).E_in + (EVs(i).C_EV - EVs(i).E_in) * rand(); % 目标电量
    EVs(i).t_in = fix(8 + 2 * randn());      % 入网时间[6-10点]
    EVs(i).t_dep = EVs(i).t_in + 2 + randi(6); % 离网时间
    EVs(i).r = 0.025;                      % SOC计算系数
    
    % 初始化状态变量
    EVs(i).E_exp = EVs(i).E_in;
    EVs(i).P_current = 0;
    EVs(i).SOC = 0;
    EVs(i).P_base = 0;
    EVs(i).E_current = EVs(i).E_in;
    EVs(i).p_incentive = 60 * rand(1);
end

% 电价参数
base_price = 30;                
% 存储历史数据（改为结构体数组字段）
for i = 1:num_EV
    EVs(i).history.SOC = zeros(1, T_total / dt + 1);
    EVs(i).history.P_plus = zeros(1, T_total / dt + 1);
    EVs(i).history.P_minus = zeros(1, T_total / dt + 1);
end
time_points = 0:dt:T_total;

%% 阶段2: 更新灵活性基准值
E_tar_max = 0.2 * [EVs.C_EV]; % 提取所有EV的容量
% 更新每个EV的目标电量
for i = 1:num_EV
    participation = calculateParticipation(EVs(i).p_incentive, base_price);
    
    [~, ~, deltaE] = incentiveTempEV(...
        EVs(i).p_incentive, 5, 50, 10, 40, E_tar_max(i));
    if rand() < participation
        EVs(i).E_tar = EVs(i).E_tar - deltaE;
    end
end

%% 2. 主时间循环（逐个EV操作）
for t_idx = 1:length(time_points)
    t = time_points(t_idx);
    fprintf('\n=== 时间 %.1f小时 ===\n', t);
    
    %% 阶段3: 计算基线功率
    for i = 1:num_EV
        if t >= EVs(i).t_in && t <= EVs(i).t_dep
            [m1, m2, m3] = calculateEVABC_single(...
                EVs(i).C_EV, EVs(i).eta, EVs(i).E_tar, EVs(i).E_in, ...
                EVs(i).t_dep, EVs(i).t_in, dt, EVs(i).r);
            
            [~, EVs(i).P_base] = EVbaseP_single(m1, m2, m3, EVs(i).p_on, EVs(i).SOC);
        else
            EVs(i).P_base = 0;
            EVs(i).P_current = 0;
        end
    end
    
    %% 阶段4: 更新SOC
    for i = 1:num_EV
        if t >= EVs(i).t_in && t <= EVs(i).t_dep
            [~, ~, m3] = calculateEVABC_single(...
                EVs(i).C_EV, EVs(i).eta, EVs(i).E_tar, EVs(i).E_in, ...
                EVs(i).t_dep, EVs(i).t_in, dt, EVs(i).r);
            
            [EVs(i).E_exp, EVs(i).E_current, EVs(i).P_current, EVs(i).SOC] = ...
                calculateEVS_single(m3, EVs(i).E_exp, EVs(i).eta, ...
                EVs(i).E_current, EVs(i).P_current, EVs(i).C_EV, ...
                EVs(i).r, EVs(i).p_on, dt);
        end
    end
    
    %% 阶段5: 计算调节潜力
    for i = 1:num_EV
        [DeltaP_plus, DeltaP_minus] = calculateEVAdjustmentPotentia(...
            EVs(i).C_EV, EVs(i).eta, EVs(i).E_tar, EVs(i).E_in, ...
            EVs(i).t_dep, EVs(i).t_in, ...
            EVs(i).p_on, 0, ... % Pmin=0
            EVs(i).P_base, EVs(i).SOC, dt);
        
        % 存储历史数据
        EVs(i).history.SOC(t_idx) = EVs(i).SOC;
        EVs(i).history.P_plus(t_idx) = DeltaP_plus;
        EVs(i).history.P_minus(t_idx) = DeltaP_minus;
    end
end

%% 3. 可视化（按EV绘制）
figure;
for i = 1:num_EV
    subplot(num_EV, 1, i);
    plot(time_points, EVs(i).history.SOC);
    title(sprintf('EV %d 的SOC演化', i));
    xlabel('时间(小时)'); ylim([0, 1]);
end