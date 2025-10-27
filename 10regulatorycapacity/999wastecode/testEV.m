%% 基于 initializeEVs 函数的 EV 集群测试（结构体数组版本）%%
clear; close all; clc;

%% 1. 初始化 EV 参数
num_EV = 10;        % EV 数量
rngSeed = 2023;     % 固定随机种子，确保结果可复现
EVs = initializeEVs(num_EV, rngSeed);  % 调用初始化函数

%% 2. 仿真参数设置
T_total = 24;       % 总时间（小时）
dt = 0.5;           % 时间步长（小时）
time_points = 0:dt:T_total; % 时间序列
base_price = 30;    % 基准电价

%% 3. 预分配历史数据存储空间
for i = 1:num_EV
    EVs(i).history.SOC = zeros(1, length(time_points));
    EVs(i).history.P_plus = zeros(1, length(time_points));
    EVs(i).history.P_minus = zeros(1, length(time_points));
end

%% 4. 计算灵活性基准值
E_tar_max = 0.2 * [EVs.C_EV]; % 计算最大目标电量变化
for i = 1:num_EV
    participation = calculateParticipation(EVs(i).p_incentive, base_price);
    [~, ~, deltaE] = incentiveTempEV(EVs(i).p_incentive, 5, 50, 10, 40, E_tar_max(i));
    ptcp_result = (rand() < participation); 
    EVs(i).ptcp = ptcp_result;
    if ptcp_result
        EVs(i).E_tar = EVs(i).E_tar - deltaE;
        if EVs(i).E_tar <= EVs(i).E_in
            EVs(i).E_tar = EVs(i).E_tar + deltaE;
        end
    end
end

%% 5. 提前计算每个 EV 在 24 小时内的基线功率 P_base
H = length(time_points); % 预测时域为 24 小时
for i = 1:num_EV
    [m1, m2, m3] = calculateEVABC_single(... 
        EVs(i).C_EV, EVs(i).eta, EVs(i).E_tar, EVs(i).E_in,...
        EVs(i).t_dep, EVs(i).t_in, dt, EVs(i).r);
    
    [S_opt, P_base_opt] = EVbaseP_single_longstep(m1, m2, m3, EVs(i).p_on, EVs(i).SOC, H);
    EVs(i).P_base_sequence = P_base_opt; % 存储 24 小时的基线功率序列
end

%% 6. 运行主时间循环
for t_idx = 1:length(time_points)
    t = time_points(t_idx);
    fprintf('\n=== 时间 %.1f 小时 ===\n', t);
    
    %% 处理每辆 EV
    for i = 1:num_EV
        %% 6.1 判断 EV 是否在线
        online_mask = (t >= EVs(i).t_in) && (t <= EVs(i).t_dep);
        
        %% 6.2 更新 SOC 状态
        if online_mask
            [~, ~, m3] = calculateEVABC_single(... 
                EVs(i).C_EV, EVs(i).eta, EVs(i).E_tar, EVs(i).E_in,...
                EVs(i).t_dep, EVs(i).t_in, dt, EVs(i).r);
            % EVs(i).m3 = m3;
            [EVs(i).E_exp, EVs(i).E_current, EVs(i).P_current, EVs(i).SOC] = ...
                calculateEVS_single(m3, EVs(i).E_exp, EVs(i).eta,...
                EVs(i).E_current, EVs(i).P_current, EVs(i).C_EV,...
                EVs(i).r, EVs(i).p_on, dt);
            
            % 使用预计算的基线功率
            EVs(i).P_base = EVs(i).P_base_sequence(t_idx);
        else
            EVs(i).P_base = 0;
            EVs(i).P_current = 0;
        end
        
        %% 6.3 计算调节潜力
       [DeltaP_plus, DeltaP_minus] = calculateEVAdjustmentPotentia(...
                EVs(i).C_EV, EVs(i).r,EVs(i).eta, EVs(i).E_tar, EVs(i).E_in,EVs(i).E_current,...
                EVs(i).t_dep, EVs(i).t_in, EVs(i).p_on, 0,...
                EVs(i).P_base, EVs(i).SOC, dt);
        
        %% 6.4 记录历史数据
        EVs(i).history.SOC(t_idx) = EVs(i).SOC;
        EVs(i).history.P_plus(t_idx) = DeltaP_plus;
        EVs(i).history.P_minus(t_idx) = DeltaP_minus;
    end
end

%% 7. 可视化 SOC 变化
for i = 1:num_EV
    figure;
    subplot(2, 1, 1);
    
    % 绘制 SOC 轨迹
    plot(time_points, EVs(i).history.SOC, 'b', 'LineWidth', 1.5);
    title(sprintf('EV%d SOC 演化', i));
    xlabel('时间 (小时)'); 
    ylabel('SOC');
    ylim([-0.1 1.1]);
    grid on;
    
    % 标记充放电时段
    hold on;
    online_period = (time_points >= EVs(i).t_in) & (time_points <= EVs(i).t_dep);
    scatter(time_points(online_period), EVs(i).history.SOC(online_period),...
        20, 'g', 'filled');
    hold off;
    
    legend('SOC 轨迹', '充放电时段', 'Location', 'best');
    
    % 绘制上下调节潜力
    subplot(2, 1, 2);
    plot(time_points, EVs(i).history.P_plus, 'r', 'LineWidth', 1.5);
    hold on;
    plot(time_points, EVs(i).history.P_minus, 'g', 'LineWidth', 1.5);
    title(sprintf('EV%d 上下调节潜力', i));
    xlabel('时间 (小时)');
    ylabel('调节潜力 (kW)');
    legend('上调潜力', '下调潜力', 'Location', 'best');
    grid on;
end

%% 8. 绘制所有 EV 的基线功率 P_base
figure;
hold on;
for i = 1:num_EV
    plot(time_points, EVs(i).P_base_sequence, 'LineWidth', 1.5);
end
title('所有 EV 的基线功率 P_base');
xlabel('时间 (小时)');
ylabel('基线功率 (kW)');
legend(arrayfun(@(i) sprintf('EV%d', i), 1:num_EV, 'UniformOutput', false), 'Location', 'best');
grid on;
hold off;
