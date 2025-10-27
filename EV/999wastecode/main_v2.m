%% 主程序（完全遵循论文模型与main.m结构）
clc; clear; close all;

%% 初始化参数
[EVs, t_sim, dt_short, dt_long, P_tar] = initializeParameters();
H = 12;  % 预测时域
num_steps = t_sim/dt_short;
results = struct(...
    'P_agg',        zeros(1, num_steps), ...
    'P_base',       zeros(1, num_steps), ...
    'S_agg',        zeros(1, num_steps), ...
    'lambda',       zeros(1, num_steps), ...
    'EV_S_original',zeros(1, num_steps));  % 修正字段名称
%% 主循环框架（严格保持main.m结构）
S_agg_current = mean([EVs.S_original]);
% [M1, M2, M3] = calculateAggregationCoefficients(EVs, dt_long);
% [P_base, S_agg] = calculateBasePower_CPLEX(...
%         EVs,...          % 电动汽车集合
%         0.5,...       % 当前时间（分钟）
%         H,...         % 预测时域
%         dt_long,...      % 控制周期
%         M1, M2, M3...   % 聚合系数
%     );
[S_agg,P_agg]=EVbaseP_aggregate(EVs, S_agg_current, H, dt_long);
for t_long = 0:dt_long:t_sim-1
    idx_long = t_long/dt_long + 1;
    
    % 长时间步操作
    [M1, M2, M3] = calculateAggregationCoefficients(EVs, dt_long);
    % [P_base, S_agg] = calculateBasePower_CPLEX(...
    %     EVs,...          % 电动汽车集合
    %     t_long,...       % 当前时间（分钟）
    %     H,...         % 预测时域
    %     dt_long,...      % 控制周期
    %     M1, M2, M3...   % 聚合系数
    % );
    [P_agg, lambda_star, S_agg_next] = aggregateEVs(EVs, P_tar(idx_long), M1, M2, M3, S_agg_current);
    
    % 短时间步循环
    for t_short = t_long:dt_short:t_long+dt_long-1
        idx_short = t_short + 1;
        
        % 更新所有EV状态
        for i = 1:length(EVs)
            EVs(i) = updateLockState(EVs(i), t_short);
            EVs(i) = generateDemandCurve(EVs(i));
            EVs(i).P_current = EVs(i).demandCurve(lambda_star);
            EVs(i) = calculateVirtualSOC(EVs(i), t_short, dt_short);
        end
        
        % 记录结果
        results.P_agg(idx_short) = sum([EVs.P_current]);
        results.lambda(idx_short) = lambda_star;
        results.S_agg(idx_short) = S_agg_current;
        % 更新记录语句
        results.EV_S_original(idx_short) = EVs(5).S_original; % 跟踪第5辆EV的原始SOC

    end
    
    S_agg_current = S_agg_next;
end

%% 专业级可视化（生成4张关键分析图）
time_h = (0:dt_short:t_sim-1)/60;

figure('Position', [100 100 1200 800])

% 图1：功率跟踪对比
subplot(2,2,1)
yyaxis left
plot(time_h, results.P_agg, 'b', 'LineWidth', 1.5)
hold on
stairs(time_h(1:dt_long/dt_short:end), P_tar, 'r--', 'LineWidth', 1.5)
ylabel('功率 (kW)')
yyaxis right
plot(time_h, results.lambda, 'g', 'LineWidth', 1)
ylabel('λ*')
title('功率跟踪与λ*动态')
legend('实际聚合功率','目标功率','λ*值')

% 图2：SOC变化分析
subplot(2,2,2)
plot(time_h, results.S_agg, 'k', 'LineWidth', 1.5)
hold on
plot(time_h, results.EV_S_original, 'm--')
title('SOC动态分析')
ylabel('SOC标幺值')
legend('聚合SOC','EV5个体SOC')

% % 图3：EV状态分布
% subplot(2,2,3)
% state_map = containers.Map({'LockON','LockOFF','ON','OFF'}, [1,2,3,4]);
% state_vec = arrayfun(@(x)state_map(x.state), EVs);
% histogram(state_vec, 'BinMethod','integers')
% set(gca, 'XTickLabel',{'LockON','LockOFF','ON','OFF'})
% title('EV状态分布')
% ylabel('车辆数量')

% % 图4：需求曲线特征
% subplot(2,2,4)
% lambda_test = -1:0.1:1;
% demand_curve = arrayfun(@(ev)ev.demandCurve(lambda_test(11)), EVs(1:5));
% bar(demand_curve)
% title('典型EV需求曲线（λ=0时）')
% ylabel('响应功率 (kW)')
% xlabel('EV编号')

%% 辅助函数（保持与main.m完全一致）
plotResults(results, dt_short)
