clear;clc;close;

% 示例参数（需根据实际数据调整）
params_ev = struct();
params_ev.start_year = 2022;
params_ev.predict_years = 9;  % 2022-2030
params_ev.N0 = 1e6;  % 2022年初始累积购买量
params_ev.M = 1e7 * ones(1,9);  % 假设每年潜在购买者1000万辆
params_ev.p_high = 0.03;  % 高情景创新系数
params_ev.p_low = 0.01;   % 低情景创新系数
params_ev.q_high = 0.35;  % 高情景模仿系数
params_ev.q_low = 0.25;   % 低情景模仿系数
params_ev.r_ac = 0.07;    % 偶然淘汰率7%
params_ev.Sur = [1, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3];  % 前10年留存率

% 调用函数
[ev_high, ev_low] = ev_scale_prediction(params_ev);
% 示例参数（需根据实际数据调整）
params_ac = struct();
params_ac.method = 'logistic';
params_ac.N_hist = [1e7, 1.2e7, 1.5e7, 1.8e7, 2.1e7];  % 2018-2022年历史保有量
params_ac.K_high = 5e7;  % 高情景饱和容量5000万台
params_ac.K_low = 4e7;   % 低情景饱和容量4000万台
params_ac.r_high = 0.15; % 高情景增长率15%
params_ac.r_low = 0.10;  % 低情景增长率10%
params_ac.r_ac = 0.05;   % 偶然淘汰率5%
params_ac.Sur = [1, 0.92, 0.85, 0.78, 0.7, 0.6, 0.5];  % 留存率
params_ac.predict_years = 8;  % 预测2023-2030年

% 调用函数
[ac_high, ac_low] = ac_scale_prediction(params_ac);

num_types = 5;  % 5类充电桩（家庭慢充、家庭快充、工作慢充、工作快充、公共快充）
D = 7;          % 7天数据
N = 144;        % 每天144个10分钟时段

% 生成模拟充电需求CE（单位：kW）
rng(2025);  % 固定随机种子保证可复现
CE = rand(num_types, D, N) * 7;  % 慢充桩功率约7kW（交流桩）
CE(4:5, :, :) = rand(num_types-3, D, N) * 50;  % 快充桩功率约50kW（直流桩）

% 计算每类充电桩的总峰值负荷CE_p（取所有天时段的最大值）
CE_p = max(reshape(CE, num_types, []), [], 2);  % 维度：[5×1]

% 调用函数（使用默认lambda1=0.2，lambda2=1.25）
[HE, ME] = ev_charging_infra_prediction(CE, CE_p);

