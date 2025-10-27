clc; clear; close all;
rng(2024);

%% ===================== 初始化参数 =====================
% 生成/加载EV参数文件
excelFile = '../0inputdata/residential_all_models.xlsx';
if ~exist(excelFile, 'file')
    generateEVParameters_real(excelFile, 100, 0.6); % 自动生成100辆EV数据
    fprintf('已生成参数模板: %s\n', excelFile);
end

% 从Excel加载参数
[EVs, t_sim, dt_short, dt_long, P_tar] = initializeFromExcel(excelFile);
fprintf('成功加载%d辆EV数据\n', length(EVs));

% plotDeltaEComparison(EVs)
%% ===================== 索引预处理 =====================
num_long_steps = t_sim / dt_long;        % 144个长步
num_short_per_long = dt_long / dt_short; % 100个短步/长步
total_steps = num_long_steps * num_short_per_long; % 14400
hours_per_step = round(24 / num_long_steps, 2);
assert(mod(num_long_steps, 24) == 0, '长时间步数必须是24的整数倍');
repeat_factor = num_long_steps / 24; % 每小时的步数

%% ===================== 创建结果存储结构 =====================
results = struct(...
    'P_agg',        zeros(1, total_steps), ...  % 聚合功率
    'P_base',       zeros(1, total_steps), ...  % 基准功率
    'S_agg',        zeros(1, total_steps), ...  % 虚拟SOC
    'lambda',       zeros(1, total_steps), ...  % λ*
    'EV_S_original', zeros(length(EVs), total_steps)... % 存储所有EV的SOC
);

%% ===================== 初始化EV目标充电功率 =====================
for i = 1:length(EVs)
    Delta_E = calculateDeltaE(EVs(i), EVs(i).p_real);
    EVs(i).E_tar = max(EVs(i).E_tar_set - Delta_E, EVs(i).E_ini); % 电量下限保护
    t_ch = 60 * (EVs(i).E_tar - EVs(i).E_ini) / EVs(i).P_N;
    EVs(i).tau_rem = t_ch;
end

%% ===================== 初始化聚合SOC =====================
S_agg_current = mean([EVs.S_original]);
[S_agg_opt, P_base_opt] = EVbaseP_aggregate_short(EVs, S_agg_current, 48, dt_long);

%% ===================== 外层循环（长时间步长） =====================
for long_idx = 1:num_long_steps
    t_long = (long_idx - 1) * dt_long; % 当前长步起始时间
    
    %% --------------------- 长时间步处理 ---------------------
    % 阶段3: 执行需求响应
    EVs = distributeBasePower(EVs, 1,P_base_opt);  % 新增调用
    [lambda_star] = aggregateEVs(EVs, P_tar(long_idx));
    [P_agg, S_agg_next] = calculateVirtualSOC_agg(EVs, dt_long);
    
    %% ===================== 内层循环（短时间步长） =====================
    for short_idx = 1:num_short_per_long
        step_idx = (long_idx - 1) * num_short_per_long + short_idx;
        t_current = t_long + (short_idx - 1) * dt_short;
        
        %% --------------------- 更新EV状态 ---------------------
        for i = 1:length(EVs)
            EV = EVs(i);
            % 状态更新
            EV = updateLockState(EV, t_current);
            EV = generateDemandCurve(EV);
            EV.P_current = EV.demandCurve(lambda_star);
            EV = calculateVirtualSOC_upgrade(EV, t_current, dt_short);
            m3 = (EV.E_tar - EV.E_ini) / (EV.eta * ((EV.t_dep - EV.t_in) / 60));
            results.m3(i) = m3;
            EVs(i) = EV;
        end
        
        %% --------------------- 记录结果 ---------------------
        results.lambda(step_idx) = lambda_star;
        results.S_agg(step_idx) = S_agg_current;
        results.P_agg(step_idx) = P_agg;
       
        results.P_tar = repelem(P_tar, num_short_per_long); % 确保与P_agg同维度
        results.P_cu(step_idx) = EVs(10).P_current;
        % 记录所有EV的SOC
        for ev_idx = 1:length(EVs)
            results.EV_S_original(ev_idx, step_idx) = EVs(ev_idx).S_original;
            results.EV_S_mod(ev_idx, step_idx) = EVs(ev_idx).S_modified;
        end
    end
    
    %% --------------------- 更新聚合SOC ---------------------
    S_agg_current = S_agg_next;
end
%% ===================== 可视化结果 =====================

%% 绘制第10台EV的SOC与lambda_star对比图
%% 绘制第10台EV的SOC与P_current对比图
selected_ev = 10; % 选择第10台EV
time_min = (0:total_steps-1) * dt_short; % 时间轴（分钟）
time_h = time_min / 60; % 转换为小时（更直观）

% 提取数据
S_original = results.EV_S_original(selected_ev, :);
S_mod = results.EV_S_mod(selected_ev, :);
P_current = results.lambda; % 第10台EV的当前功率

% 创建双轴图
figure('Position', [100 100 1200 600]); % 设置图窗大小

% 左轴：绘制SOC（原始和修正）
yyaxis left;
plot(time_h, S_original, 'b-', 'LineWidth', 1.5, 'DisplayName', '原始SOC (S_original)');
hold on;
plot(time_h, S_mod, 'r--', 'LineWidth', 1.5, 'DisplayName', '修正SOC (S_mod)');
ylabel('SOC');
ylim([-1.1 1.1]); % SOC理论范围[-1,1]，扩展显示边界
grid on;

% 右轴：绘制P_current（当前功率）
yyaxis right;
plot(time_h, P_current, 'g-.', 'LineWidth', 1.5, 'DisplayName', '当前功率 (P_current)');
ylabel('功率 (kW)');
% 根据功率实际范围调整y轴（自动适配）
P_lim = [min(P_current)-0.1, max(P_current)+0.1];
ylim(P_lim);

% 通用设置
xlabel('时间 (小时)');
title(['第', num2str(selected_ev), '台EV的SOC与当前功率对比']);
legend('Location', 'best'); % 图例位置自动优化
set(gca, 'FontSize', 12); % 设置坐标轴字体大小
hold off;

% 生成三个分析视图
% plotLambdaAndSOC(results, dt_short, selected_ev)    % 图1
% plotPowerComparison(results, dt_short)              % 图2
% plotLambdaAndAggSOC(results, dt_short)               % 图3

