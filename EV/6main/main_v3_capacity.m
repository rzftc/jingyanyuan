%% mainv3.m - 基于Excel数据接口的完整版
clc; clear; close all;
rng(2024);
%% ===================== 初始化参数 =====================
% 生成/加载EV参数文件
excelFile = '../0inputdata/evdata.xlsx';
if ~exist(excelFile, 'file')
    generateEVParameters_real(excelFile, 100, 0.6); % 自动生成100辆EV数据
    fprintf('已生成参数模板: %s\n', excelFile);
end

% 从Excel加载参数
[EVs, t_sim, dt_short, dt_long, P_tar] = initializeFromExcel(excelFile);
fprintf('成功加载%d辆EV数据\n', length(EVs));

%% ===================== 索引预处理 =====================
num_long_steps = t_sim / dt_long;        % 144个长步
num_short_per_long = dt_long / dt_short; % 100个短步/长步
total_steps = num_long_steps * num_short_per_long; % 14400
hours_per_step = round(24/num_long_steps, 2);
assert(mod(num_long_steps,24)==0, '长时间步数必须是24的整数倍');
repeat_factor = num_long_steps/24; % 每小时的步数

%% ===================== 创建结果存储结构 =====================
results = struct(...
    'P_agg',        zeros(1, total_steps), ...  % 聚合功率
    'P_base',       zeros(1, total_steps), ...  % 基准功率
    'S_agg',        zeros(1, total_steps), ...  % 虚拟SOC
    'lambda',       zeros(1, total_steps), ...  % λ*
    'EV_S_original',zeros(length(EVs), total_steps)... % 存储所有EV的SOC
);

%% ===================== 初始化EV目标充电功率 =====================
for i = 1:length(EVs)
    Delta_E = calculateDeltaE(EVs(i), EVs(i).p_real);
    EVs(i).E_tar = max(EVs(i).E_tar_set - Delta_E, EVs(i).E_ini); % 电量下限保护
end

%% ===================== 初始化聚合SOC =====================
S_agg_current = mean([EVs.S_original]);
[S_agg_opt, P_base_opt] = EVbaseP_aggregate(EVs, S_agg_current, 24, dt_long);
for i = 1:length(EVs)
    [~, P_base_hourly] = EVbaseP_aggregate(EVs(i), 0.01, 24, dt_long); 
    % 时间轴对齐（关键修改）
    EVs(i).P_base = repelem(P_base_hourly, repeat_factor); 
end

%% ===================== 外层循环（长时间步长） =====================
for long_idx = 1:num_long_steps
    t_long = (long_idx-1)*dt_long; % 当前长步起始时间
    
    %% --------------------- 长时间步处理 ---------------------
    % 阶段1: 计算聚合系数
    [M1, M2, M3, P_agg_max] = calculateAggregationCoefficients(EVs, dt_long);
    
    % 阶段3: 执行需求响应
    [~, lambda_star, S_agg_next] = aggregateEVs(...
        EVs, P_tar(long_idx), M1, M2, M3, S_agg_current...
    );
    for i = 1:length(EVs)
        % 调用调节潜力计算函数（注意保持原始参数不变）
        [EVs(i).DeltaP_pu_plus_t, EVs(i).DeltaP_pu_minus_t] = calculateEVAdjustmentPotentia(...
            EVs(i).C, EVs(i).r, EVs(i).eta, EVs(i).E_tar, EVs(i).E_ini, EVs(i).E_actual, ...
            EVs(i).t_dep, EVs(i).t_in, EVs(i).P_N, 0, EVs(i).P_base(long_idx), ...
            EVs(i).S_original, dt_long); % 参数顺序保持原始结构
    end
    % S_agg_current_short = S_agg_current;
    %% ===================== 内层循环（短时间步长） =====================
    for short_idx = 1:num_short_per_long
        step_idx = (long_idx-1)*num_short_per_long + short_idx;
        t_current = t_long + (short_idx-1)*dt_short;
        
        %% --------------------- 更新EV状态 ---------------------
        for i = 1:length(EVs)
            EV = EVs(i);
            % 状态更新流程
            EV = updateLockState(EV, t_current); 
            EV = generateDemandCurve(EV);
            EV.P_current = EV.demandCurve(lambda_star);
            EV = calculateVirtualSOC(EV, t_current, dt_short);
            m3 = (EV.E_tar - EV.E_ini)/(EV.eta*((EV.t_dep - EV.t_in)/60));
            results.m3(i) = m3;
            EVs(i) = EV;
        end
        %% --------------------- 更新短步长聚合SOC ---------------------
        % [M1, M2, M3, ~] = calculateAggregationCoefficients(EVs, dt_short);
        % [~, ~, S_agg_next_short] = aggregateEVs(...
        % EVs, P_tar(long_idx), M1, M2, M3, S_agg_current_short);
        
        %% --------------------- 记录结果 ---------------------
        % results.P_agg(step_idx)   = P_agg;
        results.lambda(step_idx) = lambda_star;
        % results.S_agg(step_idx)  = S_agg_next_short;
        
        % 记录所有EV的SOC
        for ev_idx = 1:length(EVs)
            results.EV_S_original(ev_idx, step_idx) = EVs(ev_idx).S_original;
        end

        %% --------------------- 更新短步长聚合SOC ---------------------
        % S_agg_current_short = S_agg_next_short;
    end
    
    %% --------------------- 更新聚合SOC ---------------------
    S_agg_current = S_agg_next;
end

%% ===================== 可视化结果 =====================
% plotEnhancedResults(results, dt_short);
