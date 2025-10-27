clc; 
close all; 
clear;

%% ===================== 初始化参数 =====================
[EVs, t_sim, dt_short, dt_long, P_tar] = initializeParameters();
H = 12;   % 预测时域（按长时间步长计算）

%% ===================== 创建结果存储结构 =====================
num_steps    = t_sim;  % 总时间步数（按短步长计数）
target_EV_idx = 2;     % 目标EV索引
results      = struct(...
    'P_agg',        zeros(1, num_steps), ...  % 聚合功率（短步长记录）
    'P_base',       zeros(1, num_steps), ...  % 基准功率
    'S_agg',        zeros(1, num_steps), ...  % 虚拟SOC
    'lambda',       zeros(1, num_steps), ...  % λ*
    'EV_S_original',zeros(1, num_steps) ...   % 目标EV原始SOC
);
substate_log = cell(ceil(t_sim/dt_short), 1); 
%% ===================== 初始化EV目标充电功率 =====================
for i = 1:length(EVs)
    Delta_E = calculateDeltaE(EVs(i), EVs(i).p_real);
    EVs(i).E_tar = EVs(i).E_tar - Delta_E;
    if EVs(i).E_tar < EVs(i).E_ini
        EVs(i).E_tar = EVs(i).E_tar + Delta_E;
    end
end
%% ===================== 初始化聚合SOC =====================
S_agg_current = mean([EVs.S_original]);

%% ===================== 外层循环（长时间步长） =====================
for t_long = 0 : dt_long : t_sim - 1
    % 当前长时间步索引
    idx_long = t_long / dt_long + 1;
    
    %% --------------------- 长时间步长处理 ---------------------
    % 获取当前电网需求功率
    current_P_tar = P_tar(idx_long);
    
    % 计算聚合系数
    [M1, M2, M3] = calculateAggregationCoefficients(EVs, dt_long);
    
    % 计算聚合功率和λ*
    [P_agg, lambda_star, S_agg_next] = aggregateEVs(...
        EVs, current_P_tar, M1, M2, M3, S_agg_current...
    );

    %% ===================== 内层循环（短时间步长） =====================
    for t_short = t_long : dt_short : t_long + dt_long - dt_short
        idx_short = t_short + 1;  % 短步长索引
        
        %% --------------------- 更新EV状态 ---------------------
        for i = 1:length(EVs)
            EV = EVs(i);
            % 状态更新流程
            EV = updateLockState(EV, t_short);        % 更新闭锁状态
            EV = generateDemandCurve(EV);             % 生成需求曲线
            EV.P_current = EV.demandCurve(lambda_star); % 计算当前功率
            EV = updateEVState(EV, t_short, dt_short);  % 更新SOC
            EVs(i) = EV;
        end
        
        %% --------------------- 记录结果 ---------------------
        step_idx = t_short/dt_short + 1;
        results.P_agg(idx_short)         = P_agg;
        results.lambda(idx_short)        = lambda_star;
        results.EV_S_original(idx_short) = EVs(target_EV_idx).S_original;
        substate_log{step_idx}           = EVs(target_EV_idx).substate;
    end
    
    %% --------------------- 更新聚合SOC ---------------------
    S_agg_current = S_agg_next;
end

%% ===================== 可视化结果 =====================
plotResults(results, dt_short);

 
