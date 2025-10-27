clc; clear; close all;
rng(2024);

%% ===================== 初始化参数 =====================
excelFile = 'resi_1000.xlsx';
if ~exist(excelFile, 'file')
    generateEVParameters_real(excelFile, 100, 0.6);
    fprintf('已生成参数模板: %s\n', excelFile);
end  
[EVs, t_sim, dt_short, dt_long, P_tar] = initializeFromExcel(excelFile);
fprintf('成功加载%d辆EV数据\n', length(EVs));

%% ===================== 仿真时间调整 =====================
t_start_sim = 6 * 60; % 仿真从早上6点开始 (360分钟)

%% ===================== 索引预处理 =====================
num_long_steps = t_sim / dt_long;
num_short_per_long = dt_long / dt_short;
total_steps = num_long_steps * num_short_per_long;
assert(mod(num_long_steps, 24) == 0, '长时间步数必须是24的整数倍');

%% ===================== 创建结果存储结构 =====================
num_evs = length(EVs);
results = struct(...
    'P_agg',         zeros(1, total_steps), ...
    'P_base',        zeros(1, total_steps), ...
    'S_agg',         zeros(1, total_steps), ...
    'lambda',        zeros(1, total_steps), ...
    'EV_S_original', zeros(num_evs, total_steps), ...
    'EV_S_mod',      zeros(num_evs, total_steps), ...
    'm3',            zeros(num_evs, 1), ...
    'P_tar',         zeros(1, total_steps), ...
    'P_cu',          zeros(1, total_steps) ...
    );

%% ===================== 初始化EV目标充电功率 (向量化) =====================
p_real_vec = [EVs.p_real]';
P_h_max_vec = [EVs.P_h_max]';
P_0_vec = [EVs.P_0]'; % P_0_vec 仍然被定义，但在此处不再用于 calculateParticipation
P_l_min_vec = [EVs.P_l_min]';
Delta_E_h_max_vec = [EVs.Delta_E_h_max]';
Delta_E_q_max_vec = [EVs.Delta_E_q_max]';
E_tar_set_vec = [EVs.E_tar_set]';
E_ini_vec = [EVs.E_ini]';
P_N_vec = [EVs.P_N]';

% 定义代码片段中的常量
p_min = 15;       % 原始激励范围下限
p_max = 50;       % 原始激励范围上限
p_min_prime = 10; % 调整后激励范围下限
p_max_prime = 40; % 调整后激励范围上限

% 获取车辆容量 (假设字段名为 'C')
if ~isfield(EVs, 'C')
    error('EVs结构体数组必须包含表示电池容量的字段 ''C''');
end
C_vec = [EVs.C]';
E_tar_max_val_vec = 0.2 * C_vec; % incentiveTempEV 所需的 E_tar_max

% 获取激励价格 (假设字段名为 'p_incentive')
if ~isfield(EVs, 'p_incentive')
    error('EVs结构体数组必须包含表示激励价格的字段 ''p_incentive''');
end
p_incentive_vec = [EVs.p_incentive]';

% 定义固定的基准电价向量 base_P
base_Price = 30 * ones(num_evs, 1);

if ~exist('calculateParticipation', 'file')
    error(['函数 calculateParticipation(incentive_price, base_price) 未找到。',...
        '请定义该函数。它应返回参与概率。']);
end
% 使用新定义的 base_P 调用 calculateParticipation
participation_probabilities_vec = calculateParticipation(p_incentive_vec, base_Price);
if any(participation_probabilities_vec < 0) || any(participation_probabilities_vec > 1)
    warning('calculateParticipation 函数返回的参与概率超出了[0,1]范围。');
end

% 根据概率确定实际参与情况
ptcp_vec = (rand(num_evs, 1) < participation_probabilities_vec);

% 和 E_tar_max (最大电量变化)，并返回向量 deltaE。
if ~exist('incentiveTempEV', 'file')
    error(['函数 incentiveTempEV(p, p_min, p_max, p_min_prime, p_max_prime, E_tar_max) 未找到。',...
        '请定义该函数。它应返回三个向量, 第三个是deltaE。']);
end

% 直接使用向量调用 incentiveTempEV 函数
[~, ~, deltaE_vec] = incentiveTempEV(p_incentive_vec, p_min, p_max, p_min_prime, p_max_prime, E_tar_max_val_vec);

% 使用 E_tar_set_vec (激励前的目标) 初始化
E_tar_vec_intermediate = E_tar_set_vec;

% 对参与的EV应用 deltaE
participating_indices = find(ptcp_vec);

if ~isempty(participating_indices)
    potential_E_tar_for_participants = E_tar_set_vec(participating_indices) - deltaE_vec(participating_indices);
    revert_condition_met = potential_E_tar_for_participants <= E_ini_vec(participating_indices);
    apply_reduction_indices = participating_indices(~revert_condition_met);
    E_tar_vec_intermediate(apply_reduction_indices) = potential_E_tar_for_participants(~revert_condition_met);
    % 对于 revert_condition_met 为 true 的参与者, E_tar_vec_intermediate 保持为 E_tar_set_vec (原始值),
    % 这实现了逻辑: E_tar = E_tar - deltaE; 如果 E_tar <= E_ini, 则 E_tar = E_tar + deltaE (回到原始的 E_tar_set)。
end

% 最终的 E_tar 必须至少为 E_ini_vec
E_tar_vec = max(E_tar_vec_intermediate, E_ini_vec);
fprintf('新的E_tar生成逻辑已成功应用。\n');
% --- MODIFIED E_TAR GENERATION LOGIC ENDS ---

t_ch_vec = zeros(num_evs, 1);
valid_PN_mask = P_N_vec > 0;
if any(valid_PN_mask)
    t_ch_vec(valid_PN_mask) = 60 .* (E_tar_vec(valid_PN_mask) - E_ini_vec(valid_PN_mask)) ./ P_N_vec(valid_PN_mask);
end
t_ch_vec = max(t_ch_vec, 0);

E_tar_cell = num2cell(E_tar_vec);
[EVs.E_tar] = E_tar_cell{:};
tau_rem_cell = num2cell(t_ch_vec);
[EVs.tau_rem] = tau_rem_cell{:};
E_ini_cell = num2cell(E_ini_vec);
[EVs.E_exp] = E_ini_cell{:};
[EVs.E_actual] = E_ini_cell{:};

%% ===================== 初始化聚合SOC =====================
S_agg_current = mean([EVs.S_original]);
% [~, P_base_opt] = EVbaseP(EVs, 24, dt_long);
for i=1:num_evs
    [P_base_sequence, ~] = EVbaseP(EVs(i), 24*60, dt_short);
    EVs(i).P_base = P_base_sequence;
end
%% ===================== 外层循环（长时间步长） =====================
for long_idx = 1:num_long_steps
    t_long = (long_idx - 1) * dt_long + t_start_sim;
    
    %% --------------------- 长时间步处理 ---------------------
    % EVs = distributeBasePower(EVs, 1, P_base_opt);
    [lambda_star] = aggregateEVs(EVs, P_tar(long_idx));
    [~, S_agg_next] = calculateVirtualSOC_agg(EVs, dt_long);
    
    %% ===================== 内层循环（短时间步长） =====================
    for short_idx = 1:num_short_per_long
        step_idx = (long_idx - 1) * num_short_per_long + short_idx;
        t_current = t_long + (short_idx - 1) * dt_short;
        
        temp_m3 = zeros(num_evs, 1);
        temp_S_original = zeros(num_evs, 1);
        temp_S_mod = zeros(num_evs, 1);
        temp_P_current = zeros(num_evs, 1);
        
        %% --------------------- 更新EV状态 (并行处理) ---------------------
        parfor i = 1:num_evs % 使用 parfor 加速
            EV = EVs(i);
            % 注意：t_current 现在可能大于1440，需要取模以匹配t_in/t_dep
            EV = updateLockState(EV, mod(t_current, 1440));
            
            EV_for_handle = EV;
            EV_temp_with_handle = generateDemandCurve(EV_for_handle);
            
            current_P_val = 0;
            if isfield(EV_temp_with_handle, 'demandCurve') && isa(EV_temp_with_handle.demandCurve, 'function_handle')
                current_P_val = EV_temp_with_handle.demandCurve(lambda_star);
            else
                switch EV.state
                    case 'LockON'
                        current_P_val = EV.P_N;
                    case {'LockOFF', 'OFF'}
                        current_P_val = 0;
                    otherwise
                        current_P_val = 0;
                end
            end
            EV.P_current = current_P_val;
            
            EV = calculateVirtualSOC_upgrade(EV, t_current, dt_short);
            
            if EV.t_dep > EV.t_in
                m3_val = (EV.E_tar - EV.E_ini) / (EV.eta * ((EV.t_dep - EV.t_in) / 60));
            else
                m3_val = 0;
            end
            
            temp_m3(i) = m3_val;
            temp_S_original(i) = EV.S_original;
            temp_S_mod(i) = EV.S_modified;
            temp_P_current(i) = EV.P_current;
            
            EVs(i) = EV;
            
        end % 结束 parfor i
        
        results.EV_S_original(:, step_idx) = temp_S_original;
        results.EV_S_mod(:, step_idx) = temp_S_mod;
        results.m3 = temp_m3;
        
        results.lambda(step_idx) = lambda_star;
        results.S_agg(step_idx) = S_agg_current;
        results.P_agg(step_idx) = sum(temp_P_current);
        
        results.P_cu(step_idx) = temp_P_current(10);
        
    end % 结束 short_idx
    
    %% --------------------- 更新聚合SOC ---------------------
    S_agg_current = S_agg_next;
    
end % 结束 long_idx

results.P_tar = repelem(P_tar, num_short_per_long);

