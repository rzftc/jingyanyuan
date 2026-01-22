function [HE, ME] = ev_charging_infra_prediction(CE, CE_p, lambda1, lambda2)
    % EV充电基础设施规模预测函数（高低情景）
    % 输入参数：
    %   CE: 三维数组 [num_types, D, N]，表示i类充电桩在j天k时段的需求负荷（单位：kW）
    %       - num_types: 充电桩类别数（如5类）
    %       - D: 计算天数（默认7天）
    %       - N: 每天时间间隔数（默认144，对应10分钟分辨率）
    %   CE_p: 列向量 [num_types×1]，表示i类充电桩的总峰值负荷（单位：kW，取所有天时段的最大值）
    %   lambda1: 日均平均利用率（默认0.2，即20%）
    %   lambda2: 峰值利用率参数（默认1.25，每桩对应1.25辆充电需求）
    % 输出：
    %   HE: 高情景充电桩规模（数量，单位：桩） [num_types×1]
    %   ME: 低情景充电桩规模（数量，单位：桩） [num_types×1]
    
    %% 参数校验与默认值设置
    [num_types, D_input, N_input] = size(CE);
    assert(num_types >= 1, '充电桩类别数至少为1');
    assert(D_input == 7, '计算天数D应等于7天（与理论模型一致）');
    assert(N_input == 144, '每天时间间隔数N应等于144（10分钟分辨率，24×6）');
    assert(size(CE_p,1) == num_types && size(CE_p,2) == 1, 'CE_p需为[num_types×1]的列向量');
    
    % 设置默认参数（若未输入）
    if nargin < 3 || isempty(lambda1), lambda1 = 0.2; end    % 日均利用率20%
    if nargin < 4 || isempty(lambda2), lambda2 = 1.25; end   % 峰值参数1.25
    
    %% 高情景规模计算（HE_i）：基于日均总需求和平均利用率
    % 公式：HE_i = (λ₁ × Σ（j=1到D）Σ（k=1到N）CE_(i,j,k)) / (N×D)
    HE = zeros(num_types, 1);
    for i = 1:num_types
        % 计算i类充电桩7天内所有时段的总需求（kW·时段）
        total_demand = sum(sum(CE(i, :, :)));  % 维度：[1×1]（对j和k求和）
        % 日均总需求 = 总需求 / (N×D)（单位：kW）
        avg_daily_demand = total_demand / (N_input * D_input);
        % 高情景规模：考虑平均利用率λ₁，即充电桩需满足日均需求/利用率
        HE(i) = avg_daily_demand / lambda1;
    end
    
    %% 低情景规模计算（ME_i）：基于峰值需求和峰值利用率
    % 公式：ME_i = CE_i^p / λ₂
    ME = CE_p / lambda2;
    
    %% 结果格式化（保留2位小数，符合工程精度）
    HE = round(HE, 2);
    ME = round(ME, 2);
    
    %% 输出提示
    disp('===== 充电桩高低情景规模预测结果 =====');
    disp('充电桩类别 | 高情景规模（桩） | 低情景规模（桩）');
    for i = 1:num_types
        fprintf('%d            | %.2f             | %.2f\n', i, HE(i), ME(i));
    end
end
