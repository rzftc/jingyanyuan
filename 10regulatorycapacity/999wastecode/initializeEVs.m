function EVs = initializeEVs(num_EV, rngSeed)
    % 参数校验
    if nargin < 2, rngSeed = 2023; end
    if num_EV <= 0, error('电动汽车数量必须为正整数'); end
    
    % 初始化随机数生成器（确保结果可复现）
    rng(rngSeed);
    
    % 参数范围定义
    paramRanges = struct(...
        'C_EV',   [40, 100],...    % 电池容量 (kWh)
        'eta',    [0.8, 0.95],...  % 充放电效率
        'p_on',   [3, 5],...       % 接入功率范围 (kW)
        'E_in',   [10, 40],...     % 初始电量范围 (kWh)
        't_in',   [6, 10],...      % 接入时间范围 (小时)
        't_dep',  [2, 8],...       % 停留时长范围 (小时)
        'p_incentive', [0, 60]...  % 激励电价范围
        );
    
    % 预分配结构体数组
    EVs(num_EV) = struct();
    
    % 生成EV参数（包含新增的ptcp字段）
    for i = 1:num_EV
        % 基础参数
        EVs(i).C_EV = paramRanges.C_EV(1) + diff(paramRanges.C_EV)*rand();
        EVs(i).eta = paramRanges.eta(1) + diff(paramRanges.eta)*rand();
        EVs(i).p_on = paramRanges.p_on(1) + diff(paramRanges.p_on)*rand();
        
        % 时间参数（带正态分布扰动）
        base_time = fix(paramRanges.t_in(1) + diff(paramRanges.t_in)*rand());
        EVs(i).t_in = base_time + randn()*0.5;  % 带标准差0.5的正态分布
        EVs(i).t_dep = EVs(i).t_in + paramRanges.t_dep(1) + diff(paramRanges.t_dep)*rand();
        EVs(i).p_incentive = paramRanges.p_incentive(1) + diff(paramRanges.p_incentive)*rand();
        
        % 能量参数
        EVs(i).E_in = paramRanges.E_in(1) + diff(paramRanges.E_in)*rand();
        EVs(i).E_tar = EVs(i).E_in + (EVs(i).C_EV - EVs(i).E_in)*rand();
        
        % 状态参数
        EVs(i).r = 0.025;              % 默认充电速率
        EVs(i).E_exp = EVs(i).E_in;     % 预期电量
        EVs(i).P_current = 0;          % 当前功率
        EVs(i).SOC = 0;                % 初始SOC
        EVs(i).P_base = 0;             % 基础功率
        EVs(i).E_current = EVs(i).E_in; % 当前电量
        EVs(i).ptcp = (rand() < 0.5);  % 新增参与标志（默认50%概率）
    end
end
