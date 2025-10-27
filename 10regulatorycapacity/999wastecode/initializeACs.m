function ACs = initializeACs(num_AC, rngSeed)
    % 参数校验
    if nargin < 2, rngSeed = 2023; end
    if num_AC <= 0, error('空调数量必须为正整数'); end
    
    % 初始化随机数生成器
    rng(rngSeed);
    
    % 参数范围定义（带续行符的标准格式）
    paramRanges = struct(...
        'R',    [1.5, 3.0],...    % 热阻范围 (kΩ)
        'C',    [0.3, 0.7],...    % 电容范围 (kF)
        'eta',  [0.85, 0.95],...  % 效率范围
        'Tset', [19, 25],...      % 初始温度设定范围 (℃)
        'SOC',  [0, 0.5],...       % 初始SOC范围
        'p_incentive', [0, 60]... % 激励电价范围[0,60]
        );
    
    % 预分配结构体数组
    ACs(num_AC) = struct();
    
    % 生成空调参数
    for i = 1:num_AC
        ACs(i).R = paramRanges.R(1) + diff(paramRanges.R)*rand();  % 热阻随机值
        ACs(i).C = paramRanges.C(1) + diff(paramRanges.C)*rand();  % 电容随机值
        ACs(i).eta = paramRanges.eta(1) + diff(paramRanges.eta)*rand();
        ACs(i).Tset = paramRanges.Tset(1) + diff(paramRanges.Tset)*rand();
        ACs(i).Tmax = 26;  % 温度上限
        ACs(i).Tmin = 18;  % 温度下限
        ACs(i).SOC = paramRanges.SOC(1) + diff(paramRanges.SOC)*rand();
        ACs(i).P_base = 0;  % 基础功率
        ACs(i).p_incentive = paramRanges.p_incentive(1) + ...
                            diff(paramRanges.p_incentive)*rand();
    end
end
