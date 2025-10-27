%% 修正后的初始化函数（添加E_tar字段）
function [EVs, t_sim, dt_short, dt_long, P_tar] = initializeParameters()
    num_EV    = 100;
    t_sim     = 24 * 60;
    dt        = 1;
    time_h    = (0 : dt : t_sim - 1) / 60;  % 时间轴（小时）
    dt_short  = 1;                           % 短时间步长(分钟)-状态更新
    dt_long   = 10;                          % 长时间步长(分钟)-需求曲线更新 
    assert(mod(dt_long, dt_short) == 0, '时间步长需整数倍关系');
    rng(2023);
    
    % 预定义所有字段模板（新增E_tar）
    template = struct(...
        'C', [], 'eta', [], 'r', [],...
        'P_N', [], 'P_h_max', [], 'P_l_min', [], 'P_0', [],...
        'Delta_E_h_max', [], 'Delta_E_q_max', [],...
        't_in', [], 't_dep', [],...
        'E_ini', [], 'E_tar_set', [], 'E_tar', [],...
        'E_exp', [], 'E_actual', [],...
        'S', [], 'state', [], 'demandCurve', [],...
        'S_original', 0, 'S_modified', 0,...
        'P_current', 0, 'substate', 'OFF');  % 新增substate字段
    
    % 构造P_tar曲线（可根据实际需求替换为外部数据）
    time_h_long = (0 : dt_long : t_sim - 1) / 60;
    P_tar = 1000 * sin(2 * pi * time_h_long / 24) + 500;  % 正弦波峰谷负荷
    
    EVs = repmat(template, 1, num_EV);
    
    for i = 1 : num_EV
        % 参数赋值
        EVs(i).C             = 50 + 20 * rand();
        EVs(i).eta           = 0.9;
        EVs(i).r             = 0.05;
        
        EVs(i).P_N           = 7 + 3 * rand();
        EVs(i).P_h_max       = 1.2 * EVs(i).P_N;
        EVs(i).P_l_min       = 0.5 * EVs(i).P_N;
        EVs(i).P_0           = 0.8 * EVs(i).P_N;
        EVs(i).Delta_E_h_max = 0.1 * EVs(i).C;
        EVs(i).Delta_E_q_max = 0.05 * EVs(i).C;
        EVs(i).p_real = 6;
        
        EVs(i).t_in          = randi([1, 6 * 60]);
        EVs(i).t_dep         = EVs(i).t_in + 8 * 60 + randi([-60, 60]);
        
        EVs(i).E_ini         = 0.2 * EVs(i).C;
        EVs(i).E_tar_set     = 0.8 * EVs(i).C;
        EVs(i).E_tar         = EVs(i).E_tar_set;  % 初始化E_tar为设定值
        EVs(i).E_exp         = EVs(i).E_ini;
        EVs(i).E_actual      = EVs(i).E_ini;
        
        EVs(i).S             = 0;
        EVs(i).state         = 'OFF';
        EVs(i).demandCurve   = @(lambda) 0;
        EVs(i).S_original    = 0;                % 原始SOC
        EVs(i).S_modified    = 0;                % 修正SOC
    end
end
