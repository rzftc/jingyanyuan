%% generateExampleExcel.m
function generateExampleExcel(num_devices_ac, num_devices_ev)
    % 参数初始化
    if nargin < 2
        num_devices_ac = 200;
        num_devices_ev = 100; 
    end

    % ===== 空调数据生成 =====
    acParams = struct(...
        'Tmin', 18,...
        'Tmax', 26,...
        'R_range', [1.5, 3.0],...
        'C_range', [0.3, 0.7]...
    );
    
    acData = struct(...
        'ID', num2cell(1:num_devices_ac)',...
        'R', num2cell(round(acParams.R_range(1) + diff(acParams.R_range)*rand(num_devices_ac,1), 2)),...
        'C', num2cell(round(acParams.C_range(1) + diff(acParams.C_range)*rand(num_devices_ac,1), 2)),...
        'eta', num2cell(round(0.85 + 0.1*rand(num_devices_ac,1), 2)),...
        'Tset', num2cell(round(acParams.Tmin + (acParams.Tmax-acParams.Tmin)*rand(num_devices_ac,1), 1)),...
        'Tmax', num2cell(repmat(acParams.Tmax, num_devices_ac,1)),...
        'Tmin', num2cell(repmat(acParams.Tmin, num_devices_ac,1)),...
        'SOC', num2cell(round(0.5*rand(num_devices_ac,1), 2)),...
        'p_incentive', num2cell(round(60*rand(num_devices_ac,1), 1))...
    );

    % ===== 电动汽车数据生成 =====
    evParams = struct(...
        'C_EV_range', [50, 110],...    % 电池容量(kWh)
        'eta_range', [0.8, 0.95],...   % 充电效率
        'p_on_range', [3, 5],...       % 充电功率(kW)
        't_in_range', [6, 10],...      % 接入时间(小时)
        'min_duration', 2,...          % 最小充电时长(小时)
        'max_duration', 8 );          % 最大充电时长(小时)
    
    % 预分配数组
    [C_EV, eta, p_on, t_in, t_dep, E_in, E_tar] = deal(zeros(num_devices_ev,1));
    
    for i = 1:num_devices_ev
        % 阶段1：生成独立参数
        C_EV(i) = round(evParams.C_EV_range(1) + diff(evParams.C_EV_range)*rand(), 1);
        eta(i) = round(evParams.eta_range(1) + diff(evParams.eta_range)*rand(), 2);
        t_in(i) = round(evParams.t_in_range(1) + diff(evParams.t_in_range)*rand(), 1);
        
        % 阶段2：动态生成充电参数
        while true
            % 生成充电功率和持续时间
            p_on(i) = round(evParams.p_on_range(1) + diff(evParams.p_on_range)*rand(), 1);
            duration = round(evParams.min_duration + (evParams.max_duration - evParams.min_duration)*rand(), 1);
            t_dep(i) = t_in(i) + duration;
            
            % 计算最大可充电量
            max_charge = p_on(i) * eta(i) * duration;
            
            % 生成初始电量（保留至少10%充电空间）
            E_in(i) = round(0.1*C_EV(i) + 0.3*C_EV(i)*rand(), 1);  % 10%-40%容量
            E_tar_max = min(E_in(i) + max_charge, 0.95*C_EV(i));  % 保留5%缓冲
            
            % 生成目标电量
            E_tar(i) = E_in(i) + round(0.1*(E_tar_max - E_in(i)) + 0.9*(E_tar_max - E_in(i))*rand(), 1);
            
            % 计算m3并验证约束
            m3 = (E_tar(i) - E_in(i)) / (eta(i) * duration);
            if p_on(i) > m3 && E_in(i) < E_tar(i) && E_tar(i) < C_EV(i)
                break;  % 满足所有条件退出循环
            end
        end
    end

    % 构建EV数据结构
    evData = struct(...
        'ID', num2cell(1:num_devices_ev)',...
        'C_EV', num2cell(C_EV),...
        'eta', num2cell(eta),...
        'p_on', num2cell(p_on),...
        'E_in', num2cell(E_in),...
        'E_tar', num2cell(E_tar),...
        't_in', num2cell(t_in),...
        't_dep', num2cell(t_dep),...
        'r', num2cell(repmat(0.025, num_devices_ev,1)),...
        'SOC', num2cell(round(rand(num_devices_ev,1), 2)),...
        'p_incentive', num2cell(round(60*rand(num_devices_ev,1), 1))...
    );
     % 获取当前函数所在路径的上一级目录（项目根目录）
    currentDir = fileparts(mfilename('fullpath'));
    projectRoot = fullfile(currentDir, '..');
    % 构建目标保存路径
    saveDir = fullfile(projectRoot, '0input_data');
    % 确保目录存在
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    ac_file = fullfile(saveDir, 'AC_template.xlsx');
    ev_file = fullfile(saveDir, 'EV_template.xlsx');
    % ===== 写入文件 =====
    writetable(struct2table(acData), ac_file);
    writetable(struct2table(evData), ev_file);
    disp('数据文件已生成: AC_template.xlsx 和 EV_template.xlsx');
end
