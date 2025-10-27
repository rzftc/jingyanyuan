function [EVs, t_sim, dt_short, dt_long, P_tar] = initializeFromExcel(filePath)
    % --- 参数校验 ---
    if ~exist(filePath, 'file')
        error('EV参数文件不存在: %s', filePath);
    end

    % --- 读取数据 ---
    data = readtable(filePath); % 读取原始数据表

    % --- 基础参数 ---
    t_sim = 48*60;      % 模拟总时长 (分钟) - 48小时
    dt_short = 5;     % 短时间步长 (分钟)
    dt_long = 30;       % 长时间步长 (分钟)

    % --- 创建EV结构体模板 ---
    template = struct(...
        'P_N', [], 'P_0', [], 'P_h_max', [], 'P_l_min', [],...
        'C', [], 'eta', [], 'r', [],...
        't_in', [], 't_dep', [],...
        'E_ini', [], 'E_tar_set', [], 'E_tar', [],...
        'state', [], 'demandCurve', [],...
        'S_original', 0, 'S_modified', 0,...
        'P_current', 0, 'substate', 'OFF','BaseSchedule', zeros(48*(60/dt_long),1),... % 修正BaseSchedule长度
        'Delta_E_h_max', [], 'Delta_E_q_max', [], 'p_real', [], 'tau_rem', 0, ...
        'E_exp', [], 'E_actual', [], 'EV_ID', [] ... % 补全字段以便下方引用
    );

    EVs = repmat(template, height(data),1);

    % --- 参数映射与校验 ---
    P_N_vec = zeros(height(data), 1); % 预分配 P_N 向量
    for i = 1:height(data)
        % 额定功率校验
        if data.P_N(i) <= 0
            error('EV%d: 无效额定功率%.2fkW', i, data.P_N(i));
        end
        % 电价关系校验
        if ~(data.P_l_min(i) < data.P_0(i) && data.P_0(i) < data.P_h_max(i))
            error('EV%d电价关系错误: %.2f < %.2f < %.2f 不成立',...
                  i, data.P_l_min(i), data.P_0(i), data.P_h_max(i));
        end

        % 参数映射
        EVs(i).EV_ID = data.EV_ID(i); % 映射 EV_ID
        EVs(i).P_N = data.P_N(i);
        EVs(i).P_0 = data.P_0(i);
        EVs(i).P_h_max = data.P_h_max(i);
        EVs(i).P_l_min = data.P_l_min(i);
        EVs(i).C = data.C(i);
        EVs(i).eta = data.eta(i);
        EVs(i).r = data.r(i);
        EVs(i).t_in = data.t_in(i);
        EVs(i).t_dep = data.t_dep(i);
        EVs(i).E_ini = data.E_ini(i);
        EVs(i).E_tar_set = data.E_tar_set(i);
        EVs(i).state = data.state{i}; % 假设 state 是 cellstr
        EVs(i).Delta_E_h_max = data.Delta_E_h_max(i);
        EVs(i).Delta_E_q_max = data.Delta_E_q_max(i);
        EVs(i).p_real = data.p_real(i);
        EVs(i).p_incentive = data.p_incentive(i);

        % 计算字段初始化
        EVs(i).E_tar = EVs(i).E_tar_set; % E_tar 稍后会被重新计算
        EVs(i).E_exp = EVs(i).E_ini;
        EVs(i).E_actual = EVs(i).E_ini;
        EVs(i).tau_rem = 0; % tau_rem 稍后会被重新计算

        P_N_vec(i) = EVs(i).P_N; % 存储 P_N 到向量
    end

    % --- 生成基于车辆数据的目标功率曲线 P_tar ---
    num_long_steps = t_sim / dt_long;
    P_max_potential = zeros(1, num_long_steps); % 存储每个长步的最大潜在功率
    t_in_vec = [EVs.t_in]';   % 向量化 t_in
    t_dep_vec = [EVs.t_dep]'; % 向量化 t_dep

    for j = 1:num_long_steps
        current_t_start = (j-1) * dt_long; % 当前长步开始时间
        current_t_end = j * dt_long;       % 当前长步结束时间 (用于更精确判断连接状态)

        % 找到在此长步时间段内 '在线' 的 EV
        % 条件：入网时间 <= 长步结束时间 AND 离网时间 > 长步开始时间
        connected_mask = (t_in_vec < current_t_end) & (t_dep_vec > current_t_start);

        % 计算这些在线 EV 的额定功率总和
        if any(connected_mask)
            P_max_potential(j) = sum(P_N_vec(connected_mask));
        else
            P_max_potential(j) = 0;
        end
    end

    % 设置 P_tar 为最大潜在功率的一个比例 (例如 60%)
    load_factor = 0.6; % 可调参数：目标负荷率
    P_tar = P_max_potential * load_factor;

end