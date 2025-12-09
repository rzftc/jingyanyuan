function P_signal_t = generate_hourly_regulation_signal(T_steps_total, steps_per_hour, num_hours, num_ac_total)
    % generate_hourly_regulation_signal: 生成小时级变化的电网调节指令
    %
    % 输入:
    %   T_steps_total  - 仿真总时间步数 (例如 481)
    %   steps_per_hour - 每小时的时间步数 (例如 20)
    %   num_hours      - 总小时数 (例如 24)
    %   num_ac_total   - 当前分块中的AC总数 (用于缩放指令幅度)
    %
    % 输出:
    %   P_signal_t     - (T x 1) 向量，电网总调节指令 (kW)
    %                    (正值 = VPP增加负荷/下调; 负值 = VPP降低负荷/上调)

    fprintf('  正在为 %d 台 AC 生成 %d 小时的调节指令...\n', num_ac_total, num_hours);
    P_signal_t = zeros(T_steps_total, 1);
    
    % 启发式估算：假设每台AC平均可提供 0.5kW 的调节能力
    % (这是一个非常简化的假设，用于保证指令幅度合理)
    P_potential_chunk_mag = num_ac_total * 0.3; 
    
    % 假设电网指令在总潜力的 +/- 40% 之间波动
    max_signal = P_potential_chunk_mag * 0.4;
    min_signal = -P_potential_chunk_mag * 0.4;

    for h = 1:num_hours
        % 1. 为当前小时生成 *一个* 固定的指令值
        % ( rand() * (max-min) + min )
        P_command_hourly = min_signal + (max_signal - min_signal) * rand();
        
        % 2. 确定该小时对应的所有时间步索引
        start_step = (h-1) * steps_per_hour + 1;
        end_step = min(h * steps_per_hour, T_steps_total); % 确保不超过总时长
        
        % 3. 将该指令值赋给这个小时内的所有时间步
        P_signal_t(start_step:end_step) = P_command_hourly;
    end
    
    fprintf('  指令幅度范围: %.2f kW (上调) 到 %.2f kW (下调)\n', min(P_signal_t), max(P_signal_t));
end