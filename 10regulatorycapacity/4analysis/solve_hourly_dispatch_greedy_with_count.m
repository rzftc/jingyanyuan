%% 局部函数: 下层贪心算法 (带数量约束) - solve_hourly_dispatch_greedy_with_count (同前)
function [u_ac_optimal, u_ev_optimal, total_cost] = solve_hourly_dispatch_greedy_with_count(p_ac, p_ev, c_ac, c_ev, P_req_t, n_ac_max, n_ev_max)
    % ... (函数体与上一版本完全相同) ...
    u_ac_optimal = zeros(length(p_ac), 1); u_ev_optimal = zeros(length(p_ev), 1);
    total_cost = 0;
    if P_req_t <= 0; return; end

    num_ac_avail = sum(p_ac > 1e-6); num_ev_avail = sum(p_ev > 1e-6);
    num_devices = num_ac_avail + num_ev_avail;
    if num_devices == 0; return; end

    device_data = struct('id', cell(num_devices, 1), 'type', cell(num_devices, 1), ...
                         'power', zeros(num_devices, 1), 'cost_per_kw', zeros(num_devices, 1), ...
                         'total_cost', zeros(num_devices, 1));
    idx = 1;
    for i = 1:length(p_ac); if p_ac(i) > 1e-6; device_data(idx).id=i; device_data(idx).type='AC'; device_data(idx).power=p_ac(i); device_data(idx).cost_per_kw=c_ac(i); device_data(idx).total_cost=p_ac(i)*c_ac(i); idx=idx+1; end; end
    for i = 1:length(p_ev); if p_ev(i) > 1e-6; device_data(idx).id=i; device_data(idx).type='EV'; device_data(idx).power=p_ev(i); device_data(idx).cost_per_kw=c_ev(i); device_data(idx).total_cost=p_ev(i)*c_ev(i); idx=idx+1; end; end

    [~, sorted_indices] = sort([device_data.cost_per_kw]);
    sorted_devices = device_data(sorted_indices);

    power_accumulated = 0; count_ac = 0; count_ev = 0;
    for i = 1:length(sorted_devices)
        if power_accumulated >= P_req_t; break; end
        device = sorted_devices(i);
        if strcmp(device.type, 'AC')
            if count_ac >= n_ac_max; continue; end
            count_ac = count_ac + 1; u_ac_optimal(device.id) = 1;
        else % EV
            if count_ev >= n_ev_max; continue; end
            count_ev = count_ev + 1; u_ev_optimal(device.id) = 1;
        end
        power_accumulated = power_accumulated + device.power;
        total_cost = total_cost + device.total_cost;
    end
    % 可以在这里加一个警告，如果最终 power_accumulated < P_req_t
     if power_accumulated < P_req_t - 1e-3 % 允许一点误差
          %fprintf('警告: 时间步 %d 贪心未能满足需求 %.2f kW (调度 %.2f kW) n_ac=%d/%d, n_ev=%d/%d\n', t_global, P_req_t, power_accumulated, count_ac, n_ac_max, count_ev, n_ev_max);
          % 在 parfor 中无法直接访问 t_global，可以在主循环结束后统一检查
     end
end

