load('result_1000_inc.mat');
[EVs, t_sim, dt_short, dt_long, P_tar] = initializeFromExcel(excelFile);
num_long_steps = t_sim / dt_long;
num_short_per_long = dt_long / dt_short;
total_steps = num_long_steps * num_short_per_long;
selected_ev = 10;
time_min = (0:total_steps-1) * dt_short;
time_h = time_min / 60;

S_original = results.EV_S_original(selected_ev, :);
S_mod = results.EV_S_mod(selected_ev, :);
P_current_plot = results.P_cu;

figure('Position', [100 100 1200 600]);

yyaxis left;
plot(time_h, S_original, 'b-', 'LineWidth', 1.5, 'DisplayName', '原始SOC (S_{original})');
hold on;
plot(time_h, S_mod, 'r--', 'LineWidth', 1.5, 'DisplayName', '修正SOC (S_{modified})');
ylabel('SOC');
ylim([-1.1 1.1]);
grid on;

yyaxis right;
plot(time_h, P_current_plot, 'g-.', 'LineWidth', 1.5, 'DisplayName', ['EV ', num2str(selected_ev), ' 功率 (P_{current})']);
ylabel('功率 (kW)');
P_lim = [min(P_current_plot)-1, max(P_current_plot)+1];
ylim(P_lim);

xlabel('时间 (小时)');
title(['第', num2str(selected_ev), '台EV的SOC与当前功率对比']);
legend('Location', 'best');
set(gca, 'FontSize', 12);
hold off;

% plotLambdaAndAggSOC(results, dt_short);
% plotPowerComparison(results, dt_short);