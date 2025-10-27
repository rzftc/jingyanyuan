clc;clear;close;
rng(2023, 'Threefry');                         
T_total = 24;                                   
dt = 0.05;                                      
time_points = 0:dt:T_total;                     
base_price = 30;            
p_incentive_range = linspace(0, 50, 25);        
p_min = 15; p_max = 50;                         
p_min_prime = 10; p_max_prime = 40;              

load('results_3d.mat')
cluster_Up = results_3D.EV_Up;      % 上调潜力 [时间点数×电价数]
cluster_Down = results_3D.EV_Down;  % 下调潜力 [时间点数×电价数]

%% 7. 三维上调潜力曲面（保存为PNG）
fig = figure('Position', [200 200 800 600]);
[X, Y] = meshgrid(time_points, p_incentive_range);
surf(X, Y, results_3D.EV_Up', 'EdgeColor','none')
% title('上调潜力三维分布')  % 移除标题
xlabel('时间 (小时)', 'FontSize', 16)  % 增大标签字体
ylabel('激励电价 (元)', 'FontSize', 16)
zlabel('调节潜力 (kW)', 'FontSize', 16)
colormap(jet)
cb = colorbar;  % 获取colorbar对象
cb.FontSize = 16;  % 增大colorbar字体
ax = gca;  % 获取坐标轴对象
ax.FontSize = 16;  % 增大刻度标签字体
view(30, 60)
print(fig, '3D_Up_Potential.png', '-dpng', '-r400');  % 400dpi保存

%% 8. 三维下调潜力曲面（保存为PNG）
fig = figure('Position', [200 200 800 600]);
surf(X, Y, results_3D.EV_Down', 'EdgeColor','none')
% title('下调潜力三维分布')  % 移除标题
xlabel('时间 (小时)', 'FontSize', 16)
ylabel('激励电价 (元)', 'FontSize', 16)
zlabel('调节潜力 (kW)', 'FontSize', 16)
colormap(jet)
cb = colorbar;
cb.FontSize = 16;
ax = gca;
ax.FontSize = 16;
view(30, 60)
print(fig, '3D_Down_Potential.png', '-dpng', '-r400');  % 400dpi保存

%% 9. 集群上调潜力随时间变化曲线（保存为PNG）
selected_p_idx = [3, 13, 23];  
selected_p_values = p_incentive_range(selected_p_idx);

fig = figure('Position', [200 200 1000 600]);
hold on
for i = 1:length(selected_p_idx)
    p_idx = selected_p_idx(i);
    plot(time_points, cluster_Up(:, p_idx), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('激励电价=%.1f元', selected_p_values(i)))
end
hold off
% title('集群上调潜力随时间变化曲线（不同激励电价）')  % 移除标题
xlabel('时间 (小时)', 'FontSize', 16)
ylabel('上调潜力 (kW)', 'FontSize', 16)
grid on
lgd = legend('Location', 'best');  % 获取图例对象
lgd.FontSize = 16;  % 增大图例字体
ax = gca;
ax.FontSize = 16;  % 增大刻度标签字体
up_max = max(cluster_Up(:));
ylim([0, up_max*1.1])
print(fig, 'Cluster_Up_TimeCurve.png', '-dpng', '-r400');  % 400dpi保存

%% 10. 集群上调潜力随激励电价变化曲线（保存为PNG）
selected_t = [8, 13, 18, 22];  
selected_t_idx = round(selected_t/dt) + 1;  

fig = figure('Position', [200 200 1000 600]);
hold on
for i = 1:length(selected_t)
    t_idx = selected_t_idx(i);
    plot(p_incentive_range, cluster_Up(t_idx, :), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('时间=%.1f小时', selected_t(i)))
end
hold off
% title('集群上调潜力随激励电价变化曲线（不同时间点）')  % 移除标题
xlabel('激励电价 (元)', 'FontSize', 16)
ylabel('上调潜力 (kW)', 'FontSize', 16)
grid on
lgd = legend('Location', 'best');
lgd.FontSize = 16;
ax = gca;
ax.FontSize = 16;
ylim([0, up_max*1.1])
print(fig, 'Cluster_Up_PriceCurve.png', '-dpng', '-r400');  % 400dpi保存

%% 11. 集群下调潜力随时间变化曲线（保存为PNG）
fig = figure('Position', [200 200 1000 600]);
hold on
for i = 1:length(selected_p_idx)
    p_idx = selected_p_idx(i);
    plot(time_points, cluster_Down(:, p_idx), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('激励电价=%.1f元', selected_p_values(i)))
end
hold off
% title('集群下调潜力随时间变化曲线（不同激励电价）')  % 移除标题
xlabel('时间 (小时)', 'FontSize', 16)
ylabel('下调潜力 (kW)', 'FontSize', 16)
grid on
lgd = legend('Location', 'best');
lgd.FontSize = 16;
ax = gca;
ax.FontSize = 16;
down_min = min(cluster_Down(:));
down_max = max(cluster_Down(:));
if down_min < down_max
    ylim([down_min*1.1, down_max*1.1])
else
    ylim([down_min*1.1, down_min*1.1 + 1e-6])
end
print(fig, 'Cluster_Down_TimeCurve.png', '-dpng', '-r400');  % 400dpi保存

%% 12. 集群下调潜力随激励电价变化曲线（保存为PNG）
fig = figure('Position', [200 200 1000 600]);
hold on
for i = 1:length(selected_t)
    t_idx = selected_t_idx(i);
    plot(p_incentive_range, cluster_Down(t_idx, :), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('时间=%.1f小时', selected_t(i)))
end
hold off
% title('集群下调潜力随激励电价变化曲线（不同时间点）')  % 移除标题
xlabel('激励电价 (元)', 'FontSize', 16)
ylabel('下调潜力 (kW)', 'FontSize', 16)
grid on
lgd = legend('Location', 'best');
lgd.FontSize = 16;
ax = gca;
ax.FontSize = 16;
if down_min < down_max
    ylim([down_min*1.1, down_max*1.1])
else
    ylim([down_min*1.1, down_min*1.1 + 1e-6])
end
print(fig, 'Cluster_Down_PriceCurve.png', '-dpng', '-r400');  % 400dpi保存

%% 13. 单体EV SOC随时间变化曲线（保存为PNG）
ev_index = 23;                     
p_incentive_index = 10;            
selected_p = p_incentive_range(p_incentive_index);  

SOC_data = results_3D.SOC_Individual(ev_index, :, p_incentive_index);

fig = figure('Position', [200 200 1000 500]);
plot(time_points, SOC_data, 'b-', 'LineWidth', 1.5)
% title(sprintf('EV %d SOC随时间变化曲线（激励电价=%.1f元）', ev_index, selected_p))  % 移除标题
xlabel('时间 (小时)', 'FontSize', 16)
ylabel('SOC', 'FontSize', 16)
grid on
ax = gca;
ax.FontSize = 16;  % 增大刻度标签字体
ylim([-1.5 1.5])
print(fig, 'SingleEV_SOC_Curve.png', '-dpng', '-r400');  % 400dpi保存

%% 16. 单体EV 上调调节能力曲线（保存为PNG）
Up_data = results_3D.EV_Up_Individual(ev_index, :, p_incentive_index);

fig = figure('Position', [200 200 1000 500]);
plot(time_points, Up_data, 'r-', 'LineWidth', 1.5)
% title(sprintf('EV %d 上调调节能力随时间变化曲线（激励电价=%.1f元）', ev_index, selected_p))  % 移除标题
xlabel('时间 (小时)', 'FontSize', 16)
ylabel('上调潜力 (kW)', 'FontSize', 16)
grid on
ax = gca;
ax.FontSize = 16;
up_min = min(Up_data);
up_max = max(Up_data);
if up_min < up_max
    ylim([up_min-5, up_max+5])
else
    ylim([up_min-5, up_min+5])
end
print(fig, 'SingleEV_Up_Curve.png', '-dpng', '-r400');  % 400dpi保存

%% 15. 单体EV 下调调节能力曲线（保存为PNG）
Down_data = results_3D.EV_Down_Individual(ev_index, :, p_incentive_index);

fig = figure('Position', [200 200 1000 500]);
plot(time_points, Down_data, 'g-', 'LineWidth', 1.5)
% title(sprintf('EV %d 下调调节能力随时间变化曲线（激励电价=%.1f元）', ev_index, selected_p))  % 移除标题
xlabel('时间 (小时)', 'FontSize', 16)
ylabel('下调潜力 (kW)', 'FontSize', 16)
grid on
ax = gca;
ax.FontSize = 16;
down_min_ind = min(Down_data);
down_max_ind = max(Down_data);
if down_min_ind < down_max_ind
    ylim([down_min_ind-5, down_max_ind+5])
else
    ylim([down_min_ind-5, down_min_ind+5])
end
print(fig, 'SingleEV_Down_Curve.png', '-dpng', '-r400');  % 400dpi保存

%% 16. 灵敏度分析模块（假设无绘图）
analyzeSensitivity(results_3D, p_incentive_range,...
    p_min, p_max, p_min_prime, p_max_prime);
