
clear; close all; clc;

% --- 1. 定义网络参数 (源自 test_suite_comprehensive_dispatch.m) ---
N_bus = 5;
% 资源分布比例
AcDist = [0.2; 0.2; 0.2; 0.2; 0.2]; % AC 均匀分布在所有节点
EvDist = [0; 0.4; 0; 0.6; 0];       % EV 仅分布在节点 2 和 4

% --- 2. 构建拓扑结构 (示意性 5节点 6线路) ---
% 构建一个包含环路和交叉的典型拓扑
% Line 1: 1-2
% Line 2: 2-3
% Line 3: 3-4
% Line 4: 4-5
% Line 5: 5-1
% Line 6: 2-4 (交叉线路，使EV节点互联)
sources = [1 2 3 4 5 2];
targets = [2 3 4 5 1 4];
G = graph(sources, targets);

% --- 3. 绘图设置 ---
figure('Color', 'w', 'Position', [200, 200, 800, 600]);

% 使用力导向布局，但固定节点位置以获得更好看的形状（五边形）
theta = linspace(90, 90+360, N_bus+1)'; 
theta(end) = []; % 去掉重复点
x_pos = cosd(theta);
y_pos = sind(theta);

h = plot(G, 'XData', x_pos, 'YData', y_pos, ...
    'LineWidth', 2, 'EdgeColor', [0.6 0.6 0.6], ...
    'NodeFontSize', 12, 'NodeFontWeight', 'bold');

axis off; % 关闭坐标轴
hold on;

% --- 4. 可视化资源分布 ---

% 4.1 标记纯 AC 节点 (蓝色)
idx_only_ac = find(EvDist == 0);
highlight(h, idx_only_ac, 'NodeColor', [0.0, 0.447, 0.741], ... % 蓝色
    'MarkerSize', 15, 'Marker', 'o');

% 4.2 标记 AC + EV 复合节点 (红色，强调EV接入点)
idx_ev_ac = find(EvDist > 0);
highlight(h, idx_ev_ac, 'NodeColor', [0.850, 0.325, 0.098], ... % 橙红色
    'MarkerSize', 18, 'Marker', 'hexagram');

% --- 5. 添加详细注释标签 ---
for i = 1:N_bus
    % 构建节点说明文本
    info_str = sprintf('节点 %d', i);
    
    % 添加 AC 信息
    if AcDist(i) > 0
        info_str = [info_str, sprintf('\nAC: %.0f%%', AcDist(i)*100)];
    end
    
    % 添加 EV 信息
    if EvDist(i) > 0
        info_str = [info_str, sprintf('\nEV: %.0f%%', EvDist(i)*100)];
    end
    
    % 计算标签偏移量（向外偏移）
    offset = 0.15;
    text(x_pos(i) * (1 + offset), y_pos(i) * (1 + offset), info_str, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 11, 'FontName', 'Microsoft YaHei', ...
        'EdgeColor', 'none', 'BackgroundColor', 'none', 'FontWeight', 'bold');
end

% --- 6. 添加图例 (手动构建) ---
% 绘制不可见的点用于生成图例
scatter(nan, nan, 100, [0.0, 0.447, 0.741], 'filled', 'o', 'DisplayName', '普通负荷节点 (仅AC)');
scatter(nan, nan, 100, [0.850, 0.325, 0.098], 'filled', 'h', 'DisplayName', '灵活性资源节点 (AC + EV)');
legend('Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 12, 'FontName', 'Microsoft YaHei', 'Box', 'off');

% --- 7. 保存图片 ---
% 无标题，高DPI，中文名
print(gcf, '算例电网拓扑与资源分布.png', '-dpng', '-r300');
fprintf('图片已保存为: 算例电网拓扑与资源分布.png\n');