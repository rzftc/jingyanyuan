function [dist_squared, path] = dtw(seq1, seq2)
    % 动态时间规整算法实现（线性化版本）
    % 输入:
    %   seq1 - 序列1 [n×1]
    %   seq2 - 序列2 [m×1]
    % 输出:
    %   dist_squared - 最优路径上的局部距离平方和（式0-23/0-24）
    %   path         - 对齐路径坐标矩阵 [k×2]
    
    % 强制列向量输入
    seq1 = seq1(:); seq2 = seq2(:);
    n = length(seq1); m = length(seq2);
    
    % ===== 1. 构建局部距离矩阵 =====
    D = zeros(n, m);
    for i = 1:n
        for j = 1:m
            D(i,j) = (seq1(i) - seq2(j))^2;  % 式0-5
        end
    end
    
    % ===== 2. 初始化累积距离矩阵 =====
    C = inf(n, m);
    C(1,1) = D(1,1);  % 式0-6起点
    
    % 填充第一列（垂直移动）
    for i = 2:n
        C(i,1) = C(i-1,1) + D(i,1);  % 式0-6
    end
    
    % 填充第一行（水平移动）
    for j = 2:m
        C(1,j) = C(1,j-1) + D(1,j);  % 式0-7
    end
    
    % ===== 3. 递推计算累积距离 =====
    for i = 2:n
        for j = 2:m
            C(i,j) = D(i,j) + min([C(i-1,j),   % 垂直移动
                                   C(i-1,j-1), % 对角线移动 
                                   C(i,j-1)]); % 水平移动 式0-8
        end
    end
    
    % ===== 4. 回溯最优路径 =====
    path = zeros(n+m-1, 2); % 预分配路径内存
    k = 1;
    i = n; j = m;
    path(k,:) = [i, j];
    
    while i > 1 || j > 1
        k = k + 1;
        if i == 1
            j = j - 1; % 只能水平移动
        elseif j == 1
            i = i - 1; % 只能垂直移动
        else
            [~, idx] = min([C(i-1,j), C(i-1,j-1), C(i,j-1)]);
            switch idx
                case 1 % 垂直
                    i = i - 1;
                case 2 % 对角线
                    i = i - 1;
                    j = j - 1;
                case 3 % 水平
                    j = j - 1;
            end
        end
        path(k,:) = [i, j];
    end
    
    % 整理路径坐标（按时间正序排列）
    path = flipud(path(1:k,:));
    
    % ===== 5. 计算线性化距离量纲 =====
    dist_squared = sum(D(sub2ind(size(D), path(:,1), path(:,2))));  % 式0-23/0-24
end
