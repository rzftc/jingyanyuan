function analyzeSensitivity(results, p_range, p_min, p_max, p_min_prime, p_max_prime)
    % 参数验证
    if nargin < 6
        error('缺少必要输入参数');
    end
    
    %% 维度预处理
    p_range = p_range(:).';  % 确保为行向量
    [T, P] = size(results.EV_Up);
    
    %% 改进1：增加蒙特卡洛模拟
    num_mc = 20; % 蒙特卡洛模拟次数
    mc_up = zeros(T, P, num_mc);
    mc_down = zeros(T, P, num_mc);
    
    parfor mc = 1:num_mc
        % 生成带噪声的扰动数据
        noise = 1 + 0.05*randn(T,P);
        noisy_up = results.EV_Up .* noise;
        noisy_down = results.EV_Down .* noise;
        
        % 计算梯度
        [dP_up, dP_down] = computeGradients(p_range, noisy_up, noisy_down);
        mc_up(:,:,mc) = dP_up;
        mc_down(:,:,mc) = dP_down;
    end
    
    %% 统计处理（增加时间维度平均）
    mean_up = squeeze(mean(mc_up, [1 3]));  % 平均时间和蒙特卡洛维度
    std_up = squeeze(std(mc_up, 0, [1 3]));
    mean_down = squeeze(mean(mc_down, [1 3]));
    std_down = squeeze(std(mc_down, 0, [1 3]));
    
    %% 改进2：临界点自适应标记
    critical_points = unique([p_min, p_max, p_min_prime, p_max_prime]);
    [~, idx] = min(abs(p_range - critical_points'), [], 2);
    critical_idx = unique(idx);
    
    %% 改进3：经济性分析
    cost_factor = 0.15; % 元/kW
    benefit_up = mean_up * cost_factor;
    benefit_down = abs(mean_down) * cost_factor * 0.8;
    
    %% 可视化增强
    plotSensitivityCurves(p_range, mean_up, std_up, mean_down, std_down,...
        p_min, p_max, p_min_prime, p_max_prime, critical_idx, benefit_up, benefit_down);
end