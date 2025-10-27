%% 全功能虚拟电厂调节潜力分析系统
clear; close all; clc;
%% 参数生成应用实例
% generateEVData(100, 'residential')  % 生成100辆居民区EV
% generateEVData(50, 'workplace', {'特斯拉Model Y'}) % 生成50辆工作区特斯拉
%% 1. 系统初始化
rng(2023, 'Threefry');                         % 增强随机数生成器
T_total = 24;                                   % 总时长（小时）
dt = 0.05;                                      % 时间分辨率（小时）
time_points = 0:dt:T_total;                     
base_price = 30;            
%% 2. 并行计算初始化
if isempty(gcp('nocreate'))                    % 检测现有并行池
    parpool('local',24);                    % 创建4进程池
else
    disp('已有并行池运行中');
end
%% 2. 初始化参数
generateEVData(1000, 'residential', 'all');
evFile = 'EV_居民区.xlsx';                 % 电动汽车数据文件


%% 3. 读取设备参数
EVs = initializeEVsFromExcel(evFile);           % 从Excel导入电动汽车参数
num_EV = length(EVs);                           % 获取电动汽车总数

%% 4. 激励响应模块
%% 4.1 参数设定
p_min = 15; p_max = 50;                         % 原始电价范围
p_min_prime = 10; p_max_prime = 40;              % 调整后电价范围

%% 新增模块：激励电价参数扫描
p_incentive_range = linspace(0, 50, 25);        % 生成0-50之间的25个激励电价值
all_EV_Up = zeros(length(time_points), length(p_incentive_range));
all_EV_Down = zeros(length(time_points), length(p_incentive_range));
%% 参数扫描主循环
parfor (p_idx = 1:length(p_incentive_range), 12)
    current_p = p_incentive_range(p_idx);
    fprintf('\n== 分析激励电价 %.1f 元 ==\n', current_p);
    
    %% 4.2 EV目标电量调整（动态版本）
    EVs_temp = copyEVStruct(EVs);              % 关键深拷贝操作
    
    %% 6.2 目标电量调整（独立随机数流）
    stream = RandStream('Threefry', 'Seed', 2023+p_idx);
    RandStream.setGlobalStream(stream);
    E_tar_max = 0.2 * [EVs_temp.C_EV];
    
    for i = 1:num_EV
        EVs_temp(i).p_incentive = current_p;
        participation = calculateParticipation(current_p, base_price);
        [~, ~, deltaE] = incentiveTempEV(current_p, p_min, p_max, p_min_prime, p_max_prime, E_tar_max(i));
        ptcp_result = (rand() < participation);
        EVs_temp(i).ptcp = ptcp_result;
        
        if EVs_temp(i).ptcp
            EVs_temp(i).E_tar = EVs_temp(i).E_tar - deltaE;
            if EVs_temp(i).E_tar <= EVs_temp(i).E_in
                EVs_temp(i).E_tar = EVs_temp(i).E_tar + deltaE;
            end
        end
    end
    
    %% 5. 预计算模块（完整保留）
    H = 24;                                      % 预测时域（小时）
    H_steps = H / dt;                            % 转换为时间步数
    for i = 1:num_EV
        [~, EVs_temp(i).P_base_sequence] = EVbaseP_single_longstep(...
            EVs_temp(i).C_EV, EVs_temp(i).eta,...
            EVs_temp(i).E_tar, EVs_temp(i).E_in,...
            EVs_temp(i).t_dep, EVs_temp(i).t_in, dt, EVs_temp(i).r,...
            EVs_temp(i).p_on, EVs_temp(i).SOC, H_steps+1);
    end
    
    %% 6. 主时间循环（完整保留）
    local_Up = zeros(length(time_points),1);  % 本地临时变量
    local_Down = zeros(length(time_points),1);
    local_SOC = zeros(length(time_points),num_EV);
    for t_idx = 1:length(time_points)
        t = time_points(t_idx);
        fprintf('\n时间 %.1f小时 [激励%.1f元]', t, current_p);
        
        for i = 1:num_EV
            online = (t >= EVs_temp(i).t_in) && (t <= EVs_temp(i).t_dep);
            if online
                EVs_temp(i).P_base = EVs_temp(i).P_base_sequence(t_idx);
                
                [~, ~, m3] = calculateEVABC_single(...
                    EVs_temp(i).C_EV, EVs_temp(i).eta,...
                    EVs_temp(i).E_tar, EVs_temp(i).E_in,...
                    EVs_temp(i).t_dep, EVs_temp(i).t_in, dt, EVs_temp(i).r);
                
                [EVs_temp(i).E_exp, EVs_temp(i).E_current,...
                    EVs_temp(i).P_current, EVs_temp(i).SOC] = ...
                    calculateEVS_single(m3, EVs_temp(i).E_exp,...
                    EVs_temp(i).eta, EVs_temp(i).E_current,...
                    EVs_temp(i).P_current, EVs_temp(i).C_EV,...
                    EVs_temp(i).r, EVs_temp(i).p_on, dt);
                
                [DeltaP_plus, DeltaP_minus] = calculateEVAdjustmentPotentia(...
                    EVs_temp(i).C_EV, EVs_temp(i).r, EVs_temp(i).eta,...
                    EVs_temp(i).E_tar, EVs_temp(i).E_in,...
                    EVs_temp(i).E_current, EVs_temp(i).t_dep,...
                    EVs_temp(i).t_in, EVs_temp(i).p_on, 0,...
                    EVs_temp(i).P_base, EVs_temp(i).SOC, dt);
                
                local_Up(t_idx) = sum(DeltaP_plus);  % 集群潜力累加
                local_Down(t_idx) = sum(DeltaP_minus);
                local_SOC(i, t_idx) = EVs_temp(i).SOC;
            else
                local_SOC(i, t_idx) = EVs_temp(i).SOC;
            end
        end
    end
    all_EV_Up(:, p_idx) = local_Up;  
    all_EV_Down(:, p_idx) = local_Down;

end
 results_3D.EV_Up = all_EV_Up;
 results_3D.EV_Down = all_EV_Down;
%% 7. 可视化模块（新增三维可视化）
figure('Position', [200 200 1400 600])

% 7.1 三维上调潜力曲面
subplot(1,2,1)
[X, Y] = meshgrid(time_points, p_incentive_range);
surf(X, Y, results_3D.EV_Up', 'EdgeColor','none')
title('上调潜力三维分布'), xlabel('时间 (小时)'), ylabel('激励电价 (元)')
zlabel('调节潜力 (kW)'), colormap(jet), colorbar
view(30, 60)

% 7.2 三维下调潜力曲面
subplot(1,2,2)
surf(X, Y, results_3D.EV_Down', 'EdgeColor','none')
title('下调潜力三维分布'), xlabel('时间 (小时)'), ylabel('激励电价 (元)')
zlabel('调节潜力 (kW)'), colormap(jet), colorbar
view(30, 60)

%% 8. 灵敏度分析模块（新增部分）
analyzeSensitivity(results_3D, p_incentive_range,...
    p_min, p_max, p_min_prime, p_max_prime);

%% 8. 资源释放模块
% 关闭并行计算池
% if ~isempty(gcp('nocreate'))
%     delete(gcp('nocreate'));
%     disp('并行池已关闭');
% end
%% 保留原有二维可视化
% figure('Position', [100 100 1200 800])
% subplot(2,1,1)
% plot(time_points, results.EV_Up, 'r', time_points, results.EV_Down, 'b')
% title('典型激励场景调节潜力'), grid on
% legend('上调潜力', '下调潜力')
% 
% subplot(2,1,2)
% imagesc(time_points, 1:num_EV, results.SOC_EV)
% title('SOC时空分布'), xlabel('时间 (小时)'), ylabel('车辆编号')
% colorbar




