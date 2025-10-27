%% 全功能虚拟电厂调节潜力分析系统
clear; close all; clc;

%% 1. 系统初始化
rng(2023);                                      % 固定随机种子
T_total = 24;                                   % 总时长(小时)
dt = 0.5;                                       % 时间分辨率
time_points = 0:dt:T_total;
base_price = 30;                                % 基础电价(元/kWh)

%% 2. 初始化参数
generateExampleExcel(100,100);
acFile = '../0input_data/AC_template.xlsx';
evFile = '../0input_data/EV_template.xlsx';

%% 3. 读取设备参数
ACs = initializeACsFromExcel(acFile);
EVs = initializeEVsFromExcel(evFile);
num_AC = length(ACs);                           % 空调数量
num_EV = length(EVs);                           % 电动汽车数量

%% 5. 预计算模块
% 5.1 EV基线功率预计算
for i = 1:num_EV
    [m1, m2, m3] = calculateEVABC_single(EVs(i).C_EV, EVs(i).eta,...
        EVs(i).E_tar, EVs(i).E_in, EVs(i).t_dep, EVs(i).t_in, dt, EVs(i).r);
    [~, EVs(i).P_base_sequence] = EVbaseP_single_longstep(m1, m2, m3,...
        EVs(i).p_on, EVs(i).SOC, length(time_points));
end

%% 6. 参数配置
p_min = 15; p_max = 50; 
p_min_prime = 10; p_max_prime = 45;
T_set_max = 3;                                  % 最大温度偏移
E_tar_max = 0.2 * [EVs.C_EV];                   % 计算最大目标电量变化

%% 新增：初始化存储矩阵
p_values = p_min_prime:p_max;
num_p = length(p_values);

% 调节能力存储矩阵
aggregate_AC = struct('Up',zeros(num_p,1), 'Down',zeros(num_p,1));
single_AC = struct('Up',zeros(num_p,1), 'Down',zeros(num_p,1));
aggregate_EV = struct('Up',zeros(num_p,1), 'Down',zeros(num_p,1));
single_EV = struct('Up',zeros(num_p,1), 'Down',zeros(num_p,1));

%% 主循环（按激励电价遍历）
for idx = 1:num_p
    p_incentive = p_values(idx);
    
    %% 4. 激励响应模块
    % 4.1 AC温度设定调整
    for i = 1:num_AC
        participation = calculateParticipation(p_incentive, base_price);
        [~, ~, deltaT] = incentiveTempAC(p_incentive, p_min, p_max,...
            p_min_prime, p_max_prime, T_set_max);
        
        ACs(i).ptcp = (rand() < participation); 
        
        if ACs(i).ptcp
            ACs(i).Tmax = ACs(i).Tset + deltaT;
            ACs(i).Tmin = ACs(i).Tset - deltaT;
        end
        
        base_temp = ACs(i).Tset + 4*sin(2*pi*time_points/24); 
        temp_range = ACs(i).Tmax - ACs(i).Tmin;
        noise = 0.2 * temp_range * randn(size(time_points));
        ACs(i).T_ja = min(max(base_temp + noise, ACs(i).Tmin), ACs(i).Tmax);
    end 
    
    % 4.2 EV目标电量调整
    for i = 1:num_EV
        participation = calculateParticipation(p_incentive, base_price);
        [~, ~, deltaE] = incentiveTempEV(p_incentive, p_min, p_max,...
            p_min_prime, p_max_prime, E_tar_max(i));
        ptcp_result = (rand() < participation); 
        EVs(i).ptcp = ptcp_result;
        
        if EVs(i).ptcp
            EVs(i).E_tar = EVs(i).E_tar - deltaE;
            if EVs(i).E_tar <= EVs(i).E_in
                EVs(i).E_tar = EVs(i).E_tar + deltaE;
            end
        end
    end
    
    %% 6.1 空调调节潜力计算
    AC_Up_total = 0; AC_Down_total = 0;
    for i = 1:num_AC
        if ACs(i).ptcp
            ACs(i).P_base = ACbaseP_single(ACs(i).T_ja(12), ACs(i).Tset,...
                ACs(i).R, ACs(i).eta);
            
            [alpha, beta, gamma] = calculateACABC_single(ACs(i).R, ACs(i).C,...
                ACs(i).eta, ACs(i).Tmax, ACs(i).Tmin, ACs(i).Tset, dt);
            
            [DeltaP_plus, DeltaP_minus] = calculateACAdjustmentPotentia(...
                ACs(i).P_base, 2*abs(ACs(i).P_base), 0, alpha, beta, gamma,...
                0.5, dt*(0.95+0.1*rand()));
            
            % 累计群体调节能力
            AC_Up_total = AC_Up_total + DeltaP_plus;
            AC_Down_total = AC_Down_total + DeltaP_minus;
            
            % 记录第一个单体设备
            if i == 1
                single_AC.Up(idx) = DeltaP_plus;
                single_AC.Down(idx) = DeltaP_minus;
            end
        end
    end
    aggregate_AC.Up(idx) = AC_Up_total;
    aggregate_AC.Down(idx) = AC_Down_total;
    
    %% 6.2 电动汽车调节潜力计算
    EV_Up_total = 0; EV_Down_total = 0;
    for i = 1:num_EV
        if EVs(i).ptcp
            EVs(i).P_base = EVs(i).P_base_sequence(12);
            
            [~, ~, m3] = calculateEVABC_single(EVs(i).C_EV, EVs(i).eta,...
                EVs(i).E_tar, EVs(i).E_in, EVs(i).t_dep, EVs(i).t_in, dt, EVs(i).r);
            
            [DeltaP_plus, DeltaP_minus] = calculateEVAdjustmentPotentia(...
                EVs(i).C_EV, EVs(i).r,EVs(i).eta, EVs(i).E_tar, EVs(i).E_in,EVs(i).E_current,...
                EVs(i).t_dep, EVs(i).t_in, EVs(i).p_on, 0,...
                EVs(i).P_base, EVs(i).SOC, dt);
            
            % 累计群体调节能力
            EV_Up_total = EV_Up_total + DeltaP_plus;
            EV_Down_total = EV_Down_total + DeltaP_minus;
            
            % 记录第一个单体设备
            if i == 1
                single_EV.Up(idx) = DeltaP_plus;
                single_EV.Down(idx) = DeltaP_minus;
            end
        end
    end
    aggregate_EV.Up(idx) = EV_Up_total;
    aggregate_EV.Down(idx) = EV_Down_total;
end

%% 7. 绘图模块
figure('Position', [100 100 1200 800])

% AC调节能力曲线
subplot(2,1,1)
plot(p_values, aggregate_AC.Up, 'b-o', 'LineWidth', 1.5, 'DisplayName','群体上调潜力')
hold on
plot(p_values, aggregate_AC.Down, 'r-s', 'LineWidth', 1.5, 'DisplayName','群体下调潜力')
plot(p_values, single_AC.Up, 'b--', 'LineWidth', 1, 'DisplayName','单体上调潜力')
plot(p_values, single_AC.Down, 'r--', 'LineWidth', 1, 'DisplayName','单体下调潜力')
title('空调调节能力 vs 激励电价')
xlabel('激励电价 (元/kWh)')
ylabel('调节能力 (kW)')
legend('Location','best')
grid on

% EV调节能力曲线
subplot(2,1,2)
plot(p_values, aggregate_EV.Up, 'g-o', 'LineWidth', 1.5, 'DisplayName','群体上调潜力')
hold on
plot(p_values, aggregate_EV.Down, 'm-s', 'LineWidth', 1.5, 'DisplayName','群体下调潜力')
plot(p_values, single_EV.Up, 'g--', 'LineWidth', 1, 'DisplayName','单体上调潜力')
plot(p_values, single_EV.Down, 'm--', 'LineWidth', 1, 'DisplayName','单体下调潜力')
title('电动汽车调节能力 vs 激励电价')
xlabel('激励电价 (元/kWh)')
ylabel('调节能力 (kW)')
legend('Location','best')
grid on

%% 辅助函数（需确保在路径中）
% initializeACsFromExcel, initializeEVsFromExcel, calculateEVABC_single 等函数需正常实现
