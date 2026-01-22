% main_spatial_analysis_from_excel.m
clear; clc; close all;

disp('开始从Excel进行空间特性分析的示例...');

% --- 步骤 0: 定义Excel文件名 ---
ev_excel_file = 'ev_data_example.xlsx';
ac_excel_file = 'ac_data_example.xlsx';
zones_excel_file = 'zone_definitions_example.xlsx';

% --- 步骤 1: 生成示例Excel文件 (如果不存在或需要覆盖) ---
% 将 overwrite_existing 设置为 true 会强制重新生成文件
% 通常在第一次运行时设为true，之后可以设为false或注释掉此行
generate_files_flag = false; % 或者 false，如果您不想每次都重新生成
if generate_files_flag
    disp('将生成或覆盖示例Excel文件。');
    generateSpatialAnalysisInputFiles(ev_excel_file, ac_excel_file, zones_excel_file, true);
else
    disp('检查示例Excel文件，如果不存在则生成。');
    generateSpatialAnalysisInputFiles(ev_excel_file, ac_excel_file, zones_excel_file, false);
end


% --- 步骤 2: 从Excel文件加载输入数据 ---
[ev_data, ac_data, zone_definitions] = loadSpatialAnalysisInputs(ev_excel_file, ac_excel_file, zones_excel_file);

% 检查加载的数据是否有效
if isempty(ev_data) && isempty(ac_data) && isempty(zone_definitions)
    error('所有输入数据均未能成功加载，请检查Excel文件和路径设置。示例终止。');
    return;
end
if isempty(zone_definitions)
    warning('功能区定义未能加载，空间分析可能不准确或失败。');
    % 即使功能区定义为空，也创建一个空的cell，以避免后续函数直接出错
    zone_definitions = {}; 
end


% --- 步骤 3: 分析EV空间特性 ---
% (需要确保 analyzeEVSpatialCharacteristics.m 在路径中)
if ~isempty(ev_data) && ~isempty(zone_definitions)
    disp(' '); % 添加空行以分隔输出
    disp('调用 analyzeEVSpatialCharacteristics...');
    ev_spatial_stats = analyzeEVSpatialCharacteristics(ev_data, zone_definitions);
    disp('EV空间特性分析结果:');
    disp(ev_spatial_stats);

    if isfield(ev_spatial_stats, 'Residential')
        disp(' ');
        disp('居民区 (Residential) EV统计详情:');
        disp(ev_spatial_stats.Residential);
        if isfield(ev_spatial_stats.Residential, 'charging_mode_distribution') && ~isempty(ev_spatial_stats.Residential.charging_mode_distribution)
            disp('居民区EV充电模式分布:');
            disp(ev_spatial_stats.Residential.charging_mode_distribution);
        end
        if isfield(ev_spatial_stats.Residential, 'start_time_hourly_distribution_prob') && ~isempty(ev_spatial_stats.Residential.start_time_hourly_distribution_prob)
            try % 添加try-catch以防绘图失败
                figure('Name', '居民区EV充电开始时间分布');
                bar(0:23, ev_spatial_stats.Residential.start_time_hourly_distribution_prob); % 假设概率对应0-23小时
                title('居民区EV充电开始时间分布概率');
                xlabel('小时'); ylabel('概率');
                grid on;
                xticks(0:2:23);
            catch ME_plot_ev
                warning( ME_plot_ev.message);
            end
        end
    end
else
    disp(' ');
    disp('EV数据或功能区定义不完整，跳过EV空间分析。');
end


% --- 步骤 4: 分析AC空间特性 ---
% (需要确保 analyzeACSpatialCharacteristics.m, aggregateSetpointCountsByArea.m, 
%  calculateSetpointProbabilitiesByArea.m 均在路径中)
if ~isempty(ac_data) && ~isempty(zone_definitions)
    disp(' '); % 添加空行以分隔输出
    disp('调用 analyzeACSpatialCharacteristics...');
    ac_spatial_stats = analyzeACSpatialCharacteristics(ac_data, zone_definitions);
    disp('AC空间特性分析结果:');
    disp(ac_spatial_stats);

    if isfield(ac_spatial_stats, 'Commercial')
        disp(' ');
        disp('商业区 (Commercial) AC统计详情:');
        disp(ac_spatial_stats.Commercial);
        if isfield(ac_spatial_stats.Commercial, 'setpoint_temperature_distribution') && ...
           isfield(ac_spatial_stats.Commercial.setpoint_temperature_distribution, 'P_TsetMax') && ...
           ~isempty(ac_spatial_stats.Commercial.setpoint_temperature_distribution.P_TsetMax)
            
            disp('商业区AC最高温度设定点分布:');
            commercial_tmax_dist_table = table(...
                ac_spatial_stats.Commercial.setpoint_temperature_distribution.TsetMax_values, ...
                ac_spatial_stats.Commercial.setpoint_temperature_distribution.P_TsetMax, ...
                'VariableNames', {'温度设定值', '概率'});
            disp(commercial_tmax_dist_table);
            
            try % 添加try-catch以防绘图失败
                figure('Name', '商业区AC最高温度设定点分布');
                bar(ac_spatial_stats.Commercial.setpoint_temperature_distribution.TsetMax_values, ...
                    ac_spatial_stats.Commercial.setpoint_temperature_distribution.P_TsetMax);
                title('商业区AC最高温度设定点分布概率');
                xlabel('温度设定值 (°C)'); ylabel('概率');
                grid on;
            catch ME_plot_ac
                 warning(ME_plot_ac.message);
            end
        else
            disp('商业区AC最高温度设定点数据不足或不存在。');
        end
    end
else
    disp(' ');
    disp('AC数据或功能区定义不完整，跳过AC空间分析。');
end

disp(' ');
disp('从Excel进行空间特性分析的示例结束。');