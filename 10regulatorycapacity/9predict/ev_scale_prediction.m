function [N_high, N_low] = ev_scale_prediction(params)
    % EV规模预测函数（基于Bass模型+淘汰机制）
    % 输入参数params结构体需包含：
    %   - start_year: 初始年份（如2022）
    %   - predict_years: 预测年数（例如，要预测到2030年，start_year=2022，则 predict_years = 2030 - 2022 + 1 = 9）
    %   - N0: 初始年的实际保有量（辆）
    %   - M: 潜在购买者数组（长度=predict_years，每年潜在用户数）
    %   - p_high, p_low: 高低情景创新系数 [cite: 48]
    %   - q_high, q_low: 高低情景模仿系数 [cite: 48]
    %   - r_ac: 偶然淘汰率（如0.07）[cite: 42]
    %   - Sur: 留存率数组（Sur(i)为第i年留存率，i=1,2,...）[cite: 44]
    % 输出：
    %   N_high: 高情景各年净保有量（辆），行向量，长度为 predict_years
    %   N_low: 低情景各年净保有量（辆），行向量，长度为 predict_years

    % 初始化参数
    T = params.predict_years;
    N_high_pred = zeros(1, T); % 用于存储预测结果
    N_low_pred = zeros(1, T);

    N_cum_high_hist = [params.N0]; % 存储累积购买量历史（用于计算）
    N_cum_low_hist = [params.N0];
    Sale_high_hist = []; % 存储历史销量
    Sale_low_hist = [];

    % --- 高情景计算 ---
    N_cum_high = params.N0; % 初始累积购买量
    N_net_high = params.N0; % 初始净保有量
    Sale_high = [];         % 存储销量序列（历史+预测）

    for t = 1:T
        year_index = t; % 当前是第 t 个预测年
        M_current = params.M(year_index); % 当年潜在购买者

        F_t_high = N_cum_high / M_current; % 使用上一年的累积购买量 [cite: 38]
        % Bass模型计算当年新增购买量 n(t+1) [cite: 37]
        n_t_high = M_current * (params.p_high + params.q_high * F_t_high) * (1 - F_t_high);
        Sale_high(year_index) = n_t_high; % 记录当年销量

        % 更新累积购买量（未淘汰前） N(t+1) = N(t) + n(t+1) [cite: 40]
        N_cum_high = N_cum_high + n_t_high;

        % 计算淘汰量
        E_A_high = N_cum_high * params.r_ac; % 偶然淘汰 [cite: 43]
        E_N_high = 0; % 自然淘汰 [cite: 45]
        Sales_combined_high = [Sale_high_hist, Sale_high]; % 合并历史和当前预测销量
        for i = 2:year_index % i 代表车辆已使用的年数
             sale_year_index = year_index - (i-1); % 销售发生的年份索引
             if sale_year_index >= 1 && (i-1) < length(params.Sur) && i <= length(params.Sur)
                 E_N_high = E_N_high + Sales_combined_high(sale_year_index) * (params.Sur(i-1) - params.Sur(i));
             elseif sale_year_index >= 1 && (i-1) <= length(params.Sur) && i > length(params.Sur) % 处理超过Sur数组长度的情况
                 E_N_high = E_N_high + Sales_combined_high(sale_year_index) * (params.Sur(i-1) - 0); % 假设完全淘汰
             end
        end

        % 计算净保有量 N'(t+1) = N(t+1) - EA(t) - EN(t) [cite: 47]
        % 注意：这里的N(t+1)是累积购买量，EA和EN也是基于累积量或历史销量计算
        N_net_high = N_cum_high - E_A_high - E_N_high;
        N_high_pred(t) = N_net_high; % 存储预测结果
    end

    % --- 低情景计算 (逻辑同高情景，仅替换 p_low, q_low) ---
    N_cum_low = params.N0;
    N_net_low = params.N0;
    Sale_low = [];

     for t = 1:T
        year_index = t;
        M_current = params.M(year_index);

        F_t_low = N_cum_low / M_current;
        n_t_low = M_current * (params.p_low + params.q_low * F_t_low) * (1 - F_t_low);
        Sale_low(year_index) = n_t_low;

        N_cum_low = N_cum_low + n_t_low;

        E_A_low = N_cum_low * params.r_ac;
        E_N_low = 0;
        Sales_combined_low = [Sale_low_hist, Sale_low];
         for i = 2:year_index
             sale_year_index = year_index - (i-1);
             if sale_year_index >= 1 && (i-1) < length(params.Sur) && i <= length(params.Sur)
                 E_N_low = E_N_low + Sales_combined_low(sale_year_index) * (params.Sur(i-1) - params.Sur(i));
              elseif sale_year_index >= 1 && (i-1) <= length(params.Sur) && i > length(params.Sur)
                 E_N_low = E_N_low + Sales_combined_low(sale_year_index) * (params.Sur(i-1) - 0);
             end
        end

        N_net_low = N_cum_low - E_A_low - E_N_low;
        N_low_pred(t) = N_net_low;
    end

    % 输出结果
    N_high = N_high_pred;
    N_low = N_low_pred;

    % 可选：在函数内部显示结果
    years = params.start_year : params.start_year + T - 1;
    disp('电动汽车高情景净保有量（辆）:');
    disp(array2table([years; N_high]', 'VariableNames',{'Year','Count'}));
    disp('电动汽车低情景净保有量（辆）:');
    disp(array2table([years; N_low]', 'VariableNames',{'Year','Count'}));
end