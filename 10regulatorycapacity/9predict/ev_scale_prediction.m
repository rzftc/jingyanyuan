function [N_high, N_low] = ev_scale_prediction(params)
    % EV规模预测函数（基于Bass模型+淘汰机制）
    % 输入参数params结构体需包含：
    %   - start_year: 初始年份（如2022）
    %   - predict_years: 预测年数（如2030-2022+1=9）
    %   - N0: 初始累积购买量（辆）
    %   - M: 潜在购买者数组（长度=predict_years，每年潜在用户数）
    %   - p_high, p_low: 高低情景创新系数
    %   - q_high, q_low: 高低情景模仿系数
    %   - r_ac: 偶然淘汰率（如0.07）
    %   - Sur: 留存率数组（Sur(i)为第i年留存率，i=1,2,...）
    % 输出：
    %   N_high: 高情景各年净保有量（辆）
    %   N_low: 低情景各年净保有量（辆）
    
    % 初始化参数
    T = params.predict_years;
    N_high = zeros(1, T);
    N_low = zeros(1, T);
    N_high(1) = params.N0;  % 初始净保有量（需调整为初始实际保有量）
    N_low(1) = params.N0;
    
    % 保存每年销量（用于自然淘汰计算）
    Sale_high = zeros(1, T);
    Sale_low = zeros(1, T);
    
    % 逐年预测
    for t = 1:T
        % ---------------------- 高情景计算 ----------------------
        if t == 1
            % 初始年数据（假设初始累积购买量N(t)=N0）
            F_t_high = N_high(t) / params.M(t);  % F(t) = 累积购买量/潜在购买者
            n_t1_high = params.M(t) * (params.p_high + params.q_high * F_t_high) * (1 - F_t_high);  % 当年购买量
            Sale_high(t) = n_t1_high;
            N_cum_high = N_high(t) + n_t1_high;  % 累积购买量（未淘汰）
        else
            % 前一年累积购买量（未淘汰前）
            N_cum_prev_high = N_cum_high(t-1);
            F_t_high = N_cum_prev_high / params.M(t);
            n_t_high = params.M(t) * (params.p_high + params.q_high * F_t_high) * (1 - F_t_high);
            Sale_high(t) = n_t_high;
            N_cum_high(t) = N_cum_prev_high + n_t_high;
        end
        
        % 计算淘汰量
        % 被动淘汰（偶然淘汰）
        E_A_high = N_cum_high(t) * params.r_ac;
        % 自然淘汰（需累加前i年销量*留存率变化）
        E_N_high = 0;
        for i = 2:t  % 仅考虑有历史销量的年份
            if (t - i) >= 1  % 确保t-i年有销量数据
                E_N_high = E_N_high + Sale_high(t - i) * (params.Sur(i-1) - params.Sur(i));
            end
        end
        % 净保有量
        N_high(t) = N_cum_high(t) - E_A_high - E_N_high;
    
        % ---------------------- 低情景计算（同理，仅替换p/q） ----------------------
        if t == 1
            F_t_low = N_low(t) / params.M(t);
            n_t1_low = params.M(t) * (params.p_low + params.q_low * F_t_low) * (1 - F_t_low);
            Sale_low(t) = n_t1_low;
            N_cum_low = N_low(t) + n_t1_low;
        else
            N_cum_prev_low = N_cum_low(t-1);
            F_t_low = N_cum_prev_low / params.M(t);
            n_t_low = params.M(t) * (params.p_low + params.q_low * F_t_low) * (1 - F_t_low);
            Sale_low(t) = n_t_low;
            N_cum_low(t) = N_cum_prev_low + n_t_low;
        end
        
        E_A_low = N_cum_low(t) * params.r_ac;
        E_N_low = 0;
        for i = 2:t
            if (t - i) >= 1
                E_N_low = E_N_low + Sale_low(t - i) * (params.Sur(i-1) - params.Sur(i));
            end
        end
        N_low(t) = N_cum_low(t) - E_A_low - E_N_low;
    end
    
    % 调整输出格式（年份对应）
    years = params.start_year : params.start_year + T - 1;
    disp('电动汽车高情景净保有量（辆）:');
    disp([years; N_high]');
    disp('电动汽车低情景净保有量（辆）:');
    disp([years; N_low]');
end
