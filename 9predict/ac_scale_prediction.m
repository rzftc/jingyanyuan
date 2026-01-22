function [N_high, N_low] = ac_scale_prediction(params)
    % AC_SCALE_PREDICTION 空调规模预测函数（支持Logistic和多元回归模型）
    % [版本 V5.1: 修正字段名称访问错误 + 显式参数传递]

    %% 参数校验与预处理 (同前)
    required_fields = {'method', 'N_hist', 'r_ac', 'Sur', 'predict_years'};
    for field = required_fields
        if ~isfield(params, field{1}) || isempty(params.(field{1}))
             error("缺少必选参数或参数为空: %s", field{1});
        end
    end
    if ~isnumeric(params.Sur) || isempty(params.Sur) || params.Sur(1) ~= 1
        error("Sur 数组必须是非空数值向量且 Sur(1) 必须等于 1");
    end
    if ~isnumeric(params.N_hist) || isempty(params.N_hist)
         error("N_hist 必须是非空数值向量");
    end

    N_hist = params.N_hist(:);
    T = params.predict_years;
    if T <= 0
        N_high = []; N_low = [];
        disp('预测年数为 0 或负数，不进行预测。');
        return;
    end
    N_high_pred = zeros(1, T);
    N_low_pred = zeros(1, T);
    Sale_hist = max(0, diff([0; N_hist])');
    Sale_full = [Sale_hist, zeros(1, T)];

    %% 方法选择 (调用核心函数 V5 - 使用显式参数传递)
    try
        switch lower(params.method)
            case 'logistic'
                logistic_fields = {'K_high', 'K_low', 'r_high', 'r_low'};
                 for field = logistic_fields
                    if ~isfield(params, field{1}) || isempty(params.(field{1}))
                        error("当method='logistic'时，缺少参数: %s", field{1});
                    end
                 end
                % !!! V5.1 修正: 显式传递参数 !!!
                [N_high_pred, ~] = prediction_core_v5(params, N_hist, T, Sale_full, 'high', 'logistic');
                [N_low_pred, ~] = prediction_core_v5(params, N_hist, T, Sale_full, 'low', 'logistic');

            case 'regression'
                 regression_fields = {'X_hist', 'X_pred'};
                 % ... (parameter and dimension checks as before) ...
                 % !!! V5.1 修正: 显式传递参数 !!!
                 [N_high_pred, ~] = prediction_core_v5(params, N_hist, T, Sale_full, 'high', 'regression');
                 [N_low_pred, ~] = prediction_core_v5(params, N_hist, T, Sale_full, 'low', 'regression');
            otherwise
                error("无效的预测方法: %s", params.method);
        end
    catch ME
        fprintf('预测核心函数出错: %s\n', ME.message);
        fprintf('发生在文件 %s 的第 %d 行\n', ME.stack(1).file, ME.stack(1).line);
        N_high = nan(1, T); N_low = nan(1, T); % 返回NaN表示失败
        return; % 提前退出
    end

    %% 输出结果 (同前)
    N_high = N_high_pred;
    N_low = N_low_pred;
    % ... (显示结果部分，同前) ...
     if isfield(params, 'hist_years') && length(params.hist_years) == length(N_hist)
        start_year = params.hist_years(1);
        hist_len = length(params.hist_years);
    else
        start_year = year(now) - length(N_hist);
        hist_len = length(N_hist);
        warning('未提供历史年份 (params.hist_years)，假设历史数据截止到 %d 年', start_year + hist_len - 1);
    end
    pred_years = (start_year + hist_len) : (start_year + hist_len + T - 1);

    if all(isnan(N_high)) && all(isnan(N_low))
        disp('===== 空调净保有量预测失败 =====');
    else
        disp('===== 空调净保有量预测结果（台） =====');
        valid_idx = ~isnan(N_high) & ~isnan(N_low);
        if any(valid_idx)
             results_table = table(pred_years(valid_idx)', round(N_high(valid_idx))', round(N_low(valid_idx))', 'VariableNames', {'Year', 'High_Scenario', 'Low_Scenario'});
             disp(results_table);
        else; disp('所有预测结果均为 NaN。'); end
    end
end

% --- 核心预测子函数 V5: 修正字段访问，保持V4逻辑 (保持不变) ---
function [N_pred, Sale_updated] = prediction_core_v5(params, N_hist, T, Sale_in, scenario, method)
    % --- 输入验证 ---
     if ~ismember(lower(scenario), {'high', 'low'})
        error('Scenario 参数必须是 ''high'' 或 ''low''。');
     end

    r_ac = params.r_ac;
    Sur = params.Sur;
    N_hist_len = length(N_hist);
    N_pred = zeros(1, T);
    Sale_updated = Sale_in;

    % --- 模型特定参数准备 ---
    beta = []; X_pred_scenario = []; K = 0; r_log = 0;
    if strcmpi(method, 'regression')
        % ... (regression setup as in V4) ...
         m_hist = size(params.X_hist, 2);
        if strcmp(scenario, 'high') % 选择情景因子
            X_hist_scenario = params.X_hist(:, 1:2:end); X_pred_scenario = params.X_pred(:, 1:2:end);
        else
            X_hist_scenario = params.X_hist(:, 2:2:end); X_pred_scenario = params.X_pred(:, 2:2:end);
        end
        X_design = [ones(N_hist_len, 1), X_hist_scenario];
        try beta = (X_design' * X_design) \ (X_design' * N_hist);
        catch ME_fit; error('回归模型拟合失败: %s. 请检查 X_hist 数据。', ME_fit.message); end

    elseif strcmpi(method, 'logistic')
        % !!! V5 修正: 使用显式 if 语句访问字段 !!!
        if strcmpi(scenario, 'high')
            K = params.K_high;
            r_log = params.r_high;
        else % scenario is 'low'
            K = params.K_low;
            r_log = params.r_low;
        end
        % --- 结束 V5 修正 ---
         if K <= N_hist(end)
             warning('Logistic饱和容量 K (%g) 不大于最新的历史保有量 (%g)，可能导致预测不增长或下降。', K, N_hist(end));
         end
    end

    N_prev = N_hist(end);

    for t = 1:T
        current_abs_year_index = N_hist_len + t;

        % --- Step 1: 计算总淘汰量 Total_Scrappage(t) ---
        % 1.1 自然淘汰 E_N (同 V4)
        Natural_Scrappage = 0;
        max_age_at_start_of_year = current_abs_year_index - 1;
        Sur_len = length(Sur);
        for k = 1:max_age_at_start_of_year
            sale_year_abs_index = current_abs_year_index - k;
            if sale_year_abs_index >= 1
                sur_k_minus_1 = 0; if k == 1; sur_k_minus_1 = 1; elseif (k-1) > 0 && (k-1) < Sur_len; sur_k_minus_1 = Sur(k-1); end
                sur_k = 0; if k > 0 && k < Sur_len; sur_k = Sur(k); end
                fraction_scrapped_naturally = max(0, sur_k_minus_1 - sur_k);
                Natural_Scrappage = Natural_Scrappage + Sale_updated(sale_year_abs_index) * fraction_scrapped_naturally;
            end
        end
        % 1.2 偶然淘汰 E_A (同 V4)
        Accidental_Scrappage = N_prev * r_ac;
        % 1.3 总淘汰
        Total_Scrappage = Natural_Scrappage + Accidental_Scrappage;

        % --- Step 2: 计算当期新增销量 New_Sales(t) (同 V4) ---
        New_Sales = 0;
        if strcmpi(method, 'logistic')
            Net_Increase = max(0, r_log * N_prev * (1 - N_prev / K));
            Replacement_Demand = Total_Scrappage;
            New_Sales = Net_Increase + Replacement_Demand;
        elseif strcmpi(method, 'regression')
            X_pred_design_t = [1, X_pred_scenario(t, :)];
            N_net_target_regression = X_pred_design_t * beta;
            Stock_after_scrappage = N_prev - Total_Scrappage;
            New_Sales = max(0, N_net_target_regression - Stock_after_scrappage);
        end
         % 添加一个检查，确保销量不会异常大 (例如，超过饱和容量)
         if isfinite(K) && K > 0 && New_Sales > K
             warning('计算出的年销量 (%g) 异常大，可能参数设置有误 (t=%d, N_prev=%g)', New_Sales, t, N_prev);
             New_Sales = K / 10; % 临时限制
         end
        Sale_updated(current_abs_year_index) = New_Sales;


        % --- Step 3: 计算当期末净保有量 N_net(t) (同 V4) ---
        N_net_pred = N_prev - Total_Scrappage + New_Sales;
        N_pred(t) = max(0, N_net_pred); % 确保非负

        % --- Step 4: 更新 N_prev (同 V4) ---
        N_prev = N_pred(t);
    end
end