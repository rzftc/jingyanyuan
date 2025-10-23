function [N_high, N_low] = ac_scale_prediction(params)
    % AC_SCALE_PREDICTION 空调规模预测函数（支持Logistic和多元回归模型）
    %
    % 输入参数params结构体需包含：
    %   - method: 预测方法（'logistic'或'regression'）
    %   - N_hist: 历史保有量数据（行/列向量，如[2018,2019,...,2022]年的保有量，单位：台）
    %   - X_hist: 历史影响因子矩阵（仅当method='regression'时需要，每行对应一年，列=影响因子数）
    %   - X_pred: 未来影响因子预测值（仅当method='regression'时需要，行数=预测年数，列=影响因子数）
    %   - K_high, K_low: 高低情景饱和容量（台，仅当method='logistic'时需要）
    %   - r_high, r_low: 高低情景增长率（仅当method='logistic'时需要）
    %   - r_ac: 偶然淘汰率（如0.05，5%）
    %   - Sur: 留存率数组（Sur(i)为第i年留存率，i=1,2,...，如[1,0.92,0.85,...]）
    %   - predict_years: 预测年数（如8代表预测2023-2030年）
    %
    % 输出：
    %   N_high: 高情景各年净保有量（台，行向量）
    %   N_low: 低情景各年净保有量（台，行向量）

    %% 参数校验与预处理
    % 必选参数检查
    required_fields = {'method', 'N_hist', 'r_ac', 'Sur', 'predict_years'};
    for field = required_fields
        if ~isfield(params, field{1})
            error("缺少必选参数: %s", field{1});
        end
    end

    % 历史保有量转为列向量（统一处理行/列输入）
    N_hist = params.N_hist(:);  % 转为列向量（Nx1）
    T = params.predict_years;   % 预测年数
    N_high = zeros(1, T);       % 初始化高情景结果
    N_low = zeros(1, T);        % 初始化低情景结果

    % 历史销量计算（保有量变化量）
    Sale_hist = diff([0; N_hist]);  % 销量=当年保有量-前一年保有量（(N+1)x1 → Nx1）
    Sale_hist = Sale_hist';         % 转为行向量（1xN）
    Sale = [Sale_hist, zeros(1, T)]; % 扩展销量数组（历史+预测期，1x(N+T)）

    %% 方法选择：Logistic模型或多元回归
    switch lower(params.method)
        case 'logistic'
            % 校验Logistic模型所需参数
            logistic_fields = {'K_high', 'K_low', 'r_high', 'r_low'};
            for field = logistic_fields
                if ~isfield(params, field{1})
                    error("当method='logistic'时，缺少参数: %s", field{1});
                end
            end

            % 高情景预测（Logistic模型）
            N_high = logistic_prediction(params, N_hist, T, 'high', Sale);

            % 低情景预测（Logistic模型）
            N_low = logistic_prediction(params, N_hist, T, 'low', Sale);

        case 'regression'
            % 校验回归模型所需参数
            regression_fields = {'X_hist', 'X_pred'};
            for field = regression_fields
                if ~isfield(params, field{1})
                    error("当method='regression'时，缺少参数: %s", field{1});
                end
            end

            % 校验影响因子维度
            [n_hist, m] = size(params.X_hist);
            if n_hist ~= length(N_hist)
                error("X_hist的行数（%d）与N_hist的长度（%d）不匹配", n_hist, length(N_hist));
            end
            [n_pred, ~] = size(params.X_pred);
            if n_pred ~= T
                error("X_pred的行数（%d）与predict_years（%d）不匹配", n_pred, T);
            end

            % 高情景预测（多元回归模型）
            N_high = regression_prediction(params, N_hist, T, 'high', Sale);

            % 低情景预测（多元回归模型）
            N_low = regression_prediction(params, N_hist, T, 'low', Sale);

        otherwise
            error("无效的预测方法: %s，支持'logistic'或'regression'", params.method);
    end

    %% 输出结果（年份与保有量对应）
    [~, start_year] = get_hist_years(params.N_hist);  % 获取历史数据起始年
    pred_years = (start_year + length(N_hist)) : (start_year + length(N_hist) + T - 1);
    disp('===== 空调净保有量预测结果（台） =====');
    disp('年份 | 高情景 | 低情景');
    for i = 1:T
        fprintf('%d   | %.0f  | %.0f\n', pred_years(i), N_high(i), N_low(i));
    end
end

% -------------------------------------------------------------------------
% 子函数：Logistic模型预测核心逻辑
% -------------------------------------------------------------------------
function N_pred = logistic_prediction(params, N_hist, T, scenario, Sale)
    % 输入参数：
    %   params: 参数结构体
    %   N_hist: 历史保有量（列向量）
    %   T: 预测年数
    %   scenario: 'high'或'low'情景
    %   Sale: 历史+预测期销量数组（行向量，输入时为历史销量+初始0值，输出时更新预测期销量）

    % 提取情景参数
    K = params.(sprintf('K_%s', scenario));  % 饱和容量
    r = params.(sprintf('r_%s', scenario));  % 增长率
    r_ac = params.r_ac;                      % 偶然淘汰率
    Sur = params.Sur;                        % 留存率数组（如[1,0.92,0.85,...]）

    % 初始保有量（历史最后一年）
    N_prev = N_hist(end);
    N_pred = zeros(1, T);

    for t = 1:T
        % 步骤1：Logistic增长（未淘汰前）
        N_growth = N_prev + r * N_prev * (1 - N_prev / K);

        % 步骤2：计算淘汰量
        % 被动淘汰（偶然淘汰）
        E_A = N_growth * r_ac;

        % 自然淘汰（历史+预测期销量×留存率变化）
        E_N = 0;
        total_years = length(N_hist) + t;  % 总年数（历史年数+已预测年数）
        for i = 2:total_years  % i从第2年开始（第1年无留存率变化）
            sale_year = total_years - i;  % 对应销量的年份索引（从1开始）
            if sale_year >= 1  % 确保存在销量数据（sale_year=1对应第1年销量）
                % 留存率处理（避免越界）
                if (i-1) <= length(Sur)
                    sur_prev = Sur(i-1);
                else
                    sur_prev = 0;  % 超过留存率数组长度，视为完全淘汰
                end

                if i <= length(Sur)
                    sur_curr = Sur(i);
                else
                    sur_curr = 0;  % 超过留存率数组长度，视为完全淘汰
                end

                E_N = E_N + Sale(sale_year) * (sur_prev - sur_curr);
            end
        end

        % 步骤3：净保有量（增长后 - 被动淘汰 - 自然淘汰）
        N_pred(t) = N_growth - E_A - E_N;

        % 更新预测期销量（当年增长 - 前一年保有量）
        % 注意：MATLAB函数参数按值传递，需返回更新后的Sale（若主函数需要）
        % 此处若需主函数使用更新后的销量，需将Sale作为输出参数
        Sale(length(N_hist) + t) = N_growth - N_prev;

        % 更新前一年保有量（用于下一年计算）
        N_prev = N_pred(t);
    end
end


% -------------------------------------------------------------------------
% 子函数：多元回归模型预测核心逻辑（完整实现）
% -------------------------------------------------------------------------
function N_pred = regression_prediction(params, N_hist, T, scenario, Sale)
    % 输入参数：
    %   params: 参数结构体
    %   N_hist: 历史保有量（列向量）
    %   T: 预测年数
    %   scenario: 'high'或'low'情景（用于选择影响因子列）
    %   Sale: 历史+预测期销量数组（行向量，输入时为历史销量+初始0值，输出时更新预测期销量）

    % 提取影响因子（假设X_pred的奇数列对应高情景，偶数列对应低情景）
    % 注：实际应用中需根据用户数据结构调整，此处假设X_pred的列按[高情景因子1, 低情景因子1, 高情景因子2, 低情景因子2,...]排列
    m = size(params.X_hist, 2);  % 影响因子数（每个因子有高低情景）
    if strcmp(scenario, 'high')
        % 选择高情景影响因子列（第1,3,5...列）
        X_hist_scenario = params.X_hist(:, 1:2:end);
        X_pred_scenario = params.X_pred(:, 1:2:end);
    else
        % 选择低情景影响因子列（第2,4,6...列）
        X_hist_scenario = params.X_hist(:, 2:2:end);
        X_pred_scenario = params.X_pred(:, 2:2:end);
    end

    % 步骤1：拟合多元回归模型
    % 构造设计矩阵（含截距项）
    X_design = [ones(length(N_hist), 1), X_hist_scenario];  % 历史影响因子（Nx(m+1)）
    % 最小二乘估计回归系数
    beta = (X_design' * X_design) \ (X_design' * N_hist);  % 系数向量（(m+1)x1）

    % 步骤2：预测未淘汰前的保有量
    % 构造未来设计矩阵（含截距项）
    X_pred_design = [ones(T, 1), X_pred_scenario];  % 未来影响因子（Tx(m+1)）
    N_reg = X_pred_design * beta;  % 回归预测值（未淘汰前，Tx1）
    N_reg = N_reg';  % 转为行向量（1xT）

    % 步骤3：考虑淘汰机制（被动淘汰+自然淘汰）
    N_pred = zeros(1, T);
    r_ac = params.r_ac;  % 偶然淘汰率
    Sur = params.Sur;    % 留存率数组

    for t = 1:T
        % 被动淘汰量
        E_A = N_reg(t) * r_ac;

        % 自然淘汰量（历史+预测期销量×留存率变化）
        E_N = 0;
        total_years = length(N_hist) + t;  % 总年数（历史年数+已预测年数）
        for i = 2:total_years  % i从第2年开始（第1年无留存率变化）
            sale_year = total_years - i;  % 对应销量的年份索引（从1开始）
            if sale_year >= 1  % 确保存在销量数据（sale_year=1对应第1年销量）
                % 留存率处理（避免越界，使用标准if-else语句）
                if (i-1) <= length(Sur)
                    sur_prev = Sur(i-1);
                else
                    sur_prev = 0;  % 超过留存率数组长度，视为完全淘汰
                end

                if i <= length(Sur)
                    sur_curr = Sur(i);
                else
                    sur_curr = 0;  % 超过留存率数组长度，视为完全淘汰
                end

                E_N = E_N + Sale(sale_year) * (sur_prev - sur_curr);
            end
        end

        % 净保有量（回归预测值 - 被动淘汰 - 自然淘汰）
        N_pred(t) = N_reg(t) - E_A - E_N;

        % 更新预测期销量（当年回归预测值 - 前一年净保有量）
        if t == 1
            prev_N = N_hist(end);  % 第一年的前一年是历史最后一年
        else
            prev_N = N_pred(t-1);  % 后续年份的前一年是上一年净保有量
        end
        Sale(length(N_hist) + t) = N_reg(t) - prev_N;
    end
end


% -------------------------------------------------------------------------
% 辅助函数：获取历史数据年份（需用户提供年份信息，此处为示例）
% 注：实际应用中建议params包含历史数据对应的年份（如params.years_hist）
% -------------------------------------------------------------------------
function [years_hist, start_year] = get_hist_years(N_hist)
    % 示例：假设历史数据为最近length(N_hist)年（如2018-2022年）
    years_hist = (2022 - length(N_hist) + 1) : 2022;
    start_year = years_hist(1);
end
