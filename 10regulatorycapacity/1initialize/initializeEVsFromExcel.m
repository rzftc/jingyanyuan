%% initializeEVsFromExcel.m
function EVs = initializeEVsFromExcel(filepath)
    % 读取Excel数据
    data = readtable(filepath);
    
    % 验证必要字段
    requiredFields = {'C_EV','eta','p_on','E_in','E_tar','t_in','t_dep','r'};
    validateFields(data.Properties.VariableNames, requiredFields);
    
    % 转换数据为结构体数组
    EVs = table2struct(data);
    
    % 补充默认值
    for i = 1:length(EVs)
       
            EVs(i).E_current = EVs(i).E_in;
        
        
            EVs(i).E_exp = EVs(i).E_in;
       
            EVs(i).P_current = 0;
       
            EVs(i).E_tar_original = EVs(i).E_tar;
        if isfield(EVs(i), 'SOC')
            EVs(i).SOC_original = EVs(i).SOC; % 备份原始SOC
        else
            % 如果确定 Excel 文件中可能没有 SOC 列，或者允许它缺失，
            % 你需要在这里定义一个默认的 SOC_original 值或抛出更具体的错误。
            % 但根据报错和依赖，SOC_original 是必需的。
            error('EVs结构体元素 %d 中缺少必要的 SOC 字段来创建 SOC_original。请检查Excel文件或数据生成脚本。', i);
        end
    end
    
    % 嵌套字段验证函数
    function validateFields(actualFields, requiredFields)
        missingFields = setdiff(requiredFields, actualFields);
        if ~isempty(missingFields)
            error('缺少必要字段: %s', strjoin(missingFields, ', '));
        end
    end
end