%% initializeACsFromExcel.m
function ACs = initializeACsFromExcel(filepath)
    % 读取Excel数据
    data = readtable(filepath);
    
    % 验证必要字段
    requiredFields = {'R','C','eta','Tset','Tmax','Tmin','SOC','p_incentive'};
    validateFields(data.Properties.VariableNames, requiredFields);
    
    % 转换数据为结构体数组
    ACs = table2struct(data);
    
    % 添加原始参数备份
    for i = 1:length(ACs)
            ACs(i).Tset_original = ACs(i).Tset;
            ACs(i).Tmax_original = ACs(i).Tmax;
            ACs(i).Tmin_original = ACs(i).Tmin;
    end

    % 嵌套字段验证函数
    function validateFields(actualFields, requiredFields)
        missingFields = setdiff(requiredFields, actualFields);
        if ~isempty(missingFields)
            error('缺少必要字段: %s', strjoin(missingFields, ', '));
        end
    end
end