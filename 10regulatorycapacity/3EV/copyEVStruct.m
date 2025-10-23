function newEVs = copyEVStruct(EVs)
    % 获取字段模板确保一致性
    template = EVs(1);
    required_fields = fieldnames(template);
    
    % 预分配带完整字段的结构体数组
    newEVs = repmat(template, size(EVs));
    
    % 深度拷贝每个元素
    for i = 1:numel(EVs)
        % 验证源结构体字段完整性
        if ~isequal(fieldnames(EVs(i)), required_fields)
            error('结构体元素%d字段不匹配',i);
        end
        
        % 逐字段复制
        for f = 1:length(required_fields)
            field = required_fields{f};
            value = EVs(i).(field);
            
            % 增强型拷贝逻辑
            if isobject(value)
                newEVs(i).(field) = copy(value);  % 处理MATLAB对象
            elseif isstruct(value)
                newEVs(i).(field) = copyEVStruct(value); % 递归处理嵌套结构
            else
                newEVs(i).(field) = value;        % 基本类型直接赋值
            end
        end
    end
end