function copyMFilesWithParent()
    % 1. 让用户选择源文件夹
    srcFolder = uigetdir('', '请选择源文件夹 (Source Folder)');
    if srcFolder == 0
        disp('操作已取消。');
        return;
    end

    % 2. 让用户选择目标文件夹
    destFolder = uigetdir('', '请选择目标文件夹 (Destination Folder)');
    if destFolder == 0
        disp('操作已取消。');
        return;
    end
    
    % 检查源和目标是否相同
    if strcmp(srcFolder, destFolder)
        warning('源文件夹和目标文件夹不能相同。');
        return;
    end
    
    % --- 这是您要求的关键改动 ---
    % 3. 获取源文件夹本身的名称
    %    例如，如果 srcFolder = 'C:\MyProject'，则 srcFolderName = 'MyProject'
    [~, srcFolderName, ~] = fileparts(srcFolder);
    
    % 4. 在目标文件夹内创建这个“父文件夹”
    %    这现在是所有复制操作的基础目录
    finalDestRoot = fullfile(destFolder, srcFolderName);
    
    if ~exist(finalDestRoot, 'dir')
        mkdir(finalDestRoot);
    end
    % -------------------------------

    disp('正在开始搜索和复制...');
    disp(['复制目标将包含父文件夹: ', finalDestRoot]);
    
    % 5. 递归搜索所有的 .m 文件
    files = dir(fullfile(srcFolder, '**', '*.m'));
    
    if isempty(files)
        disp('在源文件夹及其子文件夹中未找到任何 .m 文件。');
        return;
    end
    
    copiedCount = 0;
    
    % 6. 遍历并复制文件
    for i = 1:length(files)
        if files(i).isdir
            continue;
        end
        
        % A. 获取源文件的完整路径
        sourceFile = fullfile(files(i).folder, files(i).name);
        
        % B. 从文件的完整目录中，移除源根目录的部分，得到相对路径
        %    例如: 'C:\MyProject\Subfolder' - 'C:\MyProject' = '\Subfolder'
        relativePath = erase(files(i).folder, srcFolder);
        
        % C. 构建目标子文件夹的路径
        %    !!! 注意: 这里使用的是 finalDestRoot，而不是 destFolder !!!
        destSubFolder = fullfile(finalDestRoot, relativePath);
        
        % D. 如果目标子文件夹不存在，则创建它
        if ~exist(destSubFolder, 'dir')
            mkdir(destSubFolder);
        end
        
        % E. 构建最终的目标文件路径
        destFile = fullfile(destSubFolder, files(i).name);
        
        % F. 复制文件
        [status, msg] = copyfile(sourceFile, destFile);
        
        if status
            copiedCount = copiedCount + 1;
        else
            warning('复制失败: %s. 错误: %s', sourceFile, msg);
        end
    end

    % 7. 结束报告
    disp('--------------------');
    disp('操作完成！');
    disp(['总共复制了 ', num2str(copiedCount), ' 个 .m 文件。']);
    disp(['目标根文件夹: ', finalDestRoot]);

end