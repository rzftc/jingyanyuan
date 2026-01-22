function copyMFiles()
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
    
    % 3. 获取源文件夹本身的名称
    [~, srcFolderName, ~] = fileparts(srcFolder);
    
    % 4. 在目标文件夹内创建这个“父文件夹”
    %    这现在是所有复制操作的基础目录
    finalDestRoot = fullfile(destFolder, srcFolderName);
    
    % -------------------- 【修改开始】 --------------------
    % 5. 检查目标文件夹是否已存在，如果存在则删除
    if exist(finalDestRoot, 'dir')
        fprintf('目标文件夹 "%s" 已存在。\n', finalDestRoot);
        fprintf('正在删除旧文件夹及其所有内容...\n');
        try
            [status, msg, msgID] = rmdir(finalDestRoot, 's');
            if status
                fprintf('旧文件夹已成功删除。\n');
            else
                warning('删除旧文件夹失败: %s. \n错误ID: %s\n复制操作将继续，但这可能导致文件覆盖或合并。', msg, msgID);
            end
        catch ME
            % warning('删除旧文件夹时发生严重错误: %s\n复制操作将继续，但这可能导致文件覆盖或合并。', ME.message);
        end
    end
    
    % 6. 重新创建空的 finalDestRoot 文件夹
    fprintf('正在创建目标根文件夹: %s\n', finalDestRoot);
    try
        mkdir(finalDestRoot);
    catch ME
        error('无法创建目标文件夹: %s. \n请检查权限。操作终止。', ME.message);
    end
    % -------------------- 【修改结束】 --------------------

    disp('正在开始搜索和复制...');
    
    % 7. 递归搜索所有的 .m 文件
    files = dir(fullfile(srcFolder, '**', '*.m'));
    
    if isempty(files)
        disp('在源文件夹及其子文件夹中未找到任何 .m 文件。');
        return;
    end
    
    copiedCount = 0;
    
    % 8. 遍历并复制文件
    for i = 1:length(files)
        if files(i).isdir
            continue;
        end
        
        % A. 获取源文件的完整路径
        sourceFile = fullfile(files(i).folder, files(i).name);
        
        % B. 从文件的完整目录中，移除源根目录的部分，得到相对路径
        relativePath = erase(files(i).folder, srcFolder);
        
        % C. 构建目标子文件夹的路径
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

    % 9. 结束报告
    disp('--------------------');
    disp('操作完成！');
    disp(['总共复制了 ', num2str(copiedCount), ' 个 .m 文件。']);
    disp(['目标根文件夹: ', finalDestRoot]);

end