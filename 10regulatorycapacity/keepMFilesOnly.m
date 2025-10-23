function keepMFilesOnly()
% 保留所选文件夹（含子文件夹）中的 .m 文件，删除其余所有文件。
% 使用方法：保存为 keepMFilesOnly.m 后，在命令行运行 keepMFilesOnly

    %---------------- 可调参数（安全优先） ----------------%
    forceDelete      = false;  % = true 时不再弹出确认对话框，直接删除（谨慎）
    removeEmptyDirs  = false;  % = true 时在删除文件后再尝试清理空文件夹
    useRecycleBin    = true;   % = true 时启用回收站(Windows)；false 为永久删除
    previewListMax   = 30;     % 预览时最多展示多少条将被删除的文件路径
    %-----------------------------------------------------%

    % 选择根文件夹
    rootDir = uigetdir(pwd, '选择要清理的根文件夹（仅保留 .m 文件）');
    if isequal(rootDir,0)
        disp('操作已取消。');
        return;
    end

    % 删除行为设置：是否使用回收站（仅 Windows 有效）
    if useRecycleBin
        recycle('on');
    else
        recycle('off');
    end

    % 列出所有文件（包含子目录）
    listing = dir(fullfile(rootDir, '**', '*'));
    if isempty(listing)
        fprintf('目标文件夹为空：%s\n', rootDir);
        return;
    end

    % 仅保留“文件”（去掉文件夹项）
    listing = listing(~[listing.isdir]);

    % 生成完整路径
    fullpaths = arrayfun(@(f) fullfile(f.folder, f.name), listing, 'UniformOutput', false);

    % 提取扩展名并判断是否为 .m（不区分大小写）
    [~, ~, exts] = cellfun(@fileparts, fullpaths, 'UniformOutput', false);
    isMfile = cellfun(@(e) strcmpi(e, '.m'), exts);

    toKeep    = fullpaths(isMfile);
    toDelete  = fullpaths(~isMfile);

    fprintf('扫描完成：共发现文件 %d 个；其中 .m 文件 %d 个，将删除非 .m 文件 %d 个。\n', ...
        numel(fullpaths), numel(toKeep), numel(toDelete));

    if isempty(toDelete)
        disp('没有需要删除的文件，一切就绪。');
        return;
    end

    % 预览与确认
    if ~forceDelete
        previewCount = min(numel(toDelete), previewListMax);
        disp('以下为将要删除的部分文件预览：');
        disp('-------------------------------------------');
        for i = 1:previewCount
            fprintf('%s\n', toDelete{i});
        end
        if numel(toDelete) > previewCount
            fprintf('... 以及另外 %d 个文件\n', numel(toDelete) - previewCount);
        end
        disp('-------------------------------------------');

        choice = questdlg( ...
            sprintf('将删除 %d 个非 .m 文件。是否继续？', numel(toDelete)), ...
            '确认删除', ...
            '继续', '取消', '继续');

        if ~strcmp(choice, '继续')
            disp('操作已取消。');
            return;
        end
    end

    % 执行删除
    failed = {};
    for i = 1:numel(toDelete)
        fp = toDelete{i};
        try
            delete(fp); % 遵循 recycle 设置
        catch ME
            warning('无法删除：%s\n原因：%s', fp, ME.message);
            failed{end+1,1} = fp; %#ok<AGROW>
        end
    end

    deletedCount = numel(toDelete) - numel(failed);
    fprintf('删除完成：成功 %d 个，失败 %d 个。\n', deletedCount, numel(failed));

    % （可选）清理空文件夹：自底向上尝试删除
    if removeEmptyDirs
        % 获取所有目录（包括子目录）
        allDirs = dir(fullfile(rootDir, '**'));
        allDirs = allDirs([allDirs.isdir]);
        % 排除 . 和 ..
        allDirs = allDirs(~ismember({allDirs.name}, {'.','..'}));
        % 目录按“深度”从深到浅排序（优先删除更深层空目录）
        dirPaths = arrayfun(@(d) fullfile(d.folder, d.name), allDirs, 'UniformOutput', false);
        depths   = cellfun(@(p) numel(strfind(p, filesep)), dirPaths);
        [~, idx] = sort(depths, 'descend');
        dirPaths = dirPaths(idx);

        removedDirs = 0;
        for i = 1:numel(dirPaths)
            p = dirPaths{i};
            try
                % 仅删除空目录（不使用 's' 递归标志）
                rmdir(p);
                removedDirs = removedDirs + 1;
            catch
                % 非空或无权限，忽略
            end
        end
        fprintf('空文件夹清理完成：共删除空目录 %d 个。\n', removedDirs);
    end

    if ~isempty(failed)
        fprintf('\n以下文件删除失败（可能被占用或无权限）：\n');
        for i = 1:numel(failed)
            fprintf('%s\n', failed{i});
        end
    end

    disp('任务完成。建议检查关键目录是否符合预期。');
end
