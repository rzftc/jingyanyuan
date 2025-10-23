function combineMFilesToTXT(folderPath, outputFilePath)
%COMBINEMFILESTOTXT 合并指定文件夹及子文件夹下的所有.m文件到单个txt文件
%   输入参数：
%   folderPath - 需要合并的根文件夹路径
%   outputFilePath - 输出文件路径（包含文件名）

% 验证输入文件夹是否存在
if ~isfolder(folderPath)
    error('文件夹不存在：%s', folderPath);
end

% 获取并过滤隐藏文件夹（Windows系统专用方法）
[allFolders, filteredFolders] = getWindowsFilteredFolders(folderPath);

% 初始化文件列表结构体
fileList = struct('name', {}, 'folder', {}, 'date', {},...
                 'bytes', {}, 'isdir', {}, 'datenum', {},...
                 'relativePath', {});

% 收集所有.m文件（包含子文件夹）
fileList = collectMFiles(folderPath, filteredFolders, fileList);

% 创建并写入输出文件
writeCombinedFile(outputFilePath, fileList, folderPath);

fprintf('成功合并 %d 个文件到：%s\n', length(fileList), outputFilePath);
end

% 子函数：获取过滤后的文件夹列表（Windows专用）
function [allFolders, filteredFolders] = getWindowsFilteredFolders(folderPath)
% 获取所有子文件夹路径
allFolders = genpath(folderPath);
allFolders = strsplit(allFolders, ';');
allFolders = allFolders(~cellfun(@isempty, allFolders));

filteredFolders = {};
for i = 1:length(allFolders)
    currentFolder = allFolders{i};
    
    % Windows系统检测隐藏属性
    try
        [status, attrib] = fileattrib(currentFolder);
        if status && ~any(attrib.hidden)
            filteredFolders{end+1} = currentFolder;
        end
    catch
        % 跳过无访问权限的文件夹
    end
end
end

% 子函数：收集.m文件（保持不变）
function fileList = collectMFiles(rootPath, folders, fileList)
for i = 1:length(folders)
    currentFolder = folders{i};
    files = dir(fullfile(currentFolder, '*.m'));
    files = files(~[files.isdir]);
    
    if ~isempty(files)
        for f = 1:length(files)
            relPath = strrep(fullfile(files(f).folder, files(f).name),...
                            [rootPath filesep], '');
            orderedFile = struct(...
                'name', files(f).name,...
                'folder', files(f).folder,...
                'date', files(f).date,...
                'bytes', files(f).bytes,...
                'isdir', files(f).isdir,...
                'datenum', files(f).datenum,...
                'relativePath', relPath);
            fileList(end+1) = orderedFile;
        end
    end
end
end

% 子函数：写入合并文件（保持不变）
function writeCombinedFile(outputFilePath, fileList, rootPath)
fid = fopen(outputFilePath, 'w', 'n', 'UTF-8');
if fid == -1
    error('无法创建输出文件：%s', outputFilePath);
end

try
    [~, idx] = sort({fileList.relativePath});
    sortedFiles = fileList(idx);
    
    for i = 1:length(sortedFiles)
        filePath = fullfile(sortedFiles(i).folder, sortedFiles(i).name);
        fileContent = fileread(filePath);
        
        fprintf(fid, '%%%% %s %%%%\n', sortedFiles(i).relativePath);
        fprintf(fid, '%s\n\n', fileContent);
    end
    fclose(fid);
catch ME
    fclose(fid);
    rethrow(ME);
end
end
