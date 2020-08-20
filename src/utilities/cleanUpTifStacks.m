%Safely removes extra tif stacks created by the pipeline. 

%Talk to AR if you have issues.
function cleanUpTifStacks(varargin)

%optionally input a prefix to choose a project directory
if isempty(varargin)
    [SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath, configValues, movieDatabasePath]=...
        DetermineLocalFolders;
else
    Prefix = varargin{1};
    [SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath, configValues, movieDatabasePath]=...
        DetermineLocalFolders(Prefix);
end

preDir = dir(PreProcPath);
for i = 3:length(preDir)
    
    prefix = preDir(i).name;
    stacksPath = [PreProcPath, filesep, prefix, filesep, 'stacks'] ;
    stacksDir = dir([stacksPath,filesep, '*.tif']);
    if ~isempty(stacksDir)
        
        nFiles = length(stacksDir);
        if nFiles > 0
            disp(['Deleting ', num2str(nFiles - 1), ' exported .tifs from ', stacksPath, '...']);
        else
            disp(['Deleting ', num2str(nFiles), ' exported .tifs from ', stacksPath, '...']);
        end
        for k = 1:nFiles
            stackName = stacksDir(k).name;
            fullPath = [stacksPath, filesep, stackName];
%             delete(fullPath)
            java.io.File(fullPath).delete()
        end
        rmdir(stacksPath);
        
    end
end


end
