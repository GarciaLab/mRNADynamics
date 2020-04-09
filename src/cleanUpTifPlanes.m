function cleanUpTifPlanes(varargin)
%Safely removes extra tif planes now deprecated in the pipeline. 

%optionally input a prefix to choose a project directory
if isempty(varargin)
    [SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath, configValues, movieDatabasePath]=...
        DetermineLocalFolders;
else
    Prefix = varargin{1};
    [SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath, configValues, movieDatabasePath]=...
        DetermineLocalFolders(Prefix);
end

preDir = dir([PreProcPath]);
for i = 3:length(preDir)
    
    prefix = preDir(i).name;
    tifsPath = [PreProcPath, filesep, prefix, filesep];
    tifsDir = dir([tifsPath,filesep, '*.tif']);
    if ~isempty(tifsDir)
        
        nFiles = length(tifsDir);
        if nFiles > 0
            disp(['Deleting ', num2str(nFiles - 1), ' exported .tifs from ', tifsPath, '...']);
        else
            disp(['Deleting ', num2str(nFiles), ' exported .tifs from ', tifsPath, '...']);
        end
        for k = 1:nFiles
            tifName = tifsDir(k).name;
            fullPath = [tifsPath, filesep, tifName];
%             delete(fullPath)
            java.io.File(fullPath).delete();
        end
    end
end


end