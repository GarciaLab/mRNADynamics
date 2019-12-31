%Safely removes extra tifs exported by LASX. This will leave exactly one tif/jpg per
%project. 

%Talk to AR if you have issues. 
function cleanUpExportedTifs(varargin)

%optionally input a prefix to choose a project directory
if isempty(varargin)
    [SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath, configValues, movieDatabasePath]=...
        DetermineLocalFolders;
else
    Prefix = varargin{1};
    [SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath, configValues, movieDatabasePath]=...
        DetermineLocalFolders(Prefix);
end

rawDir = dir(SourcePath);
for i = 1:length(rawDir)
    datePath = [SourcePath, filesep, rawDir(i).name];
    dateDir = dir(datePath);
    for j = 1:length(dateDir)
        prefix = dateDir(j).name;
        projPath = [datePath, filesep, prefix] ;
        projDir = dir([projPath, filesep, prefix, '_Series*.tif']);
        projDirJPG = dir([projPath, filesep, prefix, '_Series*.jpg']);
        projDir = [projDir, projDirJPG];
        nFiles = length(projDir);
        if nFiles > 0
            disp(['Deleting ', num2str(nFiles - 1), ' exported .tifs from ', projPath, '...']);
        else
            disp(['Deleting ', num2str(nFiles), ' exported .tifs from ', projPath, '...']);
        end
        for k = 2:nFiles
            projName = projDir(k).name;
            fullPath = [projPath, filesep, projName];
            delete(fullPath)
        end        
    end   
end

end