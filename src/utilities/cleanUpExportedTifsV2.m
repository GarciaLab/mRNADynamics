%Safely removes extra tifs exported by LASX. This will leave exactly one tif/jpg per
%project.

%Talk to AR if you have issues.
function cleanUpExportedTifsV2(varargin)

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
        
        %% main folder
        prefix = dateDir(j).name;
        projPath = [datePath, filesep, prefix] ;
        projDir = dir([projPath, filesep, '*_Series*.tif']);
        projDirJPG = dir([projPath, filesep, '*_Series*.jpg']);
        projDirbfMemo = dir([projPath, filesep, '*.bfmemo']);
        projDir = [projDir, projDirJPG, projDirbfMemo]; %#ok<AGROW>
        nFiles = length(projDir);
        if nFiles > 0
            disp(['Deleting ', num2str(nFiles - 1), ' exported .tifs from ', projPath, '...']);
        else
            disp(['Deleting ', num2str(nFiles), ' exported .tifs from ', projPath, '...']);
        end
        for k = 1:nFiles
            projName = projDir(k).name;
            fullPath = [projPath, filesep, projName];
            delete(fullPath)
            %             java.io.File(fullPath).delete();
        end
        %% now for full embryo
        projDir = dir([projPath, filesep, 'FullEmbryo',filesep, '*.tif']);
        projDirJPG = dir([projPath,filesep, 'FullEmbryo', filesep,'*.jpg']);
        projDir = [projDir, projDirJPG]; %#ok<AGROW>
        nFiles = length(projDir);
        if nFiles > 0
            disp(['Deleting ', num2str(nFiles - 1), ' exported full embryo .tifs from ', projPath, '...']);
        else
            disp(['Deleting ', num2str(nFiles), ' exported full embryo .tifs from ', projPath, '...']);
        end
        for k = 1:nFiles
            projName = projDir(k).name;
            fullPath = [projPath, filesep, 'FullEmbryo', filesep, projName];
            delete(fullPath)
            %             java.io.File(fullPath).delete();
        end
        if exist([projPath,filesep, 'FullEmbryo', filesep, 'MetaData'], 'dir')
            rmdir([projPath,filesep, 'FullEmbryo', filesep, 'MetaData'], 's');
        end
        
        if exist([projPath,filesep, 'MetaData'], 'dir')
            disp(['Deleting metadata from ', projPath,filesep, 'MetaData', '...']);
            rmdir([projPath,filesep,  'MetaData'], 's');
        end
        
    end   %date loop
end %raw loop

end