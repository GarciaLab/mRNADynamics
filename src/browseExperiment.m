function [rawExportedDirectory,...
    processedDirectory,...
    resultsDirectory] = browseExperiment(Prefix, varargin)

displayOutput = false;

if ~isempty(varargin)
    displayOutput = varargin{2};
end

[~, ProcPath, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix);

swapColumns = @(T) movevars(T, 'folder', 'After', 'ext');




rawExportedFolder = [PreProcPath, filesep, Prefix, filesep];
% disp('Raw exported files:')
rawExportedDirectory = swapColumns(dirtab(rawExportedFolder));
% disp({rawExportedDirectory.name}');

processedFolder = [ProcPath, filesep, Prefix, '_', filesep];
% disp('Processed images:')
processedDirectory = swapColumns(dirtab(processedFolder));
% disp({processedDirectory.name}');

resultsFolder = [DropboxFolder, filesep, Prefix, filesep];
% disp('Analysis results:')
resultsDirectory = swapColumns(dirtab(resultsFolder));
% disp({resultsDirectory.name}');






%     [~,~,DropboxFolder,~, PreProcPath,...
%         ~, ~, ~,Channel1,Channel2,~,...
%         Channel3, ~, movieDatabaseFolder, movieDatabase]...
%         = readMovieDatabase(Prefix);



end