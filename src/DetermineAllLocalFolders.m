function [SourcePath, FISHPath, DefaultDropboxFolder,...
    DropboxFolder, MS2CodePath, PreProcPath, configValues,...
    movieDatabasePath, movieDatabase] = DetermineAllLocalFolders(varargin)

    optionalResults = '';
    
    if length(varargin)>1
        optionalResults = varargin{2};
    end
    
    %Get default Dropbox folder location
    CONFIG_CSV_PATH = 'ComputerFolders.csv';
    configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');    
    DefaultDropboxFolder = getConfigValue(configValues, 'DropboxFolder');
    
    %Now get rest of folder locations
    if(~isempty(varargin))   
        [SourcePath, FISHPath, DropboxFolder, MS2CodePath, PreProcPath, configValues, movieDatabasePath, movieDatabase] = DetermineLocalFolders(varargin{1}, optionalResults);
    else
        [SourcePath, FISHPath, DropboxFolder, MS2CodePath, PreProcPath, configValues, movieDatabasePath, movieDatabase] = DetermineLocalFolders;
    end
