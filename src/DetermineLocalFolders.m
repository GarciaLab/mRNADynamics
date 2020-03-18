function [rawDataPath, ProcPath, DropboxFolder, MS2CodePath, PreProcPath,...
    configValues, movieDatabasePath, movieDatabaseFolder, movieDatabase] = DetermineLocalFolders(varargin)
    %Variable arguments - always have Prefix first if doing this
    %   Prefix - sets Dropbox folder to be from specified prefix
    %   'LocalMovieDatabase', LocalDropboxFolder - loads a MovieDatabase
    %   file from a specified alternate dropbox folder string, as defined
    %   in ComputerFolders
    useLocalMovieDatabase = false; %Use default MovieDatabase by default
    if ~isempty(varargin) 
        Prefix = varargin{1}; %optionally return a different dropbox folder from the default with respect to Prefix
        if length(varargin)>1
            otherOptions = varargin(2:end);
            for j = 1:length(otherOptions)
                %Option to load a local dropbox folder's MovieDatabase
                %instead of default
                if strcmpi(otherOptions{j},'LocalMovieDatabase')
                    useLocalMovieDatabase = true;
                    LocalDropboxFolderString = otherOptions{j+1};
                end
            end
        end
    end
    optionalResults = ''; %No optional results by default

    CONFIG_CSV_PATH = 'ComputerFolders.csv';

    configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');
    
    DropboxFolder = getConfigValue(configValues, 'DropboxFolder');
    MS2CodePath = getConfigValue(configValues, 'MS2CodePath');      
    try
        DataRoot = getConfigValue(configValues, 'DataRoot');
        rawDataPath = [DataRoot,filesep,'RawDynamicsData'];
        PreProcPath = [DataRoot,filesep, 'PreProcessedData'];
        ProcPath = [DataRoot, filesep,'ProcessedData'];
    catch  
        rawDataPath = getConfigValue(configValues, 'SourcePath');
        ProcPath = getConfigValue(configValues, 'FISHPath');
        PreProcPath = getConfigValue(configValues, 'PreProcPath');
    end
    
    if useLocalMovieDatabase
        LocalDropboxFolder = getConfigValue(configValues, LocalDropboxFolderString);
        movieDatabaseFolder = LocalDropboxFolder;
    else
        movieDatabaseFolder = DropboxFolder;
    end
    movieDatabasePath = [movieDatabaseFolder,'\MovieDatabase.csv'];
    movieDatabase = csv2cell(movieDatabasePath, 'fromfile');
    
    if isempty(varargin) || isempty(varargin{1})
    %     warning('No Prefix specified. Using default Dropbox folder')
        return
    end
    PREFIX_SEPARATOR = '[\\\\/-]';
    PREFIX_REGEX = ['^.{10}', PREFIX_SEPARATOR, '.*$'];

    %% We need to look for the dropbox folder specified by the provided prefix
  
    if isempty(regexp(Prefix, PREFIX_REGEX, 'once'))
    error('Prefix %s does not match "yyyy-mm-dd[/\\-]name". Please change it accordingly.', Prefix)
    % any 10 characters will work, not only yyyy-mm-dd,
    % but we enforce this in the error msg for simplicity
    end

    if ~isempty(optionalResults)
        dropboxFolderName = getDropboxFolderFromMovieDatabase(movieDatabase, Prefix, PREFIX_SEPARATOR, optionalResults);
    else
        dropboxFolderName = getDropboxFolderFromMovieDatabase(movieDatabase, Prefix, PREFIX_SEPARATOR);
    end
    rootFolderName = getRootFolderFromMovieDatabase(movieDatabase, Prefix, PREFIX_SEPARATOR);
    % if user indicated a RootFolder in movie database, reassign paths
    % accordingly
    if ~strcmpi(rootFolderName,'noFolder')
        if strcmpi(rootFolderName,'default') || isempty(rootFolderName)
            rootFolderName = 'DataRoot';
        end
        DataRoot = getConfigValue(configValues, rootFolderName);        
        rawDataPath = [DataRoot '/RawDynamicsData'];
        PreProcPath = [DataRoot '/PreProcessedData'];
        ProcPath = [DataRoot '/ProcessedData'];
    end
    %We can use the string "Default" or "DropboxFolder" for the default Dropbox folder
    if strcmpi(dropboxFolderName,'default')
        dropboxFolderName = 'DropboxFolder';
    end

    DropboxFolder = getConfigValue(configValues, dropboxFolderName);

end
