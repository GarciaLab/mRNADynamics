function [rawDataPath, ProcPath, defaultDropboxFolder, MS2CodePath, ...
            PreProcPath, computerFoldersValues, movieDatabasePath,...
            movieDatabaseFolder, movieDatabase, allDropboxFolders, ...
            experimentsFolder] = DetermineLocalFolders(varargin)
        
% [rawDataPath, ProcPath, defaultDropboxFolder, MS2CodePath, PreProcPath, 
%   configValues, movieDatabasePath, movieDatabaseFolder, movieDatabase, 
%   allDropboxFolders, experimentsFolder] = DetermineLocalFolders(varargin)
%
% DESCRIPTION
% Extracts folder paths from the user's ComputerFolders.csv file.
%
% PARAMETERS
% None
%
% OPTIONS
% Prefix: Sets Dropbox folder to be from specified prefix. Must be the
%         first option.
% 'localMovieDatabase', localDropboxFolder: Always have Prefix as first 
%       input parameter if using this option. Loads the MovieDatabase.csv
%       file from the specified alternate dropbox folder string, as defined
%       in ComputerFolders.csv
%
% OUTPUT
% rawDataPath: Location of the RawDynamicsData folder 
% ProcPath: Location of the ProcessedData folder 
% defaultDropboxFolder: Location of the default DynamicsResults folder (as  
%                       specified by 'DropboxFolder' in
%                       ComputerFolders.csv)
% MS2CodePath: Location of the mRNADynamics code
% PreProcPath: Location of the PreProcessedData folder 
% computerFoldersValues: (formerly configValues) cell array of the contents
%                         of ComputerFolders.csv
% movieDatabasePath: Path pointing directly to the MovieDatabase.csv file 
% movieDatabaseFolder: Location of the folder containing MovieDatabase.csv
% movieDatabase: Cell array of the contents of MovieDatabase.csv
% allDropboxFolders: Locations of all results folders defined in
%                    ComputerFolders.csv.
%                    Useful if you want to find multiple MovieDatabase.csv 
%                    or DataStatus.csv files.
% experimentsFolder: Unknown, not used in this function as of 5/13/2020
%
%
% Author (contact): Hernan Garcia (hggarcia@berkeley.edu)
% Created: ????
% Last Updated: 5/13/2020 MT (previously 3/18/2020 JL)
%
% Documented by: Meghan Turner (meghan_turner@berkeley.edu)

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

        
    CONFIG_CSV_PATH =  'ComputerFolders.csv';
    
    computerFoldersValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');
    
    defaultDropboxFolder = getConfigValue(computerFoldersValues, 'DropboxFolder');
    % Find all the other results folders
    allDropboxRows = contains(computerFoldersValues(:,1),'Dropbox') | ... 
                     contains(computerFoldersValues(:,1),'Results');
    allDropboxFolders = computerFoldersValues(allDropboxRows,2);
    
    MS2CodePath = getConfigValue(computerFoldersValues, 'MS2CodePath');      
    try
        DataRoot = getConfigValue(computerFoldersValues, 'DataRoot');
        rawDataPath = [DataRoot,filesep,'RawDynamicsData'];
        PreProcPath = [DataRoot,filesep, 'PreProcessedData'];
        ProcPath = [DataRoot, filesep,'ProcessedData'];
    catch  
        rawDataPath = getConfigValue(computerFoldersValues, 'SourcePath');
        ProcPath = getConfigValue(computerFoldersValues, 'FISHPath');
        PreProcPath = getConfigValue(computerFoldersValues, 'PreProcPath');
    end
    
    % Find local MovieDatabase.csv
    if useLocalMovieDatabase
        LocalDropboxFolder = getConfigValue(computerFoldersValues, LocalDropboxFolderString);
        movieDatabaseFolder = LocalDropboxFolder;
    else
        movieDatabaseFolder = defaultDropboxFolder;
    end
    
    movieDatabasePath = [movieDatabaseFolder,filesep, 'MovieDatabase.csv'];
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
        DataRoot = getConfigValue(computerFoldersValues, rootFolderName);        
        rawDataPath = [DataRoot , filesep, 'RawDynamicsData'];
        PreProcPath = [DataRoot, filesep, 'PreProcessedData'];
        ProcPath = [DataRoot, filesep, 'ProcessedData'];
    end
    %We can use the string "Default" or "DropboxFolder" for the default Dropbox folder
    if strcmpi(dropboxFolderName,'default')
        dropboxFolderName = 'DropboxFolder';
    end

    defaultDropboxFolder = getConfigValue(computerFoldersValues, dropboxFolderName);

end
