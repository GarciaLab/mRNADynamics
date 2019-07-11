function [rawDataPath,ProcPath,DropboxFolder,MS2CodePath, PreProcPath,...
    rawDataFolder, Prefix, ExperimentType,Channel1,Channel2,OutputFolder,...
    Channel3, spotChannels, movieDatabaseFolder]...
= readMovieDatabase(Prefix, varargin)
    
    optionalResults = '';
    if ~isempty(varargin)
        optionalResults = varargin{1};
    end
        
    %Figure out the initial folders. We'll update the Drobpox one later on in the code.

    % [SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath, configValues, movieDatabasePath]=...
    %     DetermineLocalFolders;
    [rawDataPath,~,~,~, ~, ~, movieDatabasePath, movieDatabaseFolder]=...
        DetermineLocalFolders;

    %Get the Prefix if is not already present
    if isempty(Prefix)
        rawDataFolder = uigetdir(rawDataPath,'Select folder with data');


        %Get the information from the last two folders in the structure
        SlashPositions = strfind(rawDataFolder,filesep);
        Prefix = [rawDataFolder((SlashPositions(end-1)+1):(SlashPositions(end)-1)),'-',...
            rawDataFolder((SlashPositions(end)+1):(end))];
    else 
        %Obtains the subfolder using the Prefix (replaces '-' with '/' after the date,
        %knowing it takes 10 characters)
        Subfolder = [Prefix(1:10),filesep,Prefix(12:length(Prefix))];
        rawDataFolder = strcat(rawDataPath,filesep,Subfolder);
    end

    %What type of experiment are we dealing with? Get this out of MovieDatabase
    movieDatabase = csv2cell(movieDatabasePath, 'fromfile');
    movieDatabaseHeaderRow = movieDatabase(1, :);
    ExperimentTypeColumn = findColumnIndex(movieDatabaseHeaderRow, 'ExperimentType');
    Channel1Column = findColumnIndex(movieDatabaseHeaderRow, 'Channel1');
    Channel2Column = findColumnIndex(movieDatabaseHeaderRow, 'Channel2');
    Channel3Column = findColumnIndex(movieDatabaseHeaderRow, 'Channel3');


    [~, PrefixRow] = getDropboxFolderFromMovieDatabase(movieDatabase, Prefix, '[\\\\/-]', optionalResults);

    ExperimentType = movieDatabase(PrefixRow, ExperimentTypeColumn);
    Channel1 = movieDatabase(PrefixRow, Channel1Column);
    Channel2 = movieDatabase(PrefixRow, Channel2Column);
    Channel3 = movieDatabase(PrefixRow, Channel3Column);
    
    if isempty(Channel3)
        Channel3 = {'DoesNotExist'};
    end
    
[rawDataPath,ProcPath,DropboxFolder,MS2CodePath, PreProcPath, ~, ~]=...
    DetermineLocalFolders(Prefix, optionalResults);


    %Set the destination folders
    OutputFolder = [DropboxFolder, filesep, Prefix];
    
    %Determine the spot channel(s)
    spotChannels = getCoatChannel(Channel1, Channel2, Channel3);

end
