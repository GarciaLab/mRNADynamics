function configContents = InstallmRNADynamics(varargin)
% DESCRIPTION
% This function creates all the required folders to run the MS2 movie
% analysis pipeline for new installations. For existing installations, it
% will update the MATLAB startup.m script.
%
% PARAMETERS
% None
%
% OPTIONS
% 'updateStartup': Will only update MATLAB's startup.m script. This is a 
%                  redundant option, as this function will automatically
%                  update startup.m whenever it detects an existing
%                  installation
%
% OUTPUT
% configContents: A cell containing information within ComputerFolders configuration
%
% Author (contact): Hernan Garcia (hgarcia@berkeley.edu)
% Last updated: 4/17/2020 by Meghan Turner (meghan_turner@berkeley.edu)

warning('off','MATLAB:MKDIR:DirectoryExists')

addpath('./utilities');

ensureRightFolder();
cd('..');
MRNA_DYNAMICS_PATH = toSafeWindowsString(pwd);

cd('..');
ROOT_PATH = toSafeWindowsString(pwd);
cd(MRNA_DYNAMICS_PATH);

disp(['mRNADynamics absolute path is ', MRNA_DYNAMICS_PATH]);
disp(['root directory absolute path is ', ROOT_PATH]);

% default directory locations
DATA_ROOT_PATH = createDirInRoot('Data');
PREPROCESSED_DATA_PATH = createDataSubDir('PreProcessedData');
PROCESSED_DATA_PATH = createDataSubDir('ProcessedData');
RAW_DYNAMICS_DATA_PATH = createDataSubDir('RawDynamicsData'); %(old RawData folder)
DYNAMICS_RESULTS_PATH = createDataSubDir('DynamicsResults'); %(old DropboxFolder)
MOVIE_DATABASE_PATH =  [DYNAMICS_RESULTS_PATH, filesep, 'MovieDatabase.csv'];
MS2CODE_PATH = [MRNA_DYNAMICS_PATH, filesep, 'src'];
TEST_PATH =  createDirInRoot('ExpectedData');

COMPUTER_FOLDERS_PATH = [ROOT_PATH, filesep, 'ComputerFolders.csv'];

existingInstallationDetected = exist(COMPUTER_FOLDERS_PATH, 'file') &&...
        contains(COMPUTER_FOLDERS_PATH, '.csv') && ...
        exist(MOVIE_DATABASE_PATH, 'file') &&...
        exist(MOVIE_DATABASE_PATH, 'file');

shouldOnlyMakeStartupFile = ~isempty(varargin) &&...
        contains(varargin, 'updateStartup', 'IgnoreCase', true);
    
if shouldOnlyMakeStartupFile
    
    disp('Updating MATLAB startup script');
    createStartupFile();
    msgbox('Run "startup" from the command line or restart Matlab to finish the installation');

elseif existingInstallationDetected
    disp('Existing installation detected.')
    disp('Only updataing MATALB startup script')
    createStartupFile();
    msgbox('Run "startup" from the command line or restart Matlab to finish the installation');
    
elseif exist(COMPUTER_FOLDERS_PATH, 'file') &&...
        contains(COMPUTER_FOLDERS_PATH, '.xls')
    
    
    disp(['Existing installation detected (', COMPUTER_FOLDERS_PATH,...
        ' exists). Will try to update configurations to CSV format.']);
    
    configContents = migrateComputerFoldersToCSV();
    
    migrateMovieDatabaseToCSV(getConfigValue(configContents, 'DropboxFolder'));
    
    disp('Existing installation files successfully updated to CSV format.');
else
    
    disp('New installation detected. Will create required configurations and folder structures.');
    
    configContents = createFoldersConfig();
    
    createMovieDatabaseFile();
    
    createStartupFile();
            
end

warning('on','MATLAB:MKDIR:DirectoryExists');

%%
%% nested sub-functions
%%

    function subDirPath = createDataSubDir(subDirName)
        dataPath = createDirInRoot('Data');
        
        subDirPath = [dataPath, filesep, subDirName];
        mkdir(subDirPath);
    end

    function dirPath = createDirInRoot(folderName)
        dirPath = [ROOT_PATH, filesep, folderName];
        if ~exist(dirPath, 'file')
            mkdir(dirPath);
        end
    end

    function ensureRightFolder()
        while ~contains(pwd, 'src')
            cd(fileparts(mfilename('fullpath')))
        end
    end

    function createMovieDatabaseFile()
        if exist(MOVIE_DATABASE_PATH, 'file')
            warning([MOVIE_DATABASE_PATH, ' already exists. Not overwriting.']);
            return
        end
        
        contents = {
            'Date', 'ExperimentType', 'ExperimentAxis', 'CoatProtein', ...
            'StemLoop', 'APResolution','DVResolution', 'Channel1', 'Channel2',...
            'Channel3', 'Objective', 'Power Channel 1 (mW)',...
            'Power Channel 2 (mW)', 'Power Channel 3 (mW)', 'RootFolder',...
            'DataFolder', 'DropboxFolder','Comments', 'nc9', 'nc10', 'nc11',...
            'nc12', 'nc13', 'nc14', 'CF'
            };
        
        cell2csv(MOVIE_DATABASE_PATH, contents);
    end

    function contents = createFoldersConfig()
        if exist(COMPUTER_FOLDERS_PATH, 'file')
            warning([COMPUTER_FOLDERS_PATH, ' already exists. Not overwriting.']);
        else
            contents = {
                'Computer Name', getComputerName();
                'User Name', getUserName();
                'SourcePath', RAW_DYNAMICS_DATA_PATH;
                'PreProcPath', PREPROCESSED_DATA_PATH;
                'FISHPath', PROCESSED_DATA_PATH;
                'DropboxFolder', DYNAMICS_RESULTS_PATH;
                'DataRoot', DATA_ROOT_PATH;
                'MS2CodePath', MS2CODE_PATH;
                'TestPath', TEST_PATH
                };
            
            cell2csv(COMPUTER_FOLDERS_PATH, contents);
        end
        
        contents = csv2cell(COMPUTER_FOLDERS_PATH, 'fromfile');
    end

    function contents = migrateComputerFoldersToCSV
        ComputerFoldersCSVPath = [ROOT_PATH, filesep, 'ComputerFolders.csv'];
        ComputerFoldersXLSPath = [ROOT_PATH, filesep, 'ComputerFolders.xlsx'];
        contents = migrateFileToCSV(ComputerFoldersXLSPath, ComputerFoldersCSVPath, {''});
    end

    function migrateMovieDatabaseToCSV(movieDatabaseDirectory)
        MovieDatabaseCSVPath = [movieDatabaseDirectory, filesep, 'MovieDatabase.csv'];
        MovieDatabaseXLSPath = [movieDatabaseDirectory, filesep, 'MovieDatabase.xlsx'];
        migrateFileToCSV(MovieDatabaseXLSPath, MovieDatabaseCSVPath, {'date'});
    end

    function contents = migrateFileToCSV(fromXLSPath, toCSVPath, types)
        if exist(toCSVPath, 'file')
            warning([toCSVPath, ' already exists. Not overwriting.']);
        else
            if exist(fromXLSPath, 'file')
                disp([fromXLSPath, ' exists in XLS format. Will try to convert it to CSV format.']);
                [~, ~, rawData] = xlsread(fromXLSPath);
                [excelRows, excelColumns] = size(rawData);
                disp(['Excel size is ', num2str(excelRows), ' rows x ', num2str(excelColumns), ' columns'])
                
                if(excelColumns > 26)
                    warning('Warning, number of columns exceeds maximum supported. Will truncate at column 26 (Z)')
                    [~, ~, rawData] = xlsread(fromXLSPath, strcat('A1:Z', num2str(excelRows)));
                end
                disp('XLS file read OK. Starting conversion to CSV.');
                cell2csv(toCSVPath, rawData, ',', types);
                disp(['CSV File created: ', toCSVPath]);
            else
                error([fromXLSPath, ' does not exist. Cannot migrate to CSV format.']);
            end
        end
        
        contents = csv2cell(toCSVPath, 'fromfile');
    end

    function createStartupFile()
        % Add the right folders to the path.
        % This will be done as a startup file in the user's folder
        
        
        
        srcFolder = [MRNA_DYNAMICS_PATH, filesep, 'src'];
        libFolder = [MRNA_DYNAMICS_PATH, filesep, 'lib'];
        DependenciesFolder = [MRNA_DYNAMICS_PATH, filesep, 'lib/dependencies'];
        testFolder = [MRNA_DYNAMICS_PATH, filesep, 'test'];
        
        % matlab paths
        Output{1} = ['addpath(genpath(''', libFolder, filesep, 'Fiji.app', filesep, 'scripts''));'];
        Output{2} = ['path(''', PREPROCESSED_DATA_PATH, ''',path);'];
        Output{3} = ['path(''', ROOT_PATH, ''',path);'];
        Output{4} = ['path(''', MRNA_DYNAMICS_PATH, ''',path);'];
        Output{5} = ['addpath(genpath(''', srcFolder, '''));'];
        Output{6} = ['path(''', DYNAMICS_RESULTS_PATH, ''',path);'];
        Output{7} = ['addpath(genpath(''', testFolder, '''));'];
        Output{8} = ['addpath(genpath(''', DependenciesFolder, '''));'];
        Output{10} = ['cd(''', MRNA_DYNAMICS_PATH, ''');'];
        Output{11} = 'disp(''Startup script executed.'');';
        
        
        updateStartupFile = existingInstallationDetected || shouldOnlyMakeStartupFile;
        writeStartupFile(Output,updateStartupFile);
                
    end
end

function safeString = toSafeWindowsString(aPath)
% replaces any '\' with '/', since MATLAB can handle '/' in Windows Paths as well,
% so we don't need to handle linux or windows paths differently
safeString = strrep(aPath, '\', '/');
end

function userName = getUserName
%Find out the username
if isunix
    [~, userName] = system('whoami'); % exists on every unix that I know of
    % on my mac, isunix == 1
elseif ispc
    [~, userName] = system('echo %USERDOMAIN%\%USERNAME%');
    % (Not as familiar with windows,
    % found it on the net elsewhere, you might want to verify)
end

userName = strrep(userName, newline,'');
userName = toSafeWindowsString(userName);
end

function computerName = getComputerName
%Find out the computer name. This will probably fail for nonwindows but
%that's easily fixable if required
[ret, computerName] = system('hostname');
if ret ~= 0
    if ispc
        computerName = getenv('COMPUTERNAME');
    else
        computerName = getenv('HOSTNAME');
    end
end

computerName = strrep(lower(computerName), newline,'');
end

function writeStartupFile(contents,updateStartupFile)

%Create the startup.m file
StartUpPath = userpath;
if isempty(userpath)
    disp('Path for this specific user was not found. Please locate it in your documents folder.');
    disp('In Windows, got to "My Documents\Matlab"');
    StartUpPath = uigetdir;
    userpath(StartUpPath);
end

startupFile = [StartUpPath, '/startup.m']; 

if exist(startupFile, 'file') && ~updateStartupFile
    warning([startupFile, ' already exists. Not overwriting.']);
    return
end

disp('Will create startup.m with these contents: ');
disp(contents);



disp(startupFile);
fid = fopen(startupFile, 'w');
errmsg = '';
if fid < 0
    fid = fopen(startupFile, 'w');
end
while fid < 0
    disp(errmsg);
    disp('Please find your user''s Matlab folder. Maybe "My Documents\Matlab" on Windows?')
    d = uigetdir;
    [fid,errmsg] = fopen(d, 'w');
end

for i=1:length(contents)
    fprintf(fid, '%s \n', contents{i});
end
fclose(fid);

msgbox('Run "startup" from the command line to finish the installation.');


end
