% This function creates all the required folders to run the MS2 movie
% analysis pipeline. 
%
% OUTPUT
% configContents: A cell containing information within ComputerFolders configuration
%
% Author (contact): Hernan Garcia (hgarcia@berkeley.edu)
function configContents = InstallmRNADynamics(varargin)
  warning('off','MATLAB:MKDIR:DirectoryExists')

  ensureRightFolder();
  cd('..');
  MRNA_DYNAMICS_PATH = toSafeWindowsString(pwd);

  cd('..');
  ROOT_PATH = toSafeWindowsString(pwd);
  cd(MRNA_DYNAMICS_PATH);

  disp(['mRNADynamics absolute path is ', MRNA_DYNAMICS_PATH]);
  disp(['root directory absolute path is ', ROOT_PATH]);

  % default directory locations
  PREPROCESSED_DATA_PATH = createDataSubDir('PreProcessedData');
  PROCESSED_DATA_PATH = createDataSubDir('ProcessedData');
  RAW_DYNAMICS_DATA_PATH = createDataSubDir('RawDynamicsData'); %(old RawData folder)
  DYNAMICS_RESULTS_PATH = createDataSubDir('DynamicsResults'); %(old DropboxFolder)
  MOVIE_DATABASE_PATH =  [DYNAMICS_RESULTS_PATH, filesep, 'MovieDatabase.csv'];
  MS2CODE_PATH = [MRNA_DYNAMICS_PATH, filesep, 'src'];
  TEST_PATH =  createDirInRoot('ExpectedData');

  COMPUTER_FOLDERS_PATH = [ROOT_PATH, filesep, 'ComputerFolders.csv'];

  if (~isempty(varargin) & strfind('updateStartupScript', varargin))
    disp('Updating MATLAB startup script');
    createStartupFile();
  elseif exist(COMPUTER_FOLDERS_PATH)
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

    msgbox('Run "startup" from the command line or restart Matlab to finish the installation');
  end
          
  warning('on','MATLAB:MKDIR:DirectoryExists');
  startup;

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
    if ~exist(dirPath)
      mkdir(dirPath);
    end
  end

  function ensureRightFolder()
    while ~strfind(pwd, 'src')
      warning('This script must be run from the ''mRNADynamics/src'' directory.');
      cd(uigetdir);
    end
  end

  function createMovieDatabaseFile()
    if exist(MOVIE_DATABASE_PATH)
      warning([MOVIE_DATABASE_PATH, ' already exists. Not overwriting.']);
      return
    end

    contents = {
      'Date', 'ExperimentType', 'ExperimentAxis', 'CoatProtein', 'StemLoop', 'APResolution',...
      'Channel1', 'Channel2', 'Objective', 'Power Channel 1 (mW)', 'Power Channel 2 (mW)', 'DataFolder', 'DropboxFolder',...
      'Comments', 'nc9', 'nc10', 'nc11', 'nc12', 'nc13', 'nc14', 'CF'
    };

    cell2csv(MOVIE_DATABASE_PATH, contents);
  end

  function contents = createFoldersConfig()
    if exist(COMPUTER_FOLDERS_PATH)
      warning([COMPUTER_FOLDERS_PATH, ' already exists. Not overwriting.']);
    else
      contents = {
        'Computer Name', getComputerName();
        'User Name', getUserName();
        'SourcePath', RAW_DYNAMICS_DATA_PATH;
        'PreProcPath', PREPROCESSED_DATA_PATH;
        'FISHPath', PROCESSED_DATA_PATH;
        'DropboxFolder', DYNAMICS_RESULTS_PATH;
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
    if exist(toCSVPath)
      warning([toCSVPath, ' already exists. Not overwriting.']);
    else
      if(exist(fromXLSPath))
        disp([fromXLSPath, ' exists in XLS format. Will try to convert it to CSV format.']);
        [numericData, textData, rawData] = xlsread(fromXLSPath);
        [excelRows, excelColumns] = size(rawData);
        disp(['Excel size is ', num2str(excelRows), ' rows x ', num2str(excelColumns), ' columns'])
        
        if(excelColumns > 26)
            warning('Warning, number of columns exceeds maximum supported. Will truncate at column 26 (Z)')
            [numericData, textData, rawData] = xlsread(fromXLSPath, strcat('A1:Z', num2str(excelRows)));
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

    % directory constants
    Output{9} = ['ROOT_PATH = ''', ROOT_PATH, ''';'];
    Output{10} = ['MRNA_DYNAMICS_PATH = ''', MRNA_DYNAMICS_PATH, ''';'];
    Output{11} = ['PREPROCESSED_DATA_PATH = ''', PREPROCESSED_DATA_PATH, ''';'];
    Output{12} = ['PROCESSED_DATA_PATH = ''', PROCESSED_DATA_PATH, ''';'];
    Output{13} = ['RAW_DYNAMICS_DATA_PATH = ''', RAW_DYNAMICS_DATA_PATH, ''';'];
    Output{14} = ['DYNAMICS_RESULTS_PATH = ''', DYNAMICS_RESULTS_PATH, ''';'];
    Output{15} = ['MS2CODE_PATH = ''', MS2CODE_PATH, ''';'];
    Output{16} = ['MOVIE_DATABASE_PATH = ''', MOVIE_DATABASE_PATH, ''';'];
    Output{17} = ['COMPUTER_FOLDERS_PATH = ''', COMPUTER_FOLDERS_PATH, ''';'];

    Output{18} = ['disp(''Startup script executed.'');'];
    Output{19} = ['clear all;'];
    
    writeStartupFile(Output);
  end
end

function safeString = toSafeWindowsString(aPath)
  % replaces any '\' with '/', since MATLAB hand handle '/' in Windows Paths as well,
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

  userName = strrep(userName, sprintf('\n'),'');
  userName = toSafeWindowsString(userName);
end

function computerName = getComputerName
  %Find out the computer name
  [ret, computerName] = system('hostname');  
  if ret ~= 0,  
     if ispc
        computerName = getenv('COMPUTERNAME');  
     else  
        computerName = getenv('HOSTNAME');  
     end  
  end  

  computerName = strrep(lower(computerName), sprintf('\n'),'');
end

function writeStartupFile(contents)
  disp('Will create startup.m with these contents: ');
  disp(contents);

  %Create the startup.m file
  StartUpPath = userpath; 
  if isempty(userpath)
     disp('Path for this specific user was not found. Please locate it in your documents folder.');
     disp('In Windows, got to "My Documents\Matlab"');
     StartUpPath = uigetdir;
     userpath(StartUpPath);
  end

  disp([StartUpPath(1:end-1), '/startup.m']);
  fid = fopen([StartUpPath(1:end-1), '/startup.m'], 'w');
  errmsg = '';
  if fid < 0
      fid = fopen([StartUpPath, '/startup.m'], 'w');
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
end
