% This function creates all the required folders to run the MS2 movie
% analysis pipeline. 
%
% OUTPUT
% txt: A cell containing information within ComputerFolders.xlsx 
%
% Author (contact): Hernan Garcia (hgarcia@berkeley.edu)
function txt = InstallmRNADynamics
  warning('off','MATLAB:MKDIR:DirectoryExists')

  ensureRightFolder()
  mRNADynamicsPath = pwd

  cd('..')
  ROOT_PATH = pwd
  cd(mRNADynamicsPath)

  disp(['mRNADynamics absolute path is ', mRNADynamicsPath])
  disp(['root directory absolute path is ', ROOT_PATH])
  
  PREPROCESSED_DATA_PATH = createDataSubDir('PreProcessedData')
  PROCESSED_DATA_PATH = createDataSubDir('ProcessedData')
  RAW_DYNAMICS_DATA_PATH = createDataSubDir('RawDynamicsData') %(old RawData folder)
  DYNAMICS_RESULTS_PATH = createDataSubDir('DynamicsResults') %(old DropboxFolder)
  MS2CODE_PATH = mRNADynamicsPath

  createComputerFoldersConfig()
      
  createMovieDatabaseFile()

  createStartupFile()

  setJavaHomeEnvironmentVariable()

  msgbox('Run "startup" from the command line or restart Matlab to finish the installation')
          
  warning('on','MATLAB:MKDIR:DirectoryExists')

  %%
  %% nested sub-functions
  %%

  function subDirPath = createDataSubDir(subDirName)
    dataPath = [ROOT_PATH, filesep, 'Data']
    if ~exist(dataPath)
      mkdir(dataPath)
    end

    subDirPath = [dataPath, filesep, subDirName]
    mkdir(subDirPath)
  end

  function ensureRightFolder()
    while ~strfind(pwd, 'mRNADynamics')
      warning('This script must be run from the ''mRNADynamics'' directory.')
      cd(uigetdir);
    end
  end

  function copied = copyFileIfNotExists(fromPath, toPath)
    if ~exist(toPath)
        copyfile(fromPath, toPath)
        copied = true
    else
        warning([toPath, ' already exists. Not overwriting.'])
        copied = false
    end
  end

  function createMovieDatabaseFile()
    movieDatabaseInstallationFile = ['InstallationFiles', filesep, 'InstallMovieDatabase.xlsx']
    movieDatabaseFile = [DYNAMICS_RESULTS_PATH, filesep, 'MovieDatabase.xlsx']

    copyFileIfNotExists(movieDatabaseInstallationFile, movieDatabaseFile)
  end

  function createComputerFoldersConfig()
    computerFoldersInstallationFile = ['InstallationFiles', filesep, 'InstallComputerFolders.xlsx']
    computerFoldersFile = [ROOT_PATH, filesep, 'ComputerFolders.xlsx']

    blankConfigCreated = copyFileIfNotExists(computerFoldersInstallationFile, computerFoldersFile)
    
    if ispc && blankConfigCreated
      %Save the XLS file - only windows supported for now
      writeConfigFile(computerFoldersFile)
    else
        if blankConfigCreated
          disp('Warning: Macs and Linux cannot generate the XLS files.')
          disp('(1) Re-run using "txt=InstallmRNADynamics".')
          disp('(2) Type "open txt".')
          disp('(3) Copy and paste the data into a new file in Excel.')
          disp('(4) Get rid of all '' in the file.')
          disp('(5) Save as "ComputerFolders.xlsx" in folder "LivemRNAFISH".')
        end
    end
  end

  function writeConfigFile(filePath)
    [num,txt] = xlsread(filePath);

    txt{1, end+1} = getComputerName();
    txt{2, end} = getUserName();
    txt{3, end} = RAW_DYNAMICS_DATA_PATH;
    txt{4, end} = PREPROCESSED_DATA_PATH;
    txt{5, end} = PROCESSED_DATA_PATH;
    txt{6, end} = DYNAMICS_RESULTS_PATH;
    txt{7, end} = MS2CODE_PATH;
    
    xlswrite(filePath, txt);
  end
  
  function createStartupFile()
    % Add the right folders to the path.
    % This will be done as a startup file in the user's folder

    DependenciesFolder = [mRNADynamicsPath, filesep, 'dependencies']
    DeprecatedFolder = [mRNADynamicsPath, filesep, 'deprecated']
    TestClassesFolder = [mRNADynamicsPath, filesep, 'testClasses']
    TrackingFolder = [mRNADynamicsPath, filesep, 'Tracking'];
    SubfunctionsFolder = [TrackingFolder, filesep, 'subfunctions'];
    LineageCodeFolder = [mRNADynamicsPath, filesep, 'LineageCode'];

    Output{1} = ['path(''', PREPROCESSED_DATA_PATH, ''',path);'];
    Output{2} = ['addpath(genpath(''', DependenciesFolder, '''))'];
    Output{3} = ['path(''', ROOT_PATH, ''',path);'];
    Output{4} = ['path(''', mRNADynamicsPath, ''',path);'];
    Output{5} = ['path(''', TrackingFolder, ''',path);'];
    Output{6} = ['path(''', LineageCodeFolder, ''',path);'];
    Output{7} = ['path(''', SubfunctionsFolder, ''',path);'];
    Output{8} = ['addpath(genpath(''', DeprecatedFolder, '''))'];
    Output{9} = ['path(''', DYNAMICS_RESULTS_PATH, ''',path);'];
    Output{10} = ['addpath(genpath(''', TestClassesFolder, '''))'];

    writeStartupFile(Output);
  end

  function setJavaHomeEnvironmentVariable()
    %Switch Matlab over to the Java release from the FIJI packaged with this
    %repository
    if ispc
        [status,~] = dos(['setx MATLAB_JAVA ',ROOT_PATH,'\mRNADynamics\Fiji.app\java\win64\jdk1.8.0_66\jre']);
        if status
            warning('Something went wrong setting Java environment variable. Talk to Armando.')
        end
    else    
        warning(['Please note that Weka integration is not supported outside of Windows. ',
                'Talk to Armando if you need this.'])
    end
  end
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
  computerName = lower(computerName);
end

function writeStartupFile(contents)
  disp(['Will create startup.m with this contents: ', contents])

  %Create the startup.m file
  StartUpPath = userpath; 
  if isempty(userpath)
     display('Path for this specific user was not found. Please locate it in your documents folder.')
     display('In Windows, got to "My Documents\Matlab"')
     StartUpPath = uigetdir;
     userpath(StartUpPath);
  end

  disp([StartUpPath(1:end-1),filesep,'startup.m']);
  fid = fopen([StartUpPath(1:end-1),filesep,'startup.m'], 'w');
  errmsg = '';
  if fid < 0
      fid = fopen([StartUpPath,filesep,'startup.m'], 'w');
  end
  while fid < 0 
     disp(errmsg);
     disp('Please find your user''s Matlab folder. Maybe "My Documents\Matlab" on Windows?')
     d = uigetdir;
     [fid,errmsg] = fopen(d);
  end

  for i=1:length(contents)
      fprintf(fid, '%s \n', contents{i});
  end
  fclose(fid);
end
