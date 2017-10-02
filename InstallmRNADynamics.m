% This function creates all the required folders to run the MS2 movie
% analysis pipeline. 
%
% OUTPUT
% txt: A cell containing information within ComputerFolders.xlsx 
%
% Author (contact): Hernan Garcia (hgarcia@berkeley.edu)
function txt = InstallmRNADynamics
  warning('off','MATLAB:MKDIR:DirectoryExists')
  
  %make sure we are in the right folder
  while ~strfind(pwd, 'mRNADynamics')
      warning('This script must be run from the ''mRNADynamics'' directory.')
      cd(uigetdir);
  end

  mRNADynamicsPath = pwd

  cd('..')
  rootPath = pwd
  cd(mRNADynamicsPath)

  disp(['mRNADynamics absolute path is ', mRNADynamicsPath])
  disp(['root directory absolute path is ', rootPath])

  % configuration paths
  DATA_PATH = [rootPath, filesep, 'Data']
  PREPROCESSED_DATA_PATH = [DATA_PATH, filesep, 'PreProcessedData']
  PROCESSED_DATA_PATH = [DATA_PATH, filesep, 'ProcessedData']
  RAW_DYNAMICS_DATA_PATH = [DATA_PATH, filesep, 'RawDynamicsData']
  DYNAMICS_RESULTS_PATH = [DATA_PATH, filesep, 'DynamicsResults']

  INSTALLATION_FILES_DIR = 'InstallationFiles'
  INSTALL_COMPUTER_FOLDERS_PATH = [INSTALLATION_FILES_DIR, filesep, 'InstallComputerFolders.xlsx']
  COMPUTER_FOLDERS_PATH = [rootPath, filesep, 'ComputerFolders.xlsx']
  INSTALL_MOVIE_DATABASE_PATH = [INSTALLATION_FILES_DIR, filesep, 'InstallMovieDatabase.xlsx']
  MOVIE_DATABASE_PATH = [DYNAMICS_RESULTS_PATH, filesep, 'MovieDatabase.xlsx']

  MS2CODE_PATH = mRNADynamicsPath

  mkdir(DATA_PATH)
  mkdir(PREPROCESSED_DATA_PATH)
  mkdir(PROCESSED_DATA_PATH)
  mkdir(RAW_DYNAMICS_DATA_PATH) %(old RawData folder)
  mkdir(DYNAMICS_RESULTS_PATH) %(old DropboxFolder)

  cflag = 0; %ComputerFolders.xlsx has not been made yet.
  %Copy the files to the different folders:
  if ~exist(COMPUTER_FOLDERS_PATH)
      copyfile(INSTALL_COMPUTER_FOLDERS_PATH, COMPUTER_FOLDERS_PATH)
  else
      warning([COMPUTER_FOLDERS_PATH, ' already exists. Not overwriting.'])
      cflag = 1;
  end

  %Edit ComputerFolders.XLSX
  [num,txt] = xlsread(COMPUTER_FOLDERS_PATH);

  txt{1, end+1} = getComputerName();
  txt{2, end} = getUserName();
  txt{3, end} = RAW_DYNAMICS_DATA_PATH;
  txt{4, end} = PREPROCESSED_DATA_PATH;
  txt{5, end} = PROCESSED_DATA_PATH;
  txt{6, end} = DYNAMICS_RESULTS_PATH;
  txt{7, end} = MS2CODE_PATH;

  %Save the XLS file
  if ispc && ~cflag
      xlswrite(COMPUTER_FOLDERS_PATH, txt);
  else
      if ~cflag
          disp('Warning: Macs and Linux cannot generate the XLS files.')
          disp('(1) Re-run using "txt=InstallmRNADynamics".')
          disp('(2) Type "open txt".')
          disp('(3) Copy and paste the data into a new file in Excel.')
          disp('(4) Get rid of all '' in the file.')
          disp('(5) Save as "ComputerFolders.xlsx" in folder "LivemRNAFISH.')
      end
  end
      
  %Copy MovieDatabase.XLSX
  if ~exist(MOVIE_DATABASE_PATH)
      copyfile(INSTALL_MOVIE_DATABASE_PATH, MOVIE_DATABASE_PATH)
  else
      warning([MOVIE_DATABASE_PATH, ' already exists. Not overwriting.'])
  end

  %Add the right folders to the path. This will be done as a startup file in
  %the user's folder

  DependenciesFolder = [mRNADynamicsPath, filesep, 'dependencies']
  DeprecatedFolder = [mRNADynamicsPath, filesep, 'deprecated']
  TestClassesFolder = [mRNADynamicsPath, filesep, 'testClasses']
  TrackingFolder = [mRNADynamicsPath, filesep, 'Tracking'];
  SubfunctionsFolder = [TrackingFolder, filesep, 'subfunctions'];
  LineageCodeFolder = [mRNADynamicsPath, filesep, 'LineageCode'];

  Output{1} = ['path(''', PREPROCESSED_DATA_PATH, ''',path);'];
  Output{2} = ['addpath(genpath(''', DependenciesFolder, '''))'];
  Output{3} = ['path(''', rootPath, ''',path);'];
  Output{4} = ['path(''', mRNADynamicsPath, ''',path);'];
  Output{5} = ['path(''', TrackingFolder, ''',path);'];
  Output{6} = ['path(''', LineageCodeFolder, ''',path);'];
  Output{7} = ['path(''', SubfunctionsFolder, ''',path);'];
  Output{8} = ['addpath(genpath(''', DeprecatedFolder, '''))'];
  Output{9} = ['path(''', DYNAMICS_RESULTS_PATH, ''',path);'];
  Output{10} = ['addpath(genpath(''', TestClassesFolder, '''))'];

  disp(['Will create startup.m with this contents: ', Output])

  %Create the startup.m file
  StartUpPath = userpath; 
  if isempty(userpath)
     display('Path for this specific user was not found. Please locate it in your documents folder.')
     display('In Windows, got to "My Documents\Matlab"')
     StartUpPath = uigetdir;
     userpath(StartUpPath);
  end

  fid = fopen([StartUpPath(1:end-1),filesep,'startup.m'], 'a');
  errmsg = '';
  if fid < 0
      fid = fopen([StartUpPath,filesep,'startup.m'], 'a');
  end
  while fid < 0 
     disp(errmsg);
     disp('Please find your user''s Matlab folder. Maybe "My Documents\Matlab" on Windows?')
     d = uigetdir;
     [fid,errmsg] = fopen(d);
  end

  for i=1:length(Output)
      fprintf(fid, '%s \n', Output{i});
  end
  fclose(fid);

  %Switch Matlab over to the Java release from the FIJI packaged with this
  %repository
  if ispc
      [status,~] = dos(['setx MATLAB_JAVA ',mRNADynamicsParentFolder,'\mRNADynamics\Fiji.app\java\win64\jdk1.8.0_66\jre']);
      if status
          warning('Something went wrong setting Java environment variable. Talk to Armando.')
      end
  else    
      warning(['Please note that Weka integration is not supported outside of Windows. ',...
              'Talk to Armando if you need this.'])
  end

  msgbox('Run "startup" from the command line or restart Matlab to finish the installation')
          
  warning('on','MATLAB:MKDIR:DirectoryExists')
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