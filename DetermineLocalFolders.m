function [SourcePath, FISHPath, DropboxFolder, MS2CodePath, PreProcPath] =...
    DetermineLocalFolders(varargin)

CONFIG_CSV_PATH = ['..', filesep, 'ComputerFolders.csv'];

if exist(CONFIG_CSV_PATH, 'file')
  disp('Using new CSV configuration file')
  configValues = openCSVFile(CONFIG_CSV_PATH);

  SourcePath = getConfigValue(configValues, 'SourcePath')
  FISHPath = getConfigValue(configValues, 'FISHPath')
  DropboxFolder = getConfigValue(configValues, 'DropboxFolder')
  MS2CodePath = getConfigValue(configValues, 'MS2CodePath')
  PreProcPath = getConfigValue(configValues, 'PreProcPath')

  if isempty(varargin)
    warning('No Prefix specified. Using default Dropbox folder')
    return
  end
  
  %% We need to look for the dropbox folder specified by the provided prefix
  Prefix = varargin{1};
  if isempty(regexp(Prefix, '^.{10}[\\\\/].*$'))
    error(['Prefix "', Prefix, '" does not match "yyyy-mm-dd/name". Please change it accordingly.'])
    % any 10 characters will work, not only yyyy-mm-dd,
    % but we enforce this in the error msg for simplicity
  end

  MOVIE_DATABASE_CSV_PATH = [DropboxFolder,filesep,'MovieDatabase.csv'];
  if exist(MOVIE_DATABASE_CSV_PATH, 'file')
    disp(['Using new CSV MovieDatabase file: ', MOVIE_DATABASE_CSV_PATH])
    dropboxFolderName = getDropboxFolderFromMovieDatabase_CSV(...
      MOVIE_DATABASE_CSV_PATH, Prefix);
  else
    % fallback to legacy XLS format
    dropboxFolderName = getDropboxFolderFromMovieDatabase(...
      [DropboxFolder,filesep,'MovieDatabase.xlsx'], Prefix);
  end

  DropboxFolder = getConfigValue(configValues, dropboxFolderName)

  return
end

%%
% Legacy xls support starts here - for backwards compatibility only
%%

%This functions gives out the folder corresponding to each computer. If a
%Prefix is also included it will give out the corresponding DropboxFolder
%for the particular experiment. Otherwise it will give the default dropbox
%folder. Regardless, it assumes that MovieDatabase.xlsx is in the
%default Dropbox folder.

if isempty(varargin)
    %warning('No Prefix defined. Will output default Dropbox folder')
end

[Dummy,XLS]=xlsread('ComputerFolders.xlsx');

%Check if there's only one computer defined. If so, we'll just use the
%second column. If not, ask for the computer name


if size(XLS,2)>2


    %Find out which computer this is. That will determine the folder structure.
    [ret, name] = system('hostname');  
    if ret ~= 0,  
       if ispc  
          name = getenv('COMPUTERNAME');  
       else  
          name = getenv('HOSTNAME');  
       end  
    end  
    name = lower(name); 


    %Find which computer we are dealing with:
    ComputerColumn=find(strcmp(XLS(1,:),name(1:end-1)));

    % Error when computer not found
    if isempty(ComputerColumn)
        disp('%%%%%%%%%%%%%%%%%%%%%')
        disp('Computer could not be found. Check host name or update ComputerFolders.xlsx')
        disp('%%%%%%%%%%%%%%%%%%%%%')
    end


    % ES 2013-10-27: queries user name if more than one user is defined for
    % this computer
    if length(ComputerColumn) > 1
        [Dummy, username] = system('echo %username%');
        UserRow = strcmp(XLS(:, 1), 'User Name');
        ComputerColumn = find(strcmp(XLS(UserRow, :), username(1:end-1)));
    end
elseif size(XLS,2)==2
    ComputerColumn=2;
else
    error('Cannot find folders in LivemRNA\ComputerFolders.xlsx. Were all steps followed in InstallmRNADynamics.m?')
end

%Now load the corresponding folders
SourceRow=find(strcmp(XLS(:,1),'SourcePath'));
FISHRow=find(strcmp(XLS(:,1),'FISHPath'));
MS2CodeRow=find(strcmp(XLS(:,1),'MS2CodePath'));
PreProcRow=find(strcmp(XLS(:,1),'PreProcPath'));


%Assign the folders
SourcePath=XLS{SourceRow,ComputerColumn};
FISHPath=XLS{FISHRow,ComputerColumn};
MS2CodePath=XLS{MS2CodeRow,ComputerColumn};
PreProcPath=XLS{PreProcRow,ComputerColumn};


%Deal with Dropbox

%Determine the default Dropbox folder
DefaultDropboxRow=find(strcmp(XLS(:,1),'DropboxFolder'));
DefaultDropboxFolder=XLS{DefaultDropboxRow,ComputerColumn};


if isempty(varargin)
    DropboxRow=find(strcmp(XLS(:,1),'DropboxFolder'));
else
    Prefix=varargin{1};
    
    %Figure the DropboxFolder corresponding to this
    [XLSNum,XLSTxt]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);

    DataFolderColumn=find(strcmp(XLSTxt(1,:),'DataFolder'));
    DropboxFolderColumn=find(strcmp(XLSTxt(1,:),'DropboxFolder'));
    
    % Convert the prefix into the string used in the XLS file
    Dashes = strfind(Prefix, '-');
    
    %We are making the prefix now indifferent about "\" vs "/"
    PrefixRow = find(strcmp(XLSTxt(:, DataFolderColumn),...
        [Prefix(1:Dashes(3)-1), '\', Prefix(Dashes(3)+1:end)]));
    if isempty(PrefixRow)    
        PrefixRow = find(strcmp(XLSTxt(:, DataFolderColumn),...
            [Prefix(1:Dashes(3)-1), '/', Prefix(Dashes(3)+1:end)]));
        if isempty(PrefixRow)
            error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
        end
    end
    
    if isempty(PrefixRow)
        error('Data set information not found in MovieDatabase.xlsx')
    elseif length(PrefixRow)>1
        error('Two data sets seem to have the same folder information. Check MovieDatabase.xlsx')
    end
    
    
    DropboxString=XLSTxt{PrefixRow,DropboxFolderColumn};
    
    %We can use the string "Default" or "DropboxFolder" for the default
    %Dropbox folder
    if strcmpi(DropboxString,'default')
        DropboxString='DropboxFolder';
    end
    

    DropboxRow=find(strcmp(XLS(:,1),DropboxString));
end
    
if isempty(DropboxRow)
    error('Dropbox folder for this type of experiment not found. Check MovieDatabase')
end

DropboxFolder=XLS{DropboxRow,ComputerColumn};  

end

%%
% Private functions
%% 
function value = getConfigValue(configuration, propertyName)
  indexArray = strfind(configuration, propertyName);
  propertyLabelIndex = find(not(cellfun('isempty', indexArray)));

  if isempty(propertyLabelIndex)
    error(['Property ''', propertyName, ''' not found in configuration. ',...
      'Check ComputerFolders and/or MovieDatabase.'])
  end

  value = configuration{propertyLabelIndex + 1};
end

function csv = openCSVFile(filePath)
  csvCell = textscan(fopen(filePath),'%s', 'delimiter', ',');
  csv = csvCell{1};
end

function dropboxFolderName = getDropboxFolderFromMovieDatabase(movieDatabasePath, prefix)
   %% TO-DO migrate MovieDatabase from Excel to CSV as well.
  [XLSNum,XLSTxt] = xlsread(movieDatabasePath);

  DataFolderColumn = find(strcmp(XLSTxt(1,:),'DataFolder'));
  DropboxFolderColumn = find(strcmp(XLSTxt(1,:),'DropboxFolder'));

  % we build a regex from the specified prefix,
  % in the form "10characters-slashOrBackslash-restOfThePrefix"
  % this enables support for prefixes using either '\' or '/' as separator
  % also, "10characters" at the beginning is usually "yyyy-mm-dd",
  % but it could be anything else
  prefixRegex = ['^', prefix(1:10), '[\\\\/]', prefix(12:end), '$'];
  regexMatches = regexp(XLSTxt(:, DataFolderColumn), prefixRegex);

  PrefixRow = find(~cellfun(@isempty, regexMatches));

  if isempty(PrefixRow) 
    error(['Data set "', prefix, '" not found in MovieDatabase.xlsx'])
  elseif length(PrefixRow) > 1
    error(['Rows ', sprintf('%d, ', PrefixRow), 'seem to have the same folder information. CheckMovieDatabase.xlsx'])
  end

  dropboxFolderName = XLSTxt{PrefixRow, DropboxFolderColumn};

  %We can use the string "Default" or "DropboxFolder" for the default Dropbox folder
  if strcmpi(dropboxFolderName, 'default')
    dropboxFolderName = 'DropboxFolder';
  end
end

function dropboxFolderName = getDropboxFolderFromMovieDatabase_CSV(movieDatabasePath, prefix)
end