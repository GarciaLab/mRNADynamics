function txt=InstallmRNADynamics
% function txt=InstallmRNADynamics
%
% DESCRIPTION
% This function creates all the required folders to run the MS2 movie
% analysis pipeline. IMPORTANT: This needs to be run from the 'mRNADynamics' folder.
%
% ARGUMENTS
% None.
%
% OUTPUT
% txt: A cell containing information within ComputerFolders.xlsx 
%
% %Author (contact): Hernan Garcia (hgarcia@berkeley.edu)
% Created: Unknown
% Last Updated: 6/19/17
%

warning('off','MATLAB:MKDIR:DirectoryExists')

%Check that we are in the right folder
D=dir('InstallmRNADynamics.m');
if isempty(D)
    error('Run this code from ''mRNADynamics''')
end

mkdir(['..',filesep,'Data'])
mkdir(['..',filesep,'Data',filesep,'PreProcessedData'])
mkdir(['..',filesep, 'Data', filesep, 'ProcessedData'])

%Create the different folders we need
%RawDynamicsData is the old RawData folder
mkdir(['..',filesep,'Data',filesep,'RawDynamicsData'])
%DynamicsResults is the old DropboxFolder
mkdir(['..',filesep,'Data',filesep,'DynamicsResults'])

cflag = 0; %ComputerFolders.xlsx has not been made yet.
%Copy the files to the different folders:
if ~exist(['..',filesep,'ComputerFolders.xlsx'])
    copyfile(['InstallationFiles',filesep,'InstallComputerFolders.xlsx'],...
    ['..',filesep,'ComputerFolders.xlsx'])
else
    warning('ComputerFolders.xlsx already exists. Not overwriting.')
    cflag = 1;
end

%Edit ComputerFolders.XLSX
[num,txt]=xlsread(['..',filesep,'ComputerFolders.xlsx']);
%Find out the computer name
[ret, name] = system('hostname');  
if ret ~= 0,  
   if ispc  
      name = getenv('COMPUTERNAME');  
   else  
      name = getenv('HOSTNAME');  
   end  
end  
name = lower(name);
txt{1,end+1}=name;
%Find out the username
if isunix
    [~, user_name] = system('whoami'); % exists on every unix that I know of
    % on my mac, isunix == 1
elseif ispc
    [~, user_name] = system('echo %USERDOMAIN%\%USERNAME%'); % Not as familiar with windows,
                            % found it on the net elsewhere, you might want to verify
end
txt{2,end}=user_name;
%Add the RawDynamicsData folder
CurrentFolder=cd;
cd(['..',filesep,'Data',filesep,'RawDynamicsData'])
txt{3,end}=cd;
cd(CurrentFolder)
%Add PreProcessedData
cd(['..',filesep,'Data',filesep,'PreProcessedData'])
txt{4,end}=cd;
cd(CurrentFolder)
%Add ProcessedData
cd(['..',filesep,'Data',filesep,'ProcessedData'])
txt{5,end}=cd;
cd(CurrentFolder)
%Add DynamicsResults
cd(['..',filesep,'Data',filesep,'DynamicsResults'])
txt{6,end}=cd;
cd(CurrentFolder)
%Add MS2Code
txt{7,end}=cd;
%Save the XLS file
if ispc && ~cflag
    xlswrite(['..',filesep,'ComputerFolders.xlsx'],txt);
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
if ~exist(['..',filesep,'Data',filesep,'DynamicsResults',filesep,'MovieDatabase.xlsx'])
    copyfile(['InstallationFiles',filesep,'InstallMovieDatabase.xlsx'],...
        ['..',filesep,'Data',filesep,'DynamicsResults',filesep,'MovieDatabase.xlsx'])
else
    warning('MovieDatabase.xlsx already exists, we are not overwriting.')
end


%Add the right folders to the path. This will be done as a startup file in
%the user's folder

%mRNADynamics
CurrentFolder=cd;
%PreProcessedFolder
cd(['..', filesep, 'Data',filesep,'PreProcessedData']);
PreProcessedFolder=cd;
cd(CurrentFolder)
%mRNADynamics\Tracking
cd('Tracking')
TrackingFolder=cd;
%mRNADynamics\Tracking\subfunctions
cd('subfunctions')
SubfunctionsFolder=cd;
cd(CurrentFolder);
%mRNADynamics\LineageCode
cd(CurrentFolder)
cd('LineageCode')
LineageCodeFolder=cd;
cd(CurrentFolder)
%Folder up from mRNADynamics
cd('..')
mRNADynamicsParentFolder=cd;
%Data\DynamicsResults
cd(['Data',filesep,'DynamicsResults'])
DynamicsResultsFolder=cd;
cd(CurrentFolder);

Output{1}=['path(''',PreProcessedFolder,''',path);'];
Output{2}=['addpath(genpath(''',CurrentFolder,filesep,'dependencies''))'];
Output{3}=['path(''',mRNADynamicsParentFolder,''',path);'];
Output{4}=['path(''',CurrentFolder,''',path);'];
Output{5}=['path(''',TrackingFolder,''',path);'];
Output{6}=['path(''',LineageCodeFolder,''',path);'];
Output{7}=['path(''',SubfunctionsFolder,''',path);'];
Output{8}=['addpath(genpath(''',CurrentFolder,filesep,'deprecated''))'];
Output{9}=['path(''',DynamicsResultsFolder,''',path);'];
Output{10}=['addpath(genpath(''',CurrentFolder,filesep,'testClasses''))'];




%Create the startup.m file
StartUpPath=userpath; 
if isempty(userpath)
   display('Path for this specific user was not found. Please locate it in your documents folder.')
   display('In Windows, got to "My Documents\Matlab"')
   StartUpPath=uigetdir;
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


% else
%     Answer=input('WARNING: startup.m already exist. Append, overwrite or cancel? (A/O/C): ','s')
%     if strcmp(lower(Answer),'a')
%         fid = fopen([StartUpPath(1:end-1),filesep,'startup.m'], 'a');
%         for i=1:length(Output)
%             fprintf(fid, '%s \n', Output{i});
%         end
%         fclose(fid);
%     elseif strcmp(lower(Answer),'o')
%         fid = fopen([StartUpPath(1:end-1),filesep,'startup.m'], 'wt');
%         for i=1:length(Output)
%             fprintf(fid, '%s \n', Output{i});
%         end
%         fclose(fid);
%     end
% end
% 

%I had to do this because it seems to take some time for the file to be
%found by Matlab after creating it. AR: This seems to happen more often
%than not. I haven't found a function to bail out 
%the script if it hangs, so I'm defaulting to having the user run startup manually. 
% try
%     startup;
% catch
    msgbox('Run "startup" from the command line or restart Matlab to finish the installation')
% end
        
warning('on','MATLAB:MKDIR:DirectoryExists')