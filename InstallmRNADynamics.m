function InstallmRNADynamics

%This function creates all the required folders to run the FISHToolbox.
%IMPORTANT: This needs to be run from the 'mRNADynamics' folder.


%Check that we are in the right folder
D=dir('InstallmRNADynamics.m');
if isempty(D)
    error('Run this code from ''mRNADynamics''')
end


if ~exist(['..',filesep,'Data',filesep,'PreProcessedData'])
    error('Install the FISHToolbox first')
end


%Create the different folders we need
%RawDynamicsData is the old RawData folder
mkdir(['..',filesep,'Data',filesep,'RawDynamicsData'])
%DynamicsResults is the old DropboxFolder
mkdir(['..',filesep,'Data',filesep,'DynamicsResults'])

%Copy the files to the different folders:
copyfile(['InstallationFiles',filesep,'InstallComputerFolders.xlsx'],...
    ['..',filesep,'ComputerFolders.xlsx'])

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
    [Dummy, user_name] = system('whoami'); % exists on every unix that I know of
    % on my mac, isunix == 1
elseif ispc
    [Dummy, user_name] = system('echo %USERDOMAIN%\%USERNAME%'); % Not as familiar with windows,
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
xlswrite(['..',filesep,'ComputerFolders.xlsx'],txt);


%Copy MovieDatabase.XLSX
if ~exist(['..',filesep,'Data',filesep,'DynamicsResults',filesep,'MovieDatabase.xlsx'])
    copyfile(['InstallationFiles',filesep,'InstallMovieDatabase.xlsx'],...
        ['..',filesep,'Data',filesep,'DynamicsResults',filesep,'MovieDatabase.xlsx'])
else
    warning('MovieDatabase.xlsx already exist, we are not overwriting.')
end


%Add the right folders to the path
%mRNADynamics
path(cd,path);     
CurrentFolder=cd;
%mRNADynamics\Tracking
cd('Tracking')
path(cd,path);
%mRNADynamics\Tracking\subfunctions
cd('subfunctions')
path(cd,path);
cd(CurrentFolder);
%Folder up from mRNADynamics
cd('..')
path(cd,path);     
cd(CurrentFolder);





