function txt=InstallmRNADynamics

%This function creates all the required folders to run the FISHToolbox.
%IMPORTANT: This needs to be run from the 'mRNADynamics' folder.

warning('off','MATLAB:MKDIR:DirectoryExists')

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
if ispc
    xlswrite(['..',filesep,'ComputerFolders.xlsx'],txt);
else
    display('Warning: Macs and Linux cannot generate the XLS files.')
    display('(1) Re-run using "txt=InstallmRNADynamics".')
    display('(2) Type "open txt".')
    display('(3) Copy and paste the data into a new file in Excel.')
    display('(4) Get rid of all '' in the file.')
    display('(5) Save as "ComputerFolders.xlsx" in folder "LivemRNAFISH.')
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



Output{1}=['path(''',CurrentFolder,''',path);'];
Output{2}=['path(''',TrackingFolder,''',path);'];
Output{3}=['path(''',SubfunctionsFolder,''',path);'];
Output{4}=['path(''',mRNADynamicsParentFolder,''',path);'];
Output{5}=['path(''',DynamicsResultsFolder,''',path);'];
Output{6}=['path(''',LineageCodeFolder,''',path);'];


%Create the startup.m file
StartUpPath=userpath;
%if ~exist([StartUpPath(1:end-1),filesep,'startup.m'])
fid = fopen([StartUpPath(1:end-1),filesep,'startup.m'], 'a');

for i=1:length(Output)
    fprintf(fid, '%s \n', Output{i});
end
fclose(fid);
    
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
%found by Matlab after creating it
try
    startup;
catch
    display('Run "startup" to finish the installation')
end
        

warning('on','MATLAB:MKDIR:DirectoryExists')



        


