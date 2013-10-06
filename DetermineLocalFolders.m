function [SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders(varargin)



%This functions gives out the folder corresponding to each computer. If a
%Prefix is also included it will give out the corresponding DropboxFolder
%for the particular experiment. Otherwise it will give the default dropbox
%folder. Regardless, it assumes that MovieDatabase.xlsx is in the
%default Dropbox folder.

if isempty(varargin)
    %warning('No Prefix defined. Will output default Dropbox folder')
end


[Dummy,XLS]=xlsread('ComputerFolders.xlsx');

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

%Now load the corresponding folders
SourceRow=find(strcmp(XLS(:,1),'SourcePath'));
FISHRow=find(strcmp(XLS(:,1),'FISHPath'));
SchnitzRow=find(strcmp(XLS(:,1),'SchnitzcellsFolder'));
MS2CodeRow=find(strcmp(XLS(:,1),'MS2CodePath'));

%Assign the folders
SourcePath=XLS{SourceRow,ComputerColumn};
FISHPath=XLS{FISHRow,ComputerColumn};
MS2CodePath=XLS{MS2CodeRow,ComputerColumn};
SchnitzcellsFolder=XLS{SchnitzRow,ComputerColumn};


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
    
    PrefixRow=find(strcmp(XLSTxt(:,DataFolderColumn),[Prefix(1:10),'\',Prefix(12:end)])|...
        strcmp(XLSTxt(:,DataFolderColumn),[Prefix(1:10),'/',Prefix(12:end)]));
    
    if isempty(PrefixRow)
        error('Data set information not found in MovieDatabase.xlsx')
    end
    
    
    if strcmp(XLSTxt{PrefixRow,DropboxFolderColumn},'Albert+Hernan')
        DropboxString='DropboxAlbert';
    elseif strcmp(XLSTxt{PrefixRow,DropboxFolderColumn},'Jacques+Hernan')
        DropboxString='DropboxJacques';
    else
        error('Dropbox folder for this type of experiment not found. Check MovieDatabase')
    end
    
    DropboxRow=find(strcmp(XLS(:,1),DropboxString));
end
    

DropboxFolder=XLS{DropboxRow,ComputerColumn};    
    
    
