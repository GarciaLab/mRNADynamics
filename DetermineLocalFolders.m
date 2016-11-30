function [SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
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
    
%HG: I'm getting rid of the hardcoding of dropbox folders.
%     if strcmp(XLSTxt{PrefixRow,DropboxFolderColumn},'Hernan')
%         DropboxString='DropboxHernan';
%     elseif strcmp(XLSTxt{PrefixRow,DropboxFolderColumn},'Albert+Hernan')
%         DropboxString='DropboxAlbert';
%     elseif strcmp(XLSTxt{PrefixRow,DropboxFolderColumn},'Albert+Emilia')
%         DropboxString='DropboxEmilia';
%     elseif strcmp(XLSTxt{PrefixRow,DropboxFolderColumn},'Jacques+Hernan')
%         DropboxString='DropboxJacques';
%     elseif strcmp(XLSTxt{PrefixRow,DropboxFolderColumn},'HGLab')
%         DropboxString='DropboxHGLab';
%     elseif strcmp(XLSTxt{PrefixRow,DropboxFolderColumn},'Heinrich')
%         DropboxString='DropboxHeinrich';
%     elseif strcmp(XLSTxt{PrefixRow, DropboxFolderColumn}, 'Default')
%         DropboxString = 'DropboxFolder';
%     elseif strcmp(XLSTxt{PrefixRow, DropboxFolderColumn}, 'TwoColor')
%         DropboxString = 'DropboxTwoColor';
%         % ES 2013-10-06
%     else
%         error('Dropbox folder for this type of experiment not found. Check MovieDatabase')
%     end
    
    DropboxRow=find(strcmp(XLS(:,1),DropboxString));
end
    
if isempty(DropboxRow)
    error('Dropbox folder for this type of experiment not found. Check MovieDatabase')
end

DropboxFolder=XLS{DropboxRow,ComputerColumn};  

    
    
