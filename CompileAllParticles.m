function CompileAllParticles(DataType,varargin)

%Compiles all the data sets with the "READY" or "ApproveAll" labels on the XLS file.
%A DataType pointing to the relevant XLS spreadsheet needs to be included.

%Parameters:
%ForceAP  -  Force AP detection even if it's there already.
%SkipTraces - Don't output the individual traces
%SkipFluctuations - Don't generate the plots of the correlation of signal
%                   and offset


%DataType is the tab in the XLS file. The code figures out which XLS file
%and folders to use.

%Get some of the default folders
[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

%Now, get a list of all possible other Dropbox folders
[Dummy,XLS]=xlsread([MS2CodePath,filesep,'..',filesep,'ComputerFolders.XLSX']);
DropboxRows=find(~cellfun(@isempty,strfind(XLS(:,1),'Dropbox')));
DropboxFolders=XLS(DropboxRows,2);

%Look in DataStatus.XLSX in each DropboxFolder and find the tab given by
%the input variable DataType.
DataStatusToCheck=[];
for i=1:length(DropboxFolders)
    DDataStatus=dir([DropboxFolders{i},filesep,'DataStatus.*']);
    if length(DDataStatus)>1
        error(['More than one DataStatus.XLS found in folder ',DropboxFolders{i}])
    elseif length(DDataStatus)==1
        [Dummy,Sheets] = xlsfinfo([DropboxFolders{i},filesep,DDataStatus(1).name]);
        FindSheets=strcmpi(Sheets,lower(DataType));
        if sum(FindSheets)==1
            DataStatusToCheck=[DataStatusToCheck,i];
        end
    end
end
%Check that there aren't two DataStatus files with the same tab name
if length(DataStatusToCheck)>1
    error(['More than one DataStatus.XLSX found with a tab named ',DataType])
elseif length(DataStatusToCheck)==0
    error(['No DataStatus.XLSX found with a tab named ',DataType])
end

%Redefine the DropboxFolder according to the DataStatus.XLSX we'll use
DropboxFolder=DropboxFolders{DataStatusToCheck};
%Now, load the DataStatus.XLSX
D=dir([DropboxFolder,filesep,'DataStatus.*']);
[StatusNum,StatusTxt]=xlsread([DropboxFolder,filesep,D(1).name],DataType);


CompileRow=find(strcmp(StatusTxt(:,1),'AnalyzeLiveData Compile Particles'));
CompiledSets=find(strcmp(StatusTxt(CompileRow,:),'READY')|...
    strcmp(StatusTxt(CompileRow,:),'ApproveAll'));

for i=1:length(CompiledSets)
    display(['Compiling set ',num2str(i),' of ',num2str(length(CompiledSets))])
    SetName=StatusTxt{6,CompiledSets(i)};
    Quotes=strfind(SetName,'''');
    Prefix=SetName((Quotes(1)+1):(Quotes(end)-1));
    
    if isempty(varargin)
       CompileParticles(Prefix)
    else
        InputArgument=[];
        for j=1:length(varargin)
            InputArgument=[InputArgument,'''',varargin{j},''', '];
        end
        InputArgument=InputArgument(1:end-2);
            
        
        eval(['CompileParticles(Prefix,',InputArgument,');'])
    end
        
   
end
 