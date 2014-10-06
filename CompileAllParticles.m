function CompileAllParticles(DataType,varargin)

%Compiles all the data sets with the "READY" or "ApproveAll" labels on the XLS file.
%A DataType pointing to the relevant XLS spreadsheet needs to be included.

%Parameters:
%ForceAP  -  Force AP detection even if it's there already.
%SkipTraces - Don't output the individual traces
%SkipFluctuations - Don't generate the plots of the correlation of signal
%                   and offset


%% Folders

%Get some of the default folder
[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders;

if strcmp(DataType,'hbBAC')|strcmp(DataType,'Eve2')|strcmp(DataType,'snaBAC')|...
        strcmp(DataType,'snaBACNoPrimary')|strcmp(DataType,'snaBACNoShadow')|strcmp(DataType,'P2PPausing')|...
        strcmp(DataType,'hbNoPrimary')|strcmp(DataType,'kniBAC')|strcmp(DataType,'hbNoShadow')|...
        strcmp(DataType,'kniNoPrimary')|strcmp(DataType,'kniBAC')|strcmp(DataType,'kniNoShadow')
    [SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
        DetermineLocalFolders('2014-03-15-HbBACA');
    XLSName='DataStatusPausing.xlsx';
elseif strcmp(DataType,'zld')|strcmp(DataType,'Shelby')
    [SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
        DetermineLocalFolders('2013-09-09-zld-X1');
    XLSName='DataStatus.xlsx';
elseif strcmp(DataType,'MCP-GFP 5'' Data')
    error('Take care of this')
else
    error('Add this data type to the code')
end




%% Compile each data set


%Get the XLS info
[StatusNum,StatusTxt]=xlsread([DropboxFolder,filesep,XLSName],DataType);

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
 