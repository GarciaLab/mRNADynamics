function RunFISHToolbox(Prefix,Thresholds)

%This function runs the FISH Toolbox steps on our live mRNA data. If not
%threshold is provided it runs it in DoG mode for all channels.



%Get the relevant folders now:
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
    DetermineLocalFolders(Prefix);

%Figure out what type of experiment we have
[XLSNum,XLSTxt]=xlsread([DropboxFolder,filesep,'MovieDatabase.xlsx']);
DataFolderColumn=find(strcmp(XLSTxt(1,:),'DataFolder'));
ExperimentTypeColumn=find(strcmp(XLSTxt(1,:),'ExperimentType'));
Channel1Column=find(strcmp(XLSTxt(1,:),'Channel1'));
Channel2Column=find(strcmp(XLSTxt(1,:),'Channel2'));

% Convert the prefix into the string used in the XLS file
Dashes = strfind(Prefix, '-');
PrefixRow = find(strcmp(XLSTxt(:, DataFolderColumn),...
    [Prefix(1:Dashes(3)-1), '\', Prefix(Dashes(3)+1:end)]));
if isempty(PrefixRow)
    PrefixRow = find(strcmp(XLSTxt(:, DataFolderColumn),...
        [Prefix(1:Dashes(3)-1), '/', Prefix(Dashes(3)+1:end)]));
    if isempty(PrefixRow)
        error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
    end
end

ExperimentType=XLSTxt(PrefixRow,ExperimentTypeColumn);

if strcmp(ExperimentType,'1spot')
    NChannels=1;
elseif strcmp(ExperimentType,'2spot')
    NChannels=1;
elseif strcmp(ExperimentType,'2spot2color')
    NChannels=2;
elseif strcmp(ExperimentType,'inputoutput')
    NChannels=1;
else
    error('Experiment type not recognized in MovieDatabase.XLSX')
end

if ~exist('Thresholds')
    Thresholds=ones(1,NChannels)*inf;
else
    if length(Thresholds)~=NChannels
        error('Number of channels in movie does not match number of thresholds input')
    end
end



%Start the matlab workers for the FISH analysis code

%Try matlabpool only for MATLAB versions prior to 2015
year15 = datenum(2015,01,01);
[v,d] = version;
d=datenum(d);
if d < year15
    try
        matlabpool
    catch
        display('matlabpool already running')
    end
    
else
    try
        parpool
    catch
        display('matlabpool already running')
    end  
end

cd([FISHPath])%when it runs this line it cannot longer find 'params_mRNADynamics' in the new folder
analyzeDataLibrary('fad',@(x)tagged(x,'id',[Prefix,'_']),'params_mRNADynamics',Thresholds)
cd([MS2CodePath])