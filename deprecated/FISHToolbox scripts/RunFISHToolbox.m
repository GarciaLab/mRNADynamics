function RunFISHToolbox(Prefix,Thresholds)
% RunFISHToolbox(Prefix,Thresholds)
%
% DESCRIPTION
% This function runs the FISH Toolbox steps on our live mRNA data.
%
% PARAMETERS
% Prefix: Prefix of the data set to analyze
% Threshold: Threshold to be used. Should be kept at ~90-200 for lattice
%           light-sheet data, and at ~5-10 for confocal data (Leica SP8).
%           If left empty, then the code just generates the DoG files.
%
% OPTIONS
% None
%
% CONTROLS
% None
%
% OUTPUT
% None
%
% Author (contact): Unknown
% Created: Unknown
% Last Updated: Unknown
%
% Commented by: Meghan Turner (meghan_turner@berkeley.edu), 6/8/17




%Get the relevant folders now:
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
    DetermineLocalFolders(Prefix);

%Figure out what type of experiment we have
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder)

if strcmp(ExperimentType,'1spot')
    NChannels=1;
elseif strcmp(ExperimentType,'2spot')
    NChannels=1;
elseif strcmp(ExperimentType,'2spot1color')
    NChannels=1;
elseif strcmp(ExperimentType,'2spot2color')
    NChannels=2;
elseif strcmp(ExperimentType,'inputoutput')
    NChannels=1;
else
    error('Experiment type not recognized in MovieDatabase')
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