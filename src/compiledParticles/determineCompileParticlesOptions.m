function [Prefix, ForceAP, SkipTraces, SkipFluctuations, SkipFits, SkipMovie, ...
    SkipAll, ApproveAll, MinParticles, minTime, ROI, noHist, ...
    ROI1, ROI2, slimVersion, manualSingleFits, optionalResults,...
    yToManualAlignmentPrompt, minBinSize, edgeWidth] = determineCompileParticlesOptions(varargin)
%DETERMINECOMPILEPARTICLESOPTIONS Summary of this function goes here
%   Detailed explanation goes here

varargin = varargin{1};

%Look at the input parameters and use defaults if missing
Prefix='';
ForceAP=false;      %Force AP detection even if it's already there

SkipTraces=false;   %Do not output the individual traces.
SkipFluctuations=false;  %Do not generate the plots of correlations of fluctuations and offset
SkipFits=false;        %Do not run and savegenerate the 
SkipMovie=false;        %Do not generate the movie
SkipAll=false;         %Do not do other things 

ApproveAll= true;       %Only use manually approved particles
MinParticles=4;
minTime = 1;
ROI=false; % No ROI
noHist = false; 
slimVersion = false;
ROI1 = -1; % no ROI
ROI2 = -1; % no ROI
manualSingleFits = false; % no manually fitted single fits
optionalResults = ''; %different dropbox folder
yToManualAlignmentPrompt = false; %this is an option for addparticleposition
minBinSize = .7; % fraction of the median bin size allowed
edgeWidth = 2.12; %um?


% Checking Varargin 
if isempty(varargin)%looks for the folder to analyze
    FolderTemp=uigetdir(DefaultDropboxFolder,'Select folder with data to analyze');
    Dashes=strfind(FolderTemp,filesep);
    Prefix=FolderTemp((Dashes(end)+1):end);
else
    Prefix=varargin{1};
    for i=2:length(varargin)
        if strcmpi(varargin{i},'ForceAP')
            ForceAP=1;
        elseif strcmpi(varargin{i},'SkipTraces')
            SkipTraces=varargin{i+1};
        elseif strcmpi(varargin{i},'SkipFluctuations')
            SkipFluctuations=varargin{i+1};
        elseif strcmpi(varargin{i},'SkipFits')
            SkipFits=varargin{i+1};
        elseif strcmpi(varargin{i},'SkipMovie')
            SkipMovie=varargin{i+1};
        elseif strcmpi(varargin{i},'SkipAll')
            SkipTraces=1;
            SkipFluctuations=1;
            SkipFits=1;
            SkipMovie=1;
            SkipAll=1;
        elseif strcmpi(varargin{i},'ApproveAll')
            ApproveAll=1;
        elseif strcmpi(varargin{i},'KeepAll') % Opposite of the 'ApproveAll' flag, which is currently set to 1 by default
            ApproveAll=0;
        elseif strcmpi(varargin{i},'noHist')
            noHist = 1;
        elseif strcmp(varargin{i},'MinParticles')
            if ~isnumeric(varargin{i+1})
                error('Wrong input parameters. After ''MinParticles'' you should input the desired minimum number of particles per approved AP bin')
            else
                MinParticles=varargin{i+1};
            end
        elseif strcmpi(varargin{i},'minBinSize')
            if ~isnumeric(varargin{i+1})
                error('Wrong input parameters. After ''minBinSize'' you should input the desired minimum number of frames per particle.')
            else
                minBinSize=varargin{i+1};
            end
        elseif strcmpi(varargin{i},'ROI')
            ROI = 1;
            if ~isnumeric(varargin{i+1})||~isnumeric(varargin{i+2})
                error('Wrong input parameters. After ''ROI'' you should input the y-threshold of ROI ')
            else
                ROI1=varargin{i+1};
                ROI2=varargin{i+2};
            end
        elseif strcmpi(varargin{i}, 'slimVersion')
            slimVersion = 1;
        elseif strcmpi(varargin{i}, 'edgeWidth')
            edgeWidth = varargin{i+1};
        elseif strcmpi(varargin{i}, 'manualSingleFits')
            manualSingleFits = 1;
        elseif strcmpi(varargin{i}, 'optionalResults')
            optionalResults = varargin{i+1};
        elseif strcmpi(varargin{i}, 'yToManualAlignmentPrompt')
            yToManualAlignmentPrompt = 1;
        end
    end
    
end
end

