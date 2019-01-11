function [Prefix, ForceAP, SkipTraces, SkipFluctuations, SkipFits, SkipMovie, ...
    SkipAll, ApproveAll, MinParticles, minTime, ROI, intArea, noHist, ...
    doSingleFits, ROI1, ROI2, slimVersion] = determineCompileParticlesOptions(varargin)
%DETERMINECOMPILEPARTICLESOPTIONS Summary of this function goes here
%   Detailed explanation goes here

varargin = varargin{1};

%Look at the input parameters and use defaults if missing
Prefix='';
ForceAP=0;      %Force AP detection even if it's already there
SkipTraces=0;   %Do not output the individual traces.
SkipFluctuations=0;  %Do not generate the plots of correlations of fluctuations and offset
SkipFits=0;         %Do not generate the fit output (but still does the fit)
SkipMovie=0;        %Do not generate the movie
SkipAll=0;          %Do not do other things 
ApproveAll=0;       %Only use manually approved particles
MinParticles=4;
minTime = 1;
ROI=0; % No ROI
intArea = 109; %pixels. default for 220nm x 220nm zoom. for 70nm use 437 pixels. 
noHist = 0; 
doSingleFits = 0;
slimVersion = 0;
ROI1 = -1; % no ROI
ROI2 = -1; % no ROI


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
            SkipTraces=1;
        elseif strcmpi(varargin{i},'SkipFluctuations')
            SkipFluctuations=1;
        elseif strcmpi(varargin{i},'SkipFits')
            SkipFits=1;
        elseif strcmpi(varargin{i},'SkipMovie')
            SkipMovie=1;
        elseif strcmpi(varargin{i},'SkipAll')
            SkipTraces=1;
            SkipFluctuations=1;
            SkipFits=1;
            SkipMovie=1;
            SkipAll=1;
        elseif strcmpi(varargin{i},'doSingleFits')
            doSingleFits=1;
        elseif strcmpi(varargin{i},'ApproveAll')
            ApproveAll=1;
        elseif strcmpi(varargin{i},'noHist')
            noHist = 1;
        elseif strcmp(varargin{i},'MinParticles')
            if ~isnumeric(varargin{i+1})
                error('Wrong input parameters. After ''MinParticles'' you should input the desired minimum number of particles per approved AP bin')
            else
                MinParticles=varargin{i+1};
            end
        elseif strcmpi(varargin{i},'intArea')
            if ~isnumeric(varargin{i+1})
                error('Wrong input parameters. After ''intArea'' you should input the desired number of pixels for intensity integration')
            else
                intArea=varargin{i+1};
            end
        elseif strcmpi(varargin{i},'MinTime')
            if ~isnumeric(varargin{i+1})
                error('Wrong input parameters. After ''MinTime'' you should input the desired minimum number of frames per particle.')
            else
                minTime=varargin{i+1};
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
        end
    end
    
end
end

