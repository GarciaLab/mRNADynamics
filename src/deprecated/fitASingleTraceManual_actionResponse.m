% FitASingleTraceManual - this should be merged with the function
% FitInitialSlope above in the future.
function [frameIndex,Coefficients,ErrorEstimation,numberOfFramesUsedForFit] = ...
    fitASingleTraceManual(currentParticle,Particles,Spots,currentChannel,...
    schnitzcells,ElapsedTime,nuclearCycleBoundaries,correspondingNCInfo,...
    averagingLength, FramesToFit, FrameIndicesToFit, varargin)
% DESCRIPTION
% Fits a line to the initial slope of a fluorescent trace, using polyfit.
%
% ARGUMENTS
% currentParticle : particle index of interest
% Particles : a cell variable holding all the particles
% Spots : a cell variable holding all the spots
% currentChannel : index of current channel
% ElapsedTime : array of time elapsed since start of the movie
% nuclearCycleBoundaries : array of nuclearCycleBoundaries (in minutes) 
%                          corresponding to ElapsedTime 
% correspondingNCInfo : corresponding nc to the ElapsedTime array
% averagingLength : length of averaging done by the movmean function. Default is 3
% FramesToFit : Frames to fit the inital slope, manually defined from the
% upstream function, which is fitInitialSlope.
%
% OPTIONS

% Author (contact): Yang Joon Kim (yjkim90@berkeley.edu)
% This is edited from Emma Luu (emma_luu@berkeley.edu)'s code,
% fitASingleTrace.
% Created: 1/9/19
% Last Updated: 1/9/19 (YJK)
% Documented by: Yang Joon Kim (yjkim90@berkeley.edu)

%% Initializing the options
minimumLength = 3; % minimum length of trace
useDefaultTimeShift = 1;
fitApproved = 0;
%% Checking varargin
if ~isempty(varargin)
    for i=1:length(varargin)
        if strcmpi(varargin{i},'fitApproved')
            fitApproved = 1;
        end
    end
end
%% Getting particle information
% Amplitude ---------------------------------------------------------------
% getting frame information
[frame,~,ampIntegral3,~,~,~,~,~,~,~,~, ~,~, ampIntegralGauss3D,~]=...
    GetParticleTrace(currentParticle,...
    Particles{currentChannel},Spots{currentChannel});
currentLength = length(frame);

% performing moving average
if sum(isnan(ampIntegralGauss3D))
    disp('Note: Could not use the 3D guassian intensity fits so ampIntegral3 will be used')
    smoothedAmp = movmean(ampIntegral3,averagingLength);
else
    smoothedAmp = movmean(ampIntegralGauss3D,averagingLength);
end

% Time --------------------------------------------------------------------
% getting the corresponding time of the trace
currentTimeArray = ElapsedTime(frame); % Units to seconds
ncPresent = unique(correspondingNCInfo(frame));
% below subtracts 8 because the first element correspond to nc 9
timeOfFirstNC = nuclearCycleBoundaries(ncPresent(1)-8);
try
    nucleusFirstTime = ElapsedTime(...
        schnitzcells(Particles{currentChannel}(currentParticle).Nucleus).frames(1));
    timeShift = nucleusFirstTime;
catch
    timeShift = timeOfFirstNC;
end
% adjusting frameRange to have time 0 be the start of the
% first nuclear cycle the particle appears in or the first time point of
% the nucleus it is assigned to
if useDefaultTimeShift
    currentTimeArray = currentTimeArray - timeShift;
elseif useAnaphase
    currentTimeArray = currentTimeArray - timeOfFirstNC;
% leaving a structure such that you could also use this for any other
% time/event to shift...can rework it to be an if/else with ~useAnaphase as
% the condition. 11/20/18 - EL
end


%% Fitting process

% initializing outputs of the function
frameIndex = [];
Coefficients = [];
ErrorEstimation = struct('R',{},'df',{},'normr',{});
numberOfFramesUsedForFit = [];
 
% Fitting only when the trace has more than the minimumLength, which is 3
% as a default. This should be edited as an option in the future.
    if currentLength > minimumLength

        %numberOfFramesUsedForFit = zeros(1,numberOfGroups);
        % number of Frames to fit, this is given as an input
        numberOfFramesUsedForFit = length(FrameIndicesToFit);

        correspondingFrameIndices = FrameIndicesToFit;

        % saving the boundaries of correspondingFrameIndexes in frameIndex
        frameIndex(1,:)...
            = [correspondingFrameIndices(1) correspondingFrameIndices(end)];
        
        currentAmpSegment = smoothedAmp(correspondingFrameIndices);
        if fitApproved
            currentXSegment = currentTimeArray(FrameIndicesToFit);
        else
            currentXSegment = frame(FrameIndicesToFit); % Frames (not minutes)
        end
        [Coefficients(1,:),...
            ErrorEstimation(1)]...
            =  polyfit(currentXSegment,currentAmpSegment,1);
            
        currentYSegment = ...
            polyval(Coefficients(1,:),currentXSegment); 
        % this is for the predicted line
        denominator = sum((currentAmpSegment - mean(currentAmpSegment)).^2);
        normOfResiduals = ErrorEstimation(1).normr;
        RSquared = 1 - (normOfResiduals^2)/denominator;
        errorArray = ones(1,length(currentXSegment)).*...
            normOfResiduals./numberOfFramesUsedForFit(1); %EL normalized by number of points included
    end
end

