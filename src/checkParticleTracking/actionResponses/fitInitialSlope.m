function [lineFit, Coefficients, fit1E, Particles, FramesToFit, FrameIndicesToFit] =...
    fitInitialSlope(CurrentParticle, Particles, Spots, CurrentChannel, schnitzcells, ...
    ElapsedTime, anaphaseInMins, correspondingNCInfo, traceFigAxes, Frames, anaphase, ...
    averagingLength, FramesToFit, FrameIndicesToFit, lineFit)
%fitInitialSlope 
% : This function grabs two points that you click, then fits using polyfit,
% and will give you the slope, error, and fitted line for overlaid plot.
%   Detailed explanation goes here

% Input parameters that should be defined :
% AverageLength, Time window for fitting(adjustable), 
% Use GUI for defining the inputs, as well as repeating the fitting until
% it's approved.

% Plug in inputs defined in GUI, 
% averagingLength, FramesToFit

%% Part1. fitting a line using the polyfit
    % Define some time parameters
    ncPresent = unique(correspondingNCInfo(Frames));
    % below subtracts 8 because the first element corresponds to nc 9
%     priorAnaphaseInMins = anaphaseInMins(ncPresent(1)-8);
    priorAnaphase = anaphase(ncPresent(1)-8); %frame
    if ~isempty(schnitzcells)&&~isempty(Particles{CurrentChannel}(CurrentParticle).Nucleus)
        nucleusFirstFrame = ElapsedTime(...
            schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).frames(1)); %min
    else
        nucleusFirstFrame = ElapsedTime(priorAnaphase); %min
        warning('No nucleus assigned to this particle. Using anaphase from moviedatabase as the first timepoint.')
    end
    currentTimeArray = ElapsedTime(Frames) - nucleusFirstFrame;
    
    % define the Frames for fitting
    [X,Y] = ginput(2); % pick two points (left, and right)
    if ~lineFit
        pos1 = Frames(find((Frames-X(1)).^2 == min((Frames-X(1)).^2)));
        pos2 = Frames(find((Frames-X(2)).^2 == min((Frames-X(2)).^2)));
        posIndex1 = find((Frames-X(1)).^2 == min((Frames-X(1)).^2));
        posIndex2 = find((Frames-X(2)).^2 == min((Frames-X(2)).^2));
        FramesToFit = pos1:pos2; % actual frames numbers used for linear fitting
        FrameIndicesToFit = posIndex1:posIndex2; % indices of those frames in the trace
    elseif lineFit
        pos1 = currentTimeArray(find((currentTimeArray-X(1)).^2 == min((currentTimeArray-X(1)).^2)));
        pos2 = currentTimeArray(find((currentTimeArray-X(2)).^2 == min((currentTimeArray-X(2)).^2)));
        posIndex1 = find((currentTimeArray-X(1)).^2 == min((currentTimeArray-X(1)).^2));
        posIndex2 = find((currentTimeArray-X(2)).^2 == min((currentTimeArray-X(2)).^2));
        FramesToFit = pos1:pos2; % actual frames numbers used for linear fitting
        FrameIndicesToFit = posIndex1:posIndex2; % indices of those frames in the trace
    end

 %try
    % currently shifted by the first frame of the assigned nucleus
    [frameIndex,Coefficients,ErrorEstimation,nFramesForFit] = ...
        fitASingleTraceManual(CurrentParticle,Particles,Spots,CurrentChannel,...
        schnitzcells,ElapsedTime,anaphaseInMins,correspondingNCInfo,...
        averagingLength, FramesToFit,FrameIndicesToFit);

    % plotting the fitted line
    currentXSegment = ElapsedTime(Frames(frameIndex(1):frameIndex(end)))-nucleusFirstFrame; % min after the previous anaphse
    currentYSegment = polyval(Coefficients,currentXSegment);
    % error of predicted line
    %          currentAmpSegment = AmpIntegral3(frameIndex(1):frameIndex(end));
    %                       denominator = sum((currentAmpSegment - mean(currentAmpSegment)).^2);
    %              RSquared = 1 - (normOfResiduals^2)/denominator;
    %              normOfResiduals = ErrorEstimation.normr;
    %              errorArray = ones(1,length(currentXSegment)).*...
    %                  normOfResiduals./nParticlesForFit; %EL normalized by number of points included
    hold(traceFigAxes,'on')
    %              fit1E = errorbar(traceFigAxes,ElapsedTime(Frames(frameIndex(1):frameIndex(end))),...
    %                  currentYSegment,errorArray,'.-','Color','red');
    %              to = -Coefficients(2) / Coefficients(1) + priorAnaphaseInMins;
    %              to = -Coefficients(2) / Coefficients(1) + priorAnaphaseInMins;
    to = -Coefficients(2) / Coefficients(1) + nucleusFirstFrame; % frame, not minutes for now
    %              fit1ETimeAxis = [to, ElapsedTime(Frames(frameIndex(1):frameIndex(end)))] - priorAnaphaseInMins;
    fit1ETimeAxis = [to, ElapsedTime(Frames(frameIndex(1):frameIndex(end)))] - nucleusFirstFrame;
    currentYSegment = [0, currentYSegment];
    
    fit1E = plot(traceFigAxes,fit1ETimeAxis,...
        currentYSegment,'-','Color','red');
    hold(traceFigAxes,'off')
    
    lineFit = 1;

%catch
%     lineFit = 0;
%     uiwait(msgbox('A line was not fitted','Key 3 was selected'));
%     fit1E = [];
%     Coefficients = [];
% end

end

% This is a subfunction for the fitInitialSlopes.m
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
% 
% Author (contact): Yang Joon Kim (yjkim90@berkeley.edu)
% This is edited from Emma Luu (emma_luu@berkeley.edu)'s code,
% fitASingleTrace.
% Created: 1/9/19
% Last Updated: 1/23/19 (YJK)
% Documented by: Yang Joon Kim (yjkim90@berkeley.edu)

%% Initializing the options
minimumLength = 3; % minimum length of trace
useDefaultTimeShift = 1;
useAnaphase = 0;
%% Getting particle information
% Amplitude ---------------------------------------------------------------
% getting frame information
[frame,~,ampIntegral3,~,~,~,~,~,~,~,~, ~,~, ampIntegralGauss3D,~]=...
    GetParticleTrace(currentParticle,...
    Particles{currentChannel},Spots{currentChannel});
currentLength = length(frame);

% performing moving average
% Also, note that this is using the 3D-Gaussian fitted spot fluorescence,
% which is ampIntegralGauss3D, if you don't have it, you should have it by
% running CheckParticleTracking with an option.
if ~sum(isnan(ampIntegralGauss3D))
    smoothedAmp = movmean(ampIntegralGauss3D,averagingLength);
else
    error('Note: Could not use the 3D guassian intensity fits so make sure you have 3D gaussian fitted spot fluorescence. Run fit3DGaussiansToAllSpots(Prefix)')
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

% initializing outputs of function
frameIndex = [];
Coefficients = [];
ErrorEstimation = struct('R',{},'df',{},'normr',{});
numberOfFramesUsedForFit = [];
 
% Assigning states to each point ------------------------------------------
% states: increase or decrease from previous point
    if currentLength > minimumLength

        % number of Frames to fit, this is given as an input
        numberOfFramesUsedForFit = length(FrameIndicesToFit);
        
        correspondingFrameIndices = FrameIndicesToFit;

        % saving the boundaries of correspondingFrameIndexes in frameIndex
        frameIndex(1,:)...
            = [correspondingFrameIndices(1) correspondingFrameIndices(end)];
        
        currentAmpSegment = smoothedAmp(correspondingFrameIndices);
        currentXSegment = currentTimeArray(correspondingFrameIndices);
        
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