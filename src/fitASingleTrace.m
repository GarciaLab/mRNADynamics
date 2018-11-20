function [frameIndex,Coefficients,ErrorEstimation,numberOfParticlesUsedForFit] = ...
    fitASingleTrace(currentParticle,Particles,Spots,currentChannel,...
    schnitzcells,ElapsedTime,nuclearCycleBoundaries,correspondingNCInfo,...
    averagingLength,varargin)
% CompileParticles(varargin)
%
% DESCRIPTION
% Fits a piece-wise function of lines to a fluorescent trace.
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
% averagingLength : length of averaging done by the movmean function. 
%
% OPTIONS
% 'initialOnly' : only spends time fitting the polymerase 
%                 loading rate of the trace
% 'minimumLength' : adjusts the threshold number of total number of points
%                   in a trace that will be fitted, not inclusive.
%                  (default = 3 points)
% 'skipSavingTraces' : skips saving the traces into a folder (meant for
%                      CompileParticles
% 'useAnaphase' : uses the users values of the anaphase to set time 0 to be
%                 the start of the nuclear cycle of the particle. 
% 
% Author (contact): Emma Luu (emma_luu@berkeley.edu)
% Created: 11/7/18
% Last Updated: 11/20/18 (EL)
% Documented by: Emma Luu (emma_luu@berkeley.edu)

%% Initializing the options
skipSavingTraces = 0;
initialOnly = 0;
minimumLength = 3; % minimum length of trace
useDefaultTimeShift = 1;
useAnaphase = 0;

%% Checking varargin
if ~isempty(varargin)
    for i=1:length(varargin)
        if strcmpi(varargin{i},'initialOnly')
            initialOnly = 1;
        elseif strcmpi(varargin{i},'minimumLength')
            minimumLength = varargin{i+1};
        elseif strcmpi(varargin{i},'skipSavingTraces')
            skipSavingTraces = 1;
        elseif strcmpi(varargin{i},'useAnaphase')
            useAnaphase = 1;
            useDefaultTimeShift = 0;
        end
    end
end

%% Getting particle information
% Amplitude ---------------------------------------------------------------
% getting frame information
[frame,~,ampIntegral3,~,~,~,~,~,~,~,~,~,~]=...
    GetParticleTrace(currentParticle,...
    Particles{currentChannel},Spots{currentChannel});
currentLength = length(frame);

% performing moving average
smoothedAmp = movmean(ampIntegral3,averagingLength);

% Time --------------------------------------------------------------------
% getting the corresponding time of the trace
currentTimeArray = ElapsedTime(frame); % Units to seconds
ncPresent = unique(correspondingNCInfo(frame));
% below subtracts 8 because the first element correspond to nc 9
timeOfFirstNC = nuclearCycleBoundaries(ncPresent(1)-8);
nucleusFirstTime = ElapsedTime(...
    schnitzcells(Particles{currentChannel}(currentParticle).Nucleus).frames(1));
% adjusting frameRange to have time 0 be the start of the
% first nuclear cycle the particle appears in or the first time point of
% the nucleus it is assigned to
if useDefaultTimeShift
    currentTimeArray = currentTimeArray - nucleusFirstTime;
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
numberOfParticlesUsedForFit = [];
 
% Assigning states to each point ------------------------------------------
% states: increase or decrease from previous point
if currentLength > minimumLength
    pointGroupNumbers = groupingPoints(currentTimeArray,smoothedAmp);
    
    % creating the lines of best fit for each group -----------------------
    if initialOnly
        numberOfGroups = 1;
    else
        numberOfGroups = max(pointGroupNumbers);
    end
    
    numberOfParticlesUsedForFit = zeros(1,numberOfGroups);
    for currentGroup = 1:numberOfGroups
        correspondingFrameIndices = find(pointGroupNumbers == currentGroup);
        if currentGroup > 1
            % using the previous point as part of the line...
            correspondingFrameIndices = [min(correspondingFrameIndices)-1 correspondingFrameIndices];
        end
        % saving the boundaries of correspondingFrameIndexes in frameIndex
        frameIndex(currentGroup,:)...
            = [correspondingFrameIndices(1) correspondingFrameIndices(end)];
        
        currentAmpSegment = smoothedAmp(correspondingFrameIndices);
        currentXSegment = currentTimeArray(correspondingFrameIndices);
        
        [Coefficients(currentGroup,:),...
            ErrorEstimation(currentGroup)]...
            =  polyfit(currentXSegment,currentAmpSegment,1);
        
        numberOfParticlesUsedForFit(currentGroup) = ...
                sum(pointGroupNumbers==currentGroup);% Saving the number of particles included in group
            
        if currentGroup == 1
            currentYSegment = ...
                polyval(Coefficients(currentGroup,:),currentXSegment); 
                % this is for the predicted line
            denominator = sum((currentAmpSegment - mean(currentAmpSegment)).^2);
            normOfResiduals = ErrorEstimation(currentGroup).normr;
            RSquared = 1 - (normOfResiduals^2)/denominator;
            errorArray = ones(1,length(currentXSegment)).*...
                normOfResiduals./numberOfParticlesUsedForFit(currentGroup); %EL normalized by number of points included
            
            if ~skipSavingTraces
                % plot and save  the first fitted line with the
                % smoothed trace
                
                plot(currentTimeArray,smoothedAmp,'.','MarkerSize',20,...
                    'DisplayName','Avg amp'); % plot the smoothed trace and fitted line
                hold on
                % plotting the fitted line
                plot(currentXSegment, ...
                    currentYSegment,'DisplayName','fittedLine');
                errorbar(currentXSegment,currentYSegment,errorArray,'o');
                legend({'Avg Amp', ...
                    ['Fitted Line Slope: ' num2str(Coefficients(currentGroup,1))],...
                    'Norm of Residuals'},'Location','Best')
                title(['Current Particle : ' num2str(currentParticle)...
                    '   R squared :' num2str(RSquared)])
                xlabel('Time after start of nc (minutes)')
                ylabel('Integrated Intensity (A.U.)')
                hold off
                
            end
        end
        
        
    end
end



end