function [frameIndex,Coefficients,ErrorEstimation,numberOfParticlesUsedForFit] = ...
    fitASingleTrace(currentParticle,Particles,Spots,currentChannel,...
    ElapsedTime,nuclearCycleBoundaries,correspondingNCInfo,...
    averagingLength,varargin)
%% Initializing the options
skipSavingTraces = 0;
initialOnly = 0;
minimumLength = 3; % minimum length of trace

%% Checking varargin
if length(varargin)
    for i=1:length(varargin)
        if strcmpi(varargin{i},'initialOnly')
            initialOnly = 1;
        elseif strcmpi(varargin{i},'minimumLength')
            minimumLength = varargin{i+1};
        elseif strcmpi(varargin{i},'skipSavingTraces')
            skipSavingTraces = 1;
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
% adjusting frameRange to have time 0 be the start of the
% first nuclear cycle the particle appears in
currentTimeArray = currentTimeArray - timeOfFirstNC;


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