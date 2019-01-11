function [lineFit, Coefficients, fit1E, Particles] =...
    fitInitialSlope(CurrentParticle, Particles, Spots, CurrentChannel, schnitzcells, ...
    ElapsedTime, anaphaseInMins, correspondingNCInfo, traceFigAxes, Frames, anaphase, ...
    averagingLength, FramesToFit, FrameIndicesToFit)
%fitInitialSlope Summary of this function goes here
%   Detailed explanation goes here

% Input parameters that should be defined :
% AverageLength, Time window for fitting(adjustable), 
% Use GUI for defining the inputs, as well as repeating the fitting until
% it's approved.

%while (~dd==13) % 13 is same as enter key
    % Plug in inputs defined in GUI, 
    % averagingLength, FramesToFit
    
%     averagingLength = 1; % default
%     % Define the frames to fit
%     [X,Y] = ginput(2); % pick two points (left, and right)
%     pos1 = find((Frames-X(1)).^2 == min((Frames-X(1)).^2))
%     pos2 = find((Frames-X(2)).^2 == min((Frames-X(2)).^2))
%     FramesToFit = [pos1:pos2];
     try
        lineFit = 1; % this should go after the approval/rejection
        
        % currently shifted by the first frame of the assigned nucleus
        [frameIndex,Coefficients,ErrorEstimation,nFramesForFit] = ...
            fitASingleTraceManual(CurrentParticle,Particles,Spots,CurrentChannel,...
            schnitzcells,ElapsedTime,anaphaseInMins,correspondingNCInfo,...
            averagingLength, FramesToFit,FrameIndicesToFit);


        % plotting the fitted line
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
        %              currentXSegment = ElapsedTime(Frames(frameIndex(1):frameIndex(end)))-priorAnaphaseInMins;
        currentXSegment = ElapsedTime(Frames(frameIndex(1):frameIndex(end)))-nucleusFirstFrame;
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
        to = -Coefficients(2) / Coefficients(1) + nucleusFirstFrame;
        %              fit1ETimeAxis = [to, ElapsedTime(Frames(frameIndex(1):frameIndex(end)))] - priorAnaphaseInMins;
        fit1ETimeAxis = [to, ElapsedTime(Frames(frameIndex(1):frameIndex(end)))] - nucleusFirstFrame;
        currentYSegment = [0, currentYSegment];
        %              fit1ETimeAxis = ElapsedTime(Frames(frameIndex(1):frameIndex(end)));
        fit1E = plot(traceFigAxes,fit1ETimeAxis,...
            currentYSegment,'-','Color','red');
        hold(traceFigAxes,'off')


    catch
        lineFit = 0;
        uiwait(msgbox('A line was not fitted','Key 3 was selected'));
        fit1E = [];
        Coefficients = [];
     end
%end

    % save the fitted values (Slope and Time on) in Particles.mat
    % Save only after the approval, which is pressing the 'enter' key
    if ~isempty(Coefficients)
        singleTraceLoadingRate = Coefficients(1,1); %au/min
        if singleTraceLoadingRate >= 0 %some easy quality control
            singleTraceTimeOn = roots(Coefficients(1,:));
            Particles{CurrentChannel}(CurrentParticle).fittedSlope =  singleTraceLoadingRate;
            Particles{CurrentChannel}(CurrentParticle).fittedTON =  singleTraceTimeOn;
        else     
            Particles{CurrentChannel}(CurrentParticle).fittedSlope =  NaN;
            Particles{CurrentChannel}(CurrentParticle).fittedTON =  NaN; 
        end
    else
        Particles{CurrentChannel}(CurrentParticle).fittedSlope =  NaN;
        Particles{CurrentChannel}(CurrentParticle).fittedTON =  NaN;   
    end
end

