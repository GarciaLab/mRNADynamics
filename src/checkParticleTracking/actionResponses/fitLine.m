function [lineFit, Coefficients, fit1E, Particles] =...
    fitLine(CurrentParticle, Particles, Spots, CurrentChannel, schnitzcells, ...
    ElapsedTime, anaphaseInMins, correspondingNCInfo, traceFigAxes, Frames, anaphase)
%FITLINE Summary of this function goes here
%   Detailed explanation goes here

averagingLength = 3;

 try
    lineFit = 1;
    % currently shifted by the first frame of the assigned nucleus
    [frameIndex,Coefficients,ErrorEstimation,nParticlesForFit] = ...
        fitASingleTrace(CurrentParticle,Particles,Spots,CurrentChannel,...
        schnitzcells,ElapsedTime,anaphaseInMins,correspondingNCInfo,...
        averagingLength,'initialOnly','skipSavingTraces');


    % plotting the fitted line
    ncPresent = unique(correspondingNCInfo(Frames));
    % below subtracts 8 because the first element corresponds to nc 9
%     priorAnaphaseInMins = anaphaseInMins(ncPresent(1)-8);
    priorAnaphase = anaphase(ncPresent(1)-8); %frame
    if ~isempty(Particles{CurrentChannel}(CurrentParticle).Nucleus)
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
    
    % save the fitted values (Slope and Time on) in Particles.mat
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

catch
    lineFit = 0;
    uiwait(msgbox('A line was not fitted','Key 3 was selected'));
    fit1E = [];
    Coefficients = [];
end
end

