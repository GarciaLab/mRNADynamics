function [lineFit, Coefficients, fit1E] =...
    fitLine(CurrentParticle, Particles, Spots, CurrentChannel, schnitzcells, ...
    ElapsedTime, anaphaseInMins, correspondingNCInfo, traceFigAxes, Frames)
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
    priorAnaphaseInMins = anaphaseInMins(ncPresent(1)-8);
    nucleusFirstFrame = ElapsedTime(...
        schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).frames(1));
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
end
end

