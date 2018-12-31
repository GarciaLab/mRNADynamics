function [CompiledParticles, fittedLineEquations] = fitShapesToTraces(Prefix, ...
    Particles, schnitzcells, FrameInfo, ElapsedTime, CompiledParticles)
%FITSHAPESTOTRACES Summary of this function goes here
%   Detailed explanation goes here

%     try
fittedLineEquations = fitSingleTraces(Prefix,Particles,Spots,schnitzcells,FrameInfo,ElapsedTime);
ChN = 1; %only supports one channel for now
for i = 1:length(CompiledParticles{ChN})
    if ~isempty(fittedLineEquations(i).Coefficients)
        singleTraceLoadingRate = fittedLineEquations(i).Coefficients(1,1); %au/min
        if singleTraceLoadingRate >= 0 %some easy quality control
            singleTraceTimeOn = roots(fittedLineEquations(i).Coefficients(1,:));
            CompiledParticles{ChN}(i).singleTraceLoadingRate = singleTraceLoadingRate;
            CompiledParticles{ChN}(i).singleTraceTimeOn = singleTraceTimeOn;
        else     
            CompiledParticles{ChN}(i).singleTraceLoadingRate = NaN;
            CompiledParticles{ChN}(i).singleTraceTimeOn = NaN;     
        end
    else
        CompiledParticles{ChN}(i).singleTraceLoadingRate = NaN;
        CompiledParticles{ChN}(i).singleTraceTimeOn = NaN;      
    end
end
%     catch
%     %Talk to EL if you need this.     
%     end
end

