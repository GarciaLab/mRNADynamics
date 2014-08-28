function [Integral,SDIntegral]=...
    IntegrateParticle(Particle,FrameRange,ElapsedTime)

%Integrates a particle over a frame range. If FrameRange is empty it uses
%the whole range.

if isempty(FrameRange)
    FrameRange=Particle.Frame;
end

%Based on the FrameRange create a filter
FrameFilter=ismember(Particle.Frame,FrameRange);

if sum(FrameFilter)>1
    Integral=trapz(ElapsedTime(Particle.Frame(FrameFilter)),...
        Particle.Fluo(FrameFilter));

    %Estimate the error
    if length(FrameRange)==2
        SDIntegral=...
            (ElapsedTime(FrameRange(2))-ElapsedTime(FrameRange(1)))/2*...
            Particle.FluoError*sqrt(2);
    else
        ErrorTemp=[];
        %Calculate the error of the inner points
        for k=2:(length(FrameRange)-1)
             ErrorTemp(k)=...
                 (ElapsedTime(FrameRange(k+1))-ElapsedTime(FrameRange(k-1)))*...
                 Particle.FluoError;
        end
        
        ErrorTemp(1)=(ElapsedTime(FrameRange(2))-ElapsedTime(FrameRange(1)))/2*...
            Particle.FluoError;
        ErrorTemp(length(FrameRange))=(ElapsedTime(FrameRange(end))-ElapsedTime(FrameRange(end-1)))/2*...
            Particle.FluoError;

        %Now, add it all up
        SDIntegral=sqrt(sum(ErrorTemp.^2));
        
    end
else
    Integral=0;
    SDIntegral=0;
end

