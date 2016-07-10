function [Slope,SDSlope]=TraceDerivative(CompiledParticle,ElapsedTime,TimeWindow)

%Calculates the derivative of a trace as a function of time for a certain
%time window.


MinNumber=3;        %Minimum number within a window so that we can calculate
                    %a derivative

if (length(CompiledParticle.Frame)-TimeWindow)>0
    for i=1:(length(CompiledParticle.Frame)-TimeWindow)

        %Check that there MinNumber of frames within this window

        if sum(ismember(CompiledParticle.Frame,...
            [CompiledParticle.Frame(i):CompiledParticle.Frame(i)+(TimeWindow-1)]))>=MinNumber

            FramesRangeFilter=ismember(CompiledParticle.Frame,...
                [CompiledParticle.Frame(i):CompiledParticle.Frame(i)+(TimeWindow-1)]);



            [a, b, sigma_a, sigma_b] = york_fit(ElapsedTime(CompiledParticle.Frame(FramesRangeFilter)),...
                CompiledParticle.Fluo(FramesRangeFilter),...
                ones(1,sum(FramesRangeFilter))*mean(diff(ElapsedTime))/2,...
                ones(1,sum(FramesRangeFilter))*CompiledParticle.FluoError);


            Slope(i)=b;
            SDSlope(i)=sigma_b;

        else
            Slope(i)=nan;
            SDSlope(i)=nan;
        end
    end
else
    Slope=nan;
    SDSlope=nan;
end
    


                    
                    