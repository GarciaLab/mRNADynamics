function [TimeRange,Fluorescence]=FluorescenceCurveLinearSlope(ncLength,TimeStart,Rate,Delay)

%This function outputs the initial rise in a fluorescence trace given a transcription
%rate, a start time, a delay time and an end time. Note that the end time
%isn't used for this script.




% TimeStart=5;
% TimeEnd=12;
% 
% Rate=4E3;     %Rate per minute

TimeRange=linspace(0,ncLength,1000);

RateFrame=Rate*mean(diff(TimeRange));        %Rate per frame
Fluorescence=zeros(size(TimeRange));



for j=2:length(TimeRange)
    %Production part
    if (TimeRange(j)>TimeStart)
        Fluorescence(j)=Fluorescence(j-1)+RateFrame;
    end
    
    %Termination
    if ((TimeRange(j)-Delay*3)>=TimeStart)&&(Fluorescence(j-1)>1E-10)
        Fluorescence(j)=nan;
    end
end
