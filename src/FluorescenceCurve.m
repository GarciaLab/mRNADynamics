function [TimeRange,Fluorescence]=FluorescenceCurve(ncLength,TimeStart,TimeEnd,Rate,Delay)

%This function outputs a standard fluorescent trace given a transcription
%rate, a start time, a delay time and an end time.




% TimeStart=5;
% TimeEnd=12;
% 
% Rate=4E3;     %Rate per minute

TimeRange=linspace(0,ncLength,1000);

RateFrame=Rate*mean(diff(TimeRange));        %Rate per frame
Fluorescence=zeros(size(TimeRange));



for j=2:length(TimeRange)
    %Production part
    if (TimeRange(j)>TimeStart)&(TimeRange(j)<TimeEnd)
        Fluorescence(j)=Fluorescence(j-1)+RateFrame;
    elseif TimeRange(j)>=TimeEnd
        Fluorescence(j)=Fluorescence(j-1);
    end
    
    %Termination
    if ((TimeRange(j)-Delay)>=TimeStart)&(Fluorescence(j-1)>1E-10)
        Fluorescence(j)=Fluorescence(j)-RateFrame;
    end
end
