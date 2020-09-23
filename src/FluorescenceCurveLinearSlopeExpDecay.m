function [TimeRange,Fluorescence]=FluorescenceCurveLinearSlopeExpDecay(ncLength, xFit)

%This function outputs a standard fluorescent trace given a transcription
%rate, a start time, a delay time and an end time.

% TimeStart=5;
% TimeEnd=12;
% 
% Rate=4E3;     %Rate per minute

% xFit is a fitted value of coefficients that recaps the characteristics of
% the trace.
TimeStart = xFit(1);
Rate = xFit(2);
TimePeak  = xFit(3);
TimeDown = xFit(4);
Tau = xFit(5);
FluoBasal = xFit(6);

TimeRange=linspace(0,ncLength,1000);

RateFrame=Rate*mean(diff(TimeRange));        %Rate per frame
Fluorescence=zeros(size(TimeRange));

FluoMax =[];

for j=2:length(TimeRange)
    %Production part
    if (TimeRange(j)>TimeStart) && (TimeRange(j)<=TimePeak)
        Fluorescence(j)=Fluorescence(j-1)+RateFrame;
    elseif (TimeRange(j)>TimePeak) && (TimeRange(j) <= TimeDown) 
        Fluorescence(j)=Fluorescence(j-1);
        FluoMax = Fluorescence(j);
    elseif TimeRange(j) > TimeDown 
        Fluorescence(j)=FluoMax*exp(-(TimeRange(j)-TimeDown)/Tau) + FluoBasal;
    end
    
end
