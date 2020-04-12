function chi2=lsqnonlinFitFluorescenceCurve(TimeData,FluoData,Delay,ncLength,x0)

%Gives the chi square of the data to the fit form in order to do a fit with
%lsqnonlin

%Starting conditions
TimeStart0=x0(1);
TimeEnd0=x0(2);
Rate0=x0(3);


for i=1:length(TimeData)
    if TimeData(i)<=TimeStart0
        FluoPrediction(i)=0;
    elseif (TimeData(i)>TimeStart0)&(TimeData(i)<=TimeStart0+Delay)
        FluoPrediction(i)=Rate0*(TimeData(i)-TimeStart0);
    elseif (TimeData(i)>TimeStart0+Delay)&(TimeData(i)<=TimeEnd0)
        FluoPrediction(i)=Rate0*Delay;
    elseif (TimeData(i)>TimeEnd0)&(TimeData(i)<=TimeEnd0+Delay)
        FluoPrediction(i)=Rate0*(-TimeData(i)+Delay+TimeEnd0);
    elseif TimeData(i)>TimeEnd0+Delay
        FluoPrediction(i)=0;
    end
end



%[TimeRange,Fluorescence]=FluorescenceCurve(ncLength,TimeStart0,TimeEnd0,Rate0,Delay);


%Interpolate the values of the generated curve for the actual time points
%we have
%InterpData=pchip(TimeRange,Fluorescence,TimeData);

%NanFilter=~isnan(FluoData);
%chi2=(FluoData(NanFilter)-InterpData(NanFilter)').^2;

try
    chi2=(FluoData-FluoPrediction').^2;
catch
    1+1;
end

