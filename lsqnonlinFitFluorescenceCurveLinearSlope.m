function chi2=lsqnonlinFitFluorescenceCurveLinearSlope(TimeData,FluoData,Delay,ncLength,x0)

%Gives the chi square of the data to the fit form in order to do a fit with
%lsqnonlin, using a linear fit to the initial fluorescence rise. Pass to
%this function only the data points to be used in this initial fluorescence
%rise. Note that since Matlab's lsqnonlin function takes the chi square
%value as the unsquared form and performs the square itself, this function
%only returns the square root of the chi square, i.e. the residual vector.

%Starting conditions
TimeStart0=x0(1); %Time of start of fit window
Rate0=x0(2); %Initiation rate


for i=1:length(TimeData)
    if TimeData(i)<=TimeStart0
        FluoPrediction(i)=0;
    elseif (TimeData(i)>TimeStart0)
        FluoPrediction(i)=Rate0*(TimeData(i)-TimeStart0);
    end
end



%[TimeRange,Fluorescence]=FluorescenceCurve(ncLength,TimeStart0,TimeEnd0,Rate0,Delay);


%Interpolate the values of the generated curve for the actual time points
%we have
%InterpData=pchip(TimeRange,Fluorescence,TimeData);

%NanFilter=~isnan(FluoData);
%chi2=(FluoData(NanFilter)-InterpData(NanFilter)').^2;
try
    chi2=(FluoData-FluoPrediction');
catch
    1+1;
end

