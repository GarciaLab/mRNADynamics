function chi=lsqnonlinAutoFitFluorescenceCurve(TimeData,FluoData,FluoErr,Delay,Transitions,Rates,FitIdxs,x0)

%Gives the chi square of the data to the fit form in order to do a fit with
%lsqnonlin

Rates(FitIdxs)=x0;

%Get the predicted shape
[TimePrediction,FluoPrediction]=IndividualTrace(Transitions,Rates,Delay,max(TimeData)+1);
FluoPrediction=interp1(TimePrediction,FluoPrediction,TimeData);

%Calculate chi^2
try
    chi=(FluoData-FluoPrediction)./FluoErr;
catch
    1+1;
end
