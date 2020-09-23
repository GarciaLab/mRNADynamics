function chi2=lsqnonlinFitFluorescenceCurveNC14LinearSlopeExpDecay(TimeData,FluoData,x0)

% Fit the fluorescence curve in NC14 with initial slope and decay regime as
% well. (The decay regime is fitted to an exponential curve).

%Gives the chi square of the data to the fit form in order to do a fit with
%lsqnonlin

% Input : 
% TimeData : Time point (min)
% FluoData : Fluorescence data (AU)
% x0 : Starting conditions
TimeStart0=x0(1); % turn-on time
Rate0=x0(2); % initial rate
TimePeak0 = x0(3); % peak time (time at which the fluorescence peaks)
TimeDown0 = x0(4); % time point when the fluo starts to go down
Tau0 =x0(5); % half-life
Fluo_basal = x0(6); % basal fluorescence level (basal expression level given that the fluorescence doesn't go to zero even in late NC14)

% Output : chi2



% Delay (this one can be calculated from the formula : (gene length) /
% (elongation rate)
Delay = TimePeak0 - TimeStart0;

% genereate the fluorescence traces
FluoPrediction = zeros(1,length(TimeData));

for i=1:length(TimeData)
    if TimeData(i)<=TimeStart0
        FluoPrediction(i)=0;
    elseif (TimeData(i)>TimeStart0)&&(TimeData(i)<=TimeStart0+TimePeak0)
        FluoPrediction(i)=Rate0*(TimeData(i)-TimeStart0);
    elseif (TimeData(i)>TimeStart0+TimePeak0)&&(TimeData(i)<=TimeDown0)
        FluoPrediction(i)=Rate0*(TimePeak0 - TimeStart0); % at maximum fluorescence
    elseif (TimeData(i)>TimeDown0)
        FluoMax = Rate0*(TimePeak0 - TimeStart0);
        FluoPrediction(i)=FluoMax*exp(-(TimeData(i)-TimeDown0)/Tau0)+Fluo_basal;
%     elseif TimeData(i)>(-Rate0*Delay/RateOff0+TimeEnd0)
%         FluoPrediction(i)=0;
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
end

