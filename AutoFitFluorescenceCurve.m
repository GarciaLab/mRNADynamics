function [Transitions,Rates] = AutoFitFluorescenceCurve(Time,Fluo,FluoErr,GeneLength,ElongationRate)

if nargin<4
    GeneLength=6.443;
end
if nargin<5
    ElongationRate=1.54;
end

FluoErr=double(FluoErr);

Window = 4;
TimeRes = mean(diff(Time));
TimeRes = 0.5;
Delay=GeneLength/ElongationRate;
MemoryLength=ceil(Delay/TimeRes);

% Should never set window above memory length
if Window > MemoryLength
    Window=MemoryLength;
end

% Interpolate time to nice values
TimeOld=Time;
Time=min(Time):TimeRes:max(Time);
Fluo=interp1(TimeOld,Fluo,Time);

% Find first point
[~,FirstFluoIdx]=find(Fluo>0,1);
StartIdx=max(1,FirstFluoIdx-MemoryLength);
% TransitionIdxs=StartIdx:FirstFluoIdx;
Transitions=Time(StartIdx:FirstFluoIdx);
Rates=Fluo(FirstFluoIdx)/Delay * ones(size(Transitions));

% Find minimum rate to hit highest point in data
LeastMaxRate=max(Fluo)/Delay;

for i = FirstFluoIdx+1:length(Time)
    i
    Transitions(end+1)=Time(i);
    Rates(end+1)=LeastMaxRate;
    % Index at start of window (or index 1 if we are near beginning of dataset)
    MinIdx=max(1,i-Window);
    
    % Decide which of the values in Transitions to vary
    MinFitIdx=find(Transitions>=Time(MinIdx),1);
    MaxFitIdx=find(Transitions<=Time(i),1,'last');
    FitIdxs = MinFitIdx:MaxFitIdx;

%     x0=[Transitions(FitIdxs),Rates(FitIdxs)];
    x0=Rates(FitIdxs);
    %Create lower and upper bounds for the rates
    x0Lower=zeros(size(x0));
    x0Upper=inf(size(x0));
    nSteps=length(FitIdxs);
    
    % Decide which time values to evaluate chi2 against
    MaxIdx=min(i+3,length(Time));
    TimeIn=Time(MinIdx:MaxIdx);
    FluoIn=Fluo(MinIdx:MaxIdx);

    options = optimset('MaxFunEvals',2000*length(x0),'MaxIter',5000);
%     TimeIn
%     FluoIn
%     Transitions
%     Rates
%     MinIdx
%     MinFitIdx
%     MaxFitIdx
%     FitIdxs
    [xFit,resnorm,residual,exitflag,output,lambda,jacobian]=...
        lsqnonlin(@(x) lsqnonlinAutoFitFluorescenceCurve(TimeIn,...
        FluoIn,FluoErr,Delay,Transitions,Rates,FitIdxs,nSteps,x),...
        x0,x0Lower,x0Upper,options);

%     FittedTransitions=xFit(1:nSteps);
%     FittedRates=xFit(nSteps+1:end);
    
%     Transitions(FitIdxs) = FittedTransitions;
%     Rates(FitIdxs) = FittedRates;
    Rates(FitIdxs)=xFit;
end

end

