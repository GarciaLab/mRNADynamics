function [TimeInterp,Rates] = AutoFitFluorescenceCurve(Time,Fluo,FluoErr,GeneLength,ElongationRate)

if nargin<4
    GeneLength=6.443;
end
if nargin<5
    ElongationRate=1.54;
end

FluoErr=double(FluoErr);

TimeRes = 0.5;
Window = 4;
Delay=GeneLength/ElongationRate;
MemoryLength=ceil(Delay/TimeRes);

% Interpolate time to nice values
TimeInterp=min(Time):TimeRes:max(Time);
Fluo=interp1(Time,Fluo,TimeInterp);
N=length(Fluo);

% Find first point
[~,FirstFluoIdx]=find(Fluo>0,1);
% StartIdx=max(1,FirstFluoIdx-MemoryLength);
% Transitions=Time(StartIdx:FirstFluoIdx);
TransitionIdxs=1:N;
% Rates=zeros(size(TransitionIdxs));
Rates=Fluo(FirstFluoIdx)/Delay * ones(size(Fluo));

% Find minimum rate to hit highest point in data
% LeastMaxRate=max(Fluo)/Delay;

for i = 1:N-1
    i
%     Transitions(end+1)=Time(i);
%     Rates(end+1)=LeastMaxRate;
    % Index at start of window (or index 1 if we are near beginning of dataset)
%     MinIdx=max(1,i-Window);
    MinFitIdx=i;
    MaxFitIdx=min(i+Window,N);
    
    % Decide which of the values in Transitions to vary
%     MinFitIdx=find(Transitions>=Time(MinIdx),1);
%     MaxFitIdx=find(Transitions<=Time(i),1,'last');
    FitIdxs = MinFitIdx:MaxFitIdx;

%     x0=[Transitions(FitIdxs),Rates(FitIdxs)];
    x0=Rates(FitIdxs);
    %Create lower and upper bounds for the rates
    x0Lower=zeros(size(x0));
    x0Upper=inf(size(x0));
    
    % Decide which time values to evaluate chi2 against
    TimeIn=TimeInterp(MinFitIdx:MaxFitIdx);
    FluoIn=Fluo(MinFitIdx:MaxFitIdx);

    options = optimset('MaxFunEvals',2000*length(x0),'MaxIter',5000);
    [xFit,resnorm,residual,exitflag,output,lambda,jacobian]=...
        lsqnonlin(@(x) lsqnonlinAutoFitFluorescenceCurve(TimeIn,...
        FluoIn,FluoErr,Delay,TimeInterp(1:MaxFitIdx),Rates(1:MaxFitIdx),FitIdxs,x),...
        x0,x0Lower,x0Upper,options);

    Rates(FitIdxs)=xFit;
end

end

