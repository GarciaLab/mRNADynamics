function [Transitions,Rates] = AutoFitFluorescenceCurve3(Time,Fluo,FluoErr,GeneLength,ElongationRate)

if nargin<4
    GeneLength=6.443;
end
if nargin<5
    ElongationRate=1.54;
end

Window = 4;
TransitionSpacing = 1;
Skip = 1;

FluoErr=double(FluoErr);
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

% Transitions=Time(StartIdx:FirstFluoIdx);
% Rates=Fluo(FirstFluoIdx)/Delay * ones(size(Transitions));
% Transitions=fliplr(Time(FirstFluoIdx:-TransitionSpacing:StartIdx));
Transitions=Time(StartIdx:TransitionSpacing:end);
Rates=Fluo(FirstFluoIdx)/Delay * ones(size(Transitions));

% Find minimum rate to hit highest point in data
LeastMaxRate=max(Fluo)/Delay;
Rates(FirstFluoIdx+Skip:end)=LeastMaxRate;

for i = FirstFluoIdx+Skip:Skip:length(Time)
    i
    length(Time)
    % Index at start of window (or index 1 if we are near beginning of dataset)
    MinDataIdx=max(1,i-Window);
    MaxDataIdx=min(i+3,length(Time));
    % Decide which time values to evaluate for GOF
    TimeIn=Time(MinDataIdx:MaxDataIdx);
    FluoIn=Fluo(MinDataIdx:MaxDataIdx);
    
    % Decide which of the values in Transitions to vary
    MinFitIdx=find(Transitions>=Time(MinDataIdx),1);
    MaxFitIdx=find(Transitions<=Time(i),1,'last');
    FitIdxsAll = MinFitIdx:MaxFitIdx;
    
    % Initial rate
    InitialIdx = max(1, MinFitIdx-1);
    InitialRate = Rates(InitialIdx);
%     Rates(FitIdxsAll) = InitialRate;
    BestGOF = -1;
    BestComb = [];
    BestRates = [];

    for NCombinations = 0:length(FitIdxsAll)-1
        Combs = combnk(FitIdxsAll, NCombinations);
        % Include the max term so that 0 combination works
        for CombIdx = 1:max(size(Combs,1),1)
            if NCombinations == 0
                FitIdxs=[];
                x0=[];
                resnorm=lsqnonlinAutoFitFluorescenceCurve3(TimeIn,...
                    FluoIn,FluoErr,Delay,Transitions(1:MaxFitIdx),Rates(1:MaxFitIdx),...
                    FitIdxs,InitialIdx,InitialRate,x0);
                resnorm=sum(resnorm.^2);
                xFit=[];
            else
                FitIdxs=Combs(CombIdx,:);

                x0=Rates(FitIdxs);
                x0Lower=zeros(size(x0));
                x0Upper=inf(size(x0));

                options = optimset('MaxFunEvals',2000*length(x0),'MaxIter',5000);
                [xFit,resnorm,residual,exitflag,output,lambda,jacobian]=...
                    lsqnonlin(@(x) lsqnonlinAutoFitFluorescenceCurve3(TimeIn,...
                    FluoIn,FluoErr,Delay,Transitions(1:MaxFitIdx),Rates(1:MaxFitIdx),...
                    FitIdxs,InitialIdx,InitialRate,x),...
                    x0,x0Lower,x0Upper,options);
            end

        %     FittedTransitions=xFit(1:nSteps);
        %     FittedRates=xFit(nSteps+1:end);
        %     Transitions(FitIdxs) = FittedTransitions;
        %     Rates(FitIdxs) = FittedRates;

            GOF = resnorm+2*NCombinations;
%             GOF = resnorm/(length(FitIdxsAll)-NCombinations);
            if GOF < BestGOF || BestGOF == -1
                resnorm
%                 Rates(FitIdxs)=xFit;
                BestRates=SubIntoRates(Rates,FitIdxs,xFit,InitialIdx,InitialRate);
%                 BestRates=Rates;
%                 BestRates(FitIdxs)=xFit;
                figure(5)
                [PlotTime,PlotFluo]=IndividualTrace(Transitions(1:MaxFitIdx),BestRates(1:MaxFitIdx),Delay,80);
                plot(PlotTime,PlotFluo);
            end
        end
    end
    Rates=BestRates;
end

end

