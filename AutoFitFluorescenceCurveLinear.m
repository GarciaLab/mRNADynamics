function [TransitionTimes,Rates] = AutoFitFluorescenceCurveLinear(Time,Fluo,FluoErr,GeneLength,ElongationRate)

if nargin<4
    GeneLength=6.443;
end
if nargin<5
    ElongationRate=1.54;
end

% Parameters
TimeRes = 0.5;
Window = 8;
FluoErr=double(FluoErr);
Delay=GeneLength/ElongationRate;
Memory=floor(Delay/TimeRes);

% Interpolate time to nice values
TimeInterp=Time(1):TimeRes:Time(end);
Fluo=interp1(Time,Fluo,TimeInterp);

TransitionIdxs=1:length(Fluo);
Rates=zeros(size(TransitionIdxs));

for i = 1:length(Fluo)-1
    i
    % Decide which of the values in Transitions to vary
    MinFitIdx=i;
    MaxFitIdx=min(i+Window,length(TransitionIdxs));
    N=MaxFitIdx-MinFitIdx+1;
    
    % Index at start of window (or index 1 if we are near beginning of dataset)
    MinDataIdx=i;
    MaxDataIdx=min(i+Window,length(TransitionIdxs));

    % Decide which time values to evaluate for GOF
    FluoIn=Fluo(MinDataIdx:MaxDataIdx);
    
    % Previous values
    PreIdx=max(1,i-Memory+1);
    XPrevious=Rates(PreIdx:i-1)*TimeRes;
    
    % Initialize best
    BestGOF = -1;
    BestComb = [];
    BestRates = [];

    for NCombinations = 0:N-2
        Combs = combnk(1:N, NCombinations);
        % Include the max term so that 0 combination works
        for CombIdx = 1:max(1,size(Combs,1))
            if NCombinations==0
                FitIdxs=[];
                [XOut,FOut,Chi2]=lsqlinAutoFitFluorescenceCurve(FluoIn,FitIdxs,XPrevious,Memory,FluoErr);
            else
                FitIdxs=Combs(CombIdx,:);
                T=FitIdxs;
                [XOut,FOut,Chi2]=lsqlinAutoFitFluorescenceCurve(FluoIn,T,XPrevious,Memory,FluoErr);
            end
            FitRates=XOut/TimeRes;

            % Evaluate chi2
            TempRates=Rates;
            TempRates(MinFitIdx:MaxFitIdx)=FitRates;

            AIC = Chi2+2*NCombinations;
%             GOF = AIC;
%             GOF = AIC + 6*NCombinations;
            GOF = AIC + 2*NCombinations*(NCombinations+1)/(N-NCombinations-1);
%             GOF = Chi2/(N-NCombinations);
            if GOF < BestGOF || BestGOF == -1
%                 chi2
%                 BestRates=SubIntoRates(Rates,FitIdxs,xFit,InitialIdx,InitialRate);
                BestGOF=GOF;
                BestRates=TempRates;
                figure(5)
                [PlotTime,PlotFluo]=IndividualTrace(TimeInterp(1:MaxFitIdx),2*BestRates(1:MaxFitIdx),Delay,80);
                plot(PlotTime,PlotFluo);
            end
        end
    end
    Rates=BestRates;
end

% Make some adjustments due to assumptions that simplified calculations

% Subtract TimeRes from TransitionTimes because we assumed Rate(i) had its
% first finite value at Fluo(i), should be Fluo(i+1)
TransitionTimes=TimeInterp-TimeRes;

end

