function [TransitionTimes,Rates] = AutoFitFluorescenceCurveLinear(Time,Fluo,FluoErr,GeneLength,ElongationRate,MaxAllowedRate)

if nargin<4
    GeneLength=6.443;
end
if nargin<5
    ElongationRate=1.54;
end
if nargin<6
    MaxAllowedRate=-1;
end

% Parameters
TimeRes = 0.5;
Window = 9;
FluoErr=double(FluoErr);
Delay=GeneLength/ElongationRate;
Memory=floor(Delay/TimeRes);

% Interpolate time to nice values
TimeInterp=Time(1):TimeRes:Time(end);
Fluo=interp1(Time,Fluo,TimeInterp);

TransitionIdxs=1:length(Fluo);
Rates=zeros(size(TransitionIdxs));

for i = 1:length(Fluo)-1
    % Just a place for a conditional breakpoint
%     if i==28
%         abcdef=1;
%     end

    % Decide which of the values in Transitions to vary
    MinFitIdx=i;
    MaxFitIdx=min(i+Window-1,length(TransitionIdxs));
    N=MaxFitIdx-MinFitIdx+1;
    
    % Index at start of window (or index 1 if we are near beginning of dataset)
    MinDataIdx=i;
    MaxDataIdx=min(i+Window-1,length(TransitionIdxs));

    % Decide which time values to evaluate for GOF
    FluoIn=Fluo(MinDataIdx:MaxDataIdx);
    
    % Previous values
    PreIdx=max(1,i-Memory);
    XPrevious=Rates(PreIdx:i-1)*TimeRes;
    
    % Initialize best
    BestGOF = -1;
    BestComb = [];
    BestRates = [];

    for NCombinations = 0:N-2%min(1,N-2)%N-2
        Combs = combnk(1:N, NCombinations);
        % Include the max term so that 0 combination works
        for CombIdx = 1:max(1,size(Combs,1))
            if NCombinations==0
                FitIdxs=[];
                [XOut,FOut,Chi2]=lsqlinAutoFitFluorescenceCurve(FluoIn,FitIdxs,XPrevious,Delay,TimeRes,FluoErr);
            else
                FitIdxs=Combs(CombIdx,:);
                [XOut,FOut,Chi2]=lsqlinAutoFitFluorescenceCurve(FluoIn,FitIdxs,XPrevious,Delay,TimeRes,FluoErr);
            end
            FitRates=XOut/TimeRes;

            % Evaluate chi2
            TempRates=Rates;
            TempRates(MinFitIdx:MaxFitIdx)=FitRates;

            AIC = Chi2+2*NCombinations;
%             GOF = AIC;
%             GOF = AIC + 4*NCombinations;
            GOF = AIC + 2*NCombinations*(NCombinations+1)/(N-NCombinations-1);
%             GOF = Chi2/(N-NCombinations);
%             GOF = Chi2;
            if GOF < BestGOF || BestGOF == -1
%                 chi2
%                 BestRates=SubIntoRates(Rates,FitIdxs,xFit,InitialIdx,InitialRate);
                BestGOF=GOF;
                BestRates=TempRates;
%                 figure(5)
%                 clf;
%                 hold on;
%                 [PlotTime,PlotFluo]=IndividualTrace(TimeInterp(1:MaxFitIdx)-0.5,BestRates(1:MaxFitIdx),Delay,TimeInterp(MaxFitIdx));
%                 plot(PlotTime,PlotFluo,'b');
%                 plot(TimeInterp(MinFitIdx:MaxFitIdx),FOut,'r')
%                 plot(TimeInterp(1:MaxFitIdx),Fluo(1:MaxFitIdx),'g')
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

