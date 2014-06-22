function [TransitionTimes,Rates] = AutoFitFluorescenceCurveLinear(Time,Fluo,FluoError,GeneLength,ElongationRate,MaxAllowedRate)
% AutoFitFluorescenceCurveLinear: Fit a set of underlying rates of
% initiation to given fluorescence data using the multi-rate fitting
% algorithm.
% Optimality is determined by minimizing the adjusted Akaike information
% criteria in a sliding window
%
% Time : Time of input data
% Fluo : Fluorescence of input data
% FluoError : Uncertainty of fluorescence (a scalar, for now)
% GeneLength : Length of the gene in bp
% ElongationRate : Rate of elongation in bp/min
% MaxAllowedRate : Maximum rate to use in fitting (CURRENTLY UNIMPLEMENTED)

%%% DEFAULT INPUTS
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
TimeRes = 0.5; % grid for rate changes
Window = 9; % number of timesteps to look ahead while fitting current rate
Delay=GeneLength/ElongationRate; % Elongation time (in minutes, not indices)
Memory=floor(Delay/TimeRes); % Elongation time (in indices, rounded down)

FluoError=double(FluoError); % FluoError stored as single, needs to be double
% Interpolate data to grid
TimeInterp=Time(1):TimeRes:Time(end);
Fluo=interp1(Time,Fluo,TimeInterp);

% Rates : list of all fitted rates so far
Rates=zeros(size(TimeInterp));

for i = 1:length(Fluo)-1
    % Window of rates to vary
    MinFitIdx=i;
    MaxFitIdx=min(i+Window-1,length(TimeInterp));
    N=MaxFitIdx-MinFitIdx+1;

    % Fit Fluo within window
    FluoIn=Fluo(MinFitIdx:MaxFitIdx);
    
    % Account for previous rates within memory time
    PreIdx=max(1,i-Memory);
    % Convert rates to X
    XPrevious=Rates(PreIdx:i-1)*TimeRes;
    
    % Initialize best tracking
    BestGOF = -1;
    BestComb = [];
    BestRates = [];

    % To use AICc, can have at most floor((N-2)/2) free rates
    for NFree = 0:floor((N-2)/2)
        Combs = combnk(1:N, NFree);
        % The max term in the next line is necessary for NFree=0
        for CombIdx = 1:max(1,size(Combs,1))
            %%% Do the fit with these rates free
            if NFree==0
                FitIdxs=[];
                [XOut,FOut,Chi2]=lsqlinAutoFitFluorescenceCurve(FluoIn,FitIdxs,XPrevious,Delay,TimeRes,FluoError);
            else
                FitIdxs=Combs(CombIdx,:);
                [XOut,FOut,Chi2]=lsqlinAutoFitFluorescenceCurve(FluoIn,FitIdxs,XPrevious,Delay,TimeRes,FluoError);
            end
            % Convert output X to rates
            FitRates=XOut/TimeRes;

            % Sub fit rates into all rates
            TempRates=Rates;
            TempRates(MinFitIdx:MaxFitIdx)=FitRates;

            %%% Calculate GOF
            % Number of parameters = 2*NCombinations (1 for index, 1 for
            % rate)
            k = 2*NFree;
            % Akaike information criterion
            AIC = Chi2+2*k;
            % Corrected AIC for finite sample size
            AICc = AIC + 2*k*(k+1)/(N-k-1);
            GOF = AICc;
            if GOF < BestGOF || BestGOF == -1
                BestGOF=GOF;
                BestRates=TempRates;
                % Plot fit in real time
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
    % Update Rates from BestRates
    Rates=BestRates;
end

% To simplify calculations, we had Rate(i) have its first finite value at
% Fluo(i). Should be Fluo(i+1), so subtract TimeRes from Time
TransitionTimes=TimeInterp-TimeRes;

end

