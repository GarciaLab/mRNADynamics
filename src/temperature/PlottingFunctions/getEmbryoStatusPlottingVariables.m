function [PlottedParams, PlottedParamSEs,ylab,OutputString,GlobalPlotYmax,GlobalPlotYmin,LogPlotYmin] = ...
    getEmbryoStatusPlottingVariables(this, parameter)
%% Load relevant parameters into memory
SchnitzCount = getEmbryoStatsParameters(this, 'SchnitzCount');
FractionSickNuclei = getEmbryoStatsParameters(this, 'FractionSickNuclei');
FractionRejectedNuclei = getEmbryoStatsParameters(this, 'FractionRejectedNuclei');
FractionCompleteNuclei = getEmbryoStatsParameters(this, 'FractionCompleteNuclei');
FractionFirstLastFrameNuclei = getEmbryoStatsParameters(this, 'FractionFirstLastFrameNuclei');
[TotalXDistanceTraveled, TotalXDistanceTraveledStd] = getEmbryoStatsParameters(this, 'MeanTotalXDistanceTraveled');
[TotalYDistanceTraveled, TotalYDistanceTraveledStd] = getEmbryoStatsParameters(this, 'MeanTotalYDistanceTraveled');
[TotalDistanceTraveled, TotalDistanceTraveledStd] = getEmbryoStatsParameters(this, 'MeanTotalDistanceTraveled');
[DistanceTraveledPerSecond, DistanceTraveledPerSecondStd] = getEmbryoStatsParameters(this, 'MeanDistanceTraveledPerSecond');
[TotalXDisplacement, TotalXDisplacementStd] = getEmbryoStatsParameters(this, 'MeanTotalXDisplacement');
[TotalYDisplacement, TotalYDisplacementStd] = getEmbryoStatsParameters(this, 'MeanTotalYDisplacement');
[TotalDisplacement, TotalDisplacementStd] = getEmbryoStatsParameters(this, 'MeanTotalDisplacement');
[DisplacementPerSecond, DisplacementPerSecondStd] = getEmbryoStatsParameters(this, 'MeanDisplacementPerSecond');
[NCDivisionTimes, NCDivisionTimeStd] = getEmbryoStatsParameters(this, 'NCDivisionInfo');

%%
if strcmpi(parameter, 'FractionSickNuclei')
    OutputString = 'FractionSickNuclei';
    PlottedParams = FractionSickNuclei;
    PlottedParamSEs = NaN(size(PlottedParams));
    ylab = 'sick nuclear fraction';
    GlobalPlotYmax = 1.05;
    GlobalPlotYmin = -0.05;
    LogPlotYmin = 0.05;
elseif strcmpi(parameter, 'FractionHealthyNuclei')
    OutputString = 'FractionHealthyNuclei';
    PlottedParams = 1-FractionSickNuclei;
    PlottedParamSEs = NaN(size(PlottedParams));
    ylab = 'healthy nuclear fraction';
    GlobalPlotYmax = 1.05;
    GlobalPlotYmin = -0.05;
    LogPlotYmin = 0.05;
elseif strcmpi(parameter, 'FractionRejectedNuclei')
    OutputString = 'FractionRejectedNuclei';
    PlottedParams = FractionRejectedNuclei;
    PlottedParamSEs = NaN(size(PlottedParams));
    ylab = 'rejected nuclear fraction';
    GlobalPlotYmax = 1.05;
    GlobalPlotYmin = -0.05;
    LogPlotYmin = 0.05;
elseif strcmpi(parameter, 'FractionApprovedNuclei')
    OutputString = 'FractionApprovedNuclei';
    PlottedParams = 1-FractionRejectedNuclei;
    PlottedParamSEs = NaN(size(PlottedParams));
    ylab = 'approved nuclear fraction';
    GlobalPlotYmax = 1.05;
    GlobalPlotYmin = -0.05;
    LogPlotYmin = 0.05;
elseif strcmpi(parameter, 'FractionCompleteNuclei')
    OutputString = 'FractionCompleteNuclei';
    PlottedParams = FractionCompleteNuclei;
    PlottedParamSEs = NaN(size(PlottedParams));
    ylab = 'complete cycle nuclear fraction';
    GlobalPlotYmax = 1.05;
    GlobalPlotYmin = -0.05;
    LogPlotYmin = 0.05;
elseif strcmpi(parameter, 'FractionFirstLastFrameNuclei')
    OutputString = 'FractionFirstLastFrameNuclei';
    PlottedParams = FractionFirstLastFrameNuclei;
    PlottedParamSEs = NaN(size(PlottedParams));
    ylab = 'first/last frame nuclear fraction';
    GlobalPlotYmax = 1.05;
    GlobalPlotYmin = -0.05;
    LogPlotYmin = 0.05;
elseif strcmpi(parameter, 'TotalXDistanceTraveled')
    OutputString = 'TotalXDistanceTraveled';
    PlottedParams = TotalXDistanceTraveled;
    PlottedParamSEs = TotalXDistanceTraveledStd;
    ylab = '';
    GlobalPlotYmax = 1.05;
    GlobalPlotYmin = -0.05;
    LogPlotYmin = 0.05;
  
    
    
    
    
    
    
elseif strcmpi(parameter, 'timeons') | strcmpi(parameter, 'timeon')
    OutputString = 'TimeOn';
    PlottedParams = TimeOns;
    PlottedParamSEs = SETimeOns;
    ylab = 't_{on} (min)';
    [GlobalPlotYmax, GlobalPlotYmin] = getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = min(20, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(.1, GlobalPlotYmin);
    
elseif strcmpi(parameter, 'timeoffs') | strcmpi(parameter, 'timeoff')
    OutputString = 'TimeOff';
    PlottedParams = TimeOffs;
    PlottedParamSEs = SETimeOffs;
    ylab = 't_{off} (min)';
    [GlobalPlotYmax, GlobalPlotYmin] = getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = min(60, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(.1, GlobalPlotYmin);
elseif strcmpi(parameter, 'elongationtimes') | strcmpi(parameter, 'elongationtime')
    OutputString = 'ElongationTime';
    PlottedParams = ElongationTimes;
    PlottedParamSEs = SEElongationTimes;
    ylab = 'Elongation time (min)';
    [GlobalPlotYmax, GlobalPlotYmin] = getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = min(30, GlobalPlotYmax);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(1, GlobalPlotYmin);
elseif strcmpi(parameter, 'elongationrates') | strcmpi(parameter, 'elongationrate')
    OutputString = 'ElongationRate';
    PlottedParams = this.GeneLength/ElongationTimes;
    PlottedParamSEs = PlottedParams.*SEElongationTimes./ElongationTimes;
    ylab = 'Elongation rate (bp/min)';
    [GlobalPlotYmax, GlobalPlotYmin] = getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = min(10000, GlobalPlotYmax);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(100, GlobalPlotYmin);
elseif strcmpi(parameter, 'loadingrates') | strcmpi(parameter, 'loadingrate') | ...
        strcmpi(parameter, 'initiationrates') | strcmpi(parameter, 'initiationrate') |...
        strcmpi(parameter, 'meaninitiationrates') | strcmpi(parameter, 'meaninitiationrate') 
    OutputString = 'LoadingRate';
    PlottedParams = InitiationRates;
    PlottedParamSEs = SEInitiationRates;
    ylab = 'Loading Rate (AU/min)';
    [GlobalPlotYmax, GlobalPlotYmin] = getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = min(5000, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(10, GlobalPlotYmin);
elseif strcmpi(parameter, 'transcriptionwindows') | strcmpi(parameter, 'transcriptionwindow')
    OutputString = 'TranscriptionWindow';
    PlottedParams = TimeOffs-TimeOns;
    PlottedParamSEs = sqrt(SETimeOffs.^2+SETimeOns.^2);
    ylab = 'Transcriotion Window Duration (min)';
    [GlobalPlotYmax, GlobalPlotYmin] = getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = min(60, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(.1, GlobalPlotYmin);
elseif strcmpi(parameter, 'maxfluos') | strcmpi(parameter, 'maxfluo')
    OutputString = 'MaxFluos';
    PlottedParams = MaxFluos;
    PlottedParamSEs = SEMaxFluos;
    ylab = 'Max Fluoresence (AU)';
    [GlobalPlotYmax, GlobalPlotYmin] = getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = GlobalPlotYmax*1.05;
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(10, GlobalPlotYmin);
    
elseif strcmpi(parameter, 'plateauheights') | strcmpi(parameter, 'plateauheight')
    OutputString = 'PlateauHeights';
    PlottedParams = InitiationRates.*ElongationTimes;
    PlottedParamSEs = sqrt((InitiationRates.^2).*(SEElongationTimes.^2)+...
        (ElongationTimes.^2).*(SEInitiationRates.^2) );
    ylab = 'Fitted max polymerases loaded (AU)';
    [GlobalPlotYmax, GlobalPlotYmin] = getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    
    GlobalPlotYmax = min(12000, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(10, GlobalPlotYmin);
elseif strcmpi(parameter, 'unloadingrates') | strcmpi(parameter, 'unloadingrate')
    OutputString = 'UnloadingRate';
    PlottedParams = -UnloadingRates;
    PlottedParamSEs = SEUnloadingRates;
    ylab = '- unloading Rate (AU/min)';
    [GlobalPlotYmax, GlobalPlotYmin] = getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = min(5000, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(10, GlobalPlotYmin);
elseif strcmpi(parameter, 'mrnaproductions') | strcmpi(parameter, 'mrnaproduction') | strcmpi(parameter, 'totalmrnaproductions') | strcmpi(parameter, 'totalmrnaproduction')
    OutputString = 'TotalmRNAProduction';
    FractionOnsTemp = FractionOns;
    FractionOnsTemp(SchnitzCounts < this.MinimumSchnitzCount) = NaN;
    PlottedParams = FractionOnsTemp.*InitiationRates.*(TimeOffs-TimeOns);
    PlottedParamSEs = sqrt(((FractionOnsTemp.^2).*(TimeOffs-TimeOns).^2).*(SEInitiationRates.^2)+...
        (FractionOnsTemp.^2).*(InitiationRates.^2).*(SETimeOffs.^2+SETimeOns.^2));
    ylab = 'mean mRNA produced per cell (AU)';
    [GlobalPlotYmax, GlobalPlotYmin] = getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    
    GlobalPlotYmax = min(100000, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(10, GlobalPlotYmin);
elseif strcmpi(parameter, 'unweightedmrnaproductions') | strcmpi(parameter, 'unweightedmrnaproduction') | strcmpi(parameter, 'unweightedtotalmrnaproductions') | strcmpi(parameter, 'unweightedtotalmrnaproduction')
    OutputString = 'UnweightedTotalmRNAProduction';
    
    PlottedParams = InitiationRates.*(TimeOffs-TimeOns);
    PlottedParamSEs = sqrt(((TimeOffs-TimeOns).^2).*(SEInitiationRates.^2)+...
        (InitiationRates.^2).*(SETimeOffs.^2+SETimeOns.^2));
    ylab = 'mean mRNA produced per cell (AU)';
    [GlobalPlotYmax, GlobalPlotYmin] = getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    
    GlobalPlotYmax = min(100000, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(10, GlobalPlotYmin);
else
    error('Invalid choice of parameter.')
end

