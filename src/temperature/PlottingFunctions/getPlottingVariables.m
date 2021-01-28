function [PlottedParams, PlottedParamSEs,R2s, ylab,OutputString,GlobalPlotYmax,GlobalPlotYmin,LogPlotYmin] = ...
    getPlottingVariables(this, parameter,  TraceType, R2bound)
%% Load relevant parameters into memory
InitiationRates = getTrapezoidParameters(this, 'MeanInitiationRates', TraceType, false);
SEInitiationRates = getTrapezoidParameters(this, 'MeanInitiationRates', TraceType, true);
TimeOns = getTrapezoidParameters(this, 'TimeOns', TraceType, false);
SETimeOns = getTrapezoidParameters(this, 'TimeOns', TraceType, true);
TimeOffs = getTrapezoidParameters(this, 'TimeOffs', TraceType, false);
SETimeOffs = getTrapezoidParameters(this, 'TimeOffs', TraceType, true);
ElongationTimes = getTrapezoidParameters(this, 'ElongationTimes', TraceType, false);
SEElongationTimes = getTrapezoidParameters(this, 'ElongationTimes', TraceType, true);
UnloadingRates = getTrapezoidParameters(this, 'UnloadingRates', TraceType, false);
SEUnloadingRates = getTrapezoidParameters(this, 'UnloadingRates', TraceType, true);
MeanSpotFluos = getTrapezoidParameters(this, 'MeanSpotFluos', TraceType, false);
SEMeanSpotFluos = getTrapezoidParameters(this, 'MeanSpotFluos', TraceType, true);
R2s = getTrapezoidParameters(this, 'R2s', TraceType);
FractionOns = this.FractionOns;
SchnitzCounts = this.SchnitzCounts;
[MaxFluos, SEMaxFluos] = getMaxFluoMatForPlotting(this, TraceType);
%%
if strcmpi(parameter, 'fractionons') | strcmpi(parameter, 'fractionon') 
    OutputString = 'FractionOn';
    PlottedParams = FractionOns;
    PlottedParams(SchnitzCounts < this.MinimumSchnitzCount) = NaN;
    PlottedParamSEs = NaN(size(PlottedParams));
    R2s = ones(size(PlottedParams));
    R2s(isnan(PlottedParams)) = NaN;
    R2s(SchnitzCounts < this.MinimumSchnitzCount)  = 0;
    ylab = 'Fraction Competent';
    GlobalPlotYmax = 1.05;
    GlobalPlotYmin = -0.05;
    LogPlotYmin = 0.05;
elseif strcmpi(parameter, 'meanspotfluos') | strcmpi(parameter, 'meanspotfluo') 
    OutputString = 'MeanSpotFluos';
    PlottedParams = MeanSpotFluos;
    PlottedParamSEs = SEMeanSpotFluos;
    R2s = ones(size(PlottedParams));
    R2s(isnan(PlottedParams)) = NaN;
    R2s(SchnitzCounts < this.MinimumSchnitzCount)  = 0;
    ylab = 'mean spot fluorescence (AU)';
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = GlobalPlotYmax*1.05;
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(10, GlobalPlotYmin);
elseif strcmpi(parameter, 'timeons') | strcmpi(parameter, 'timeon')
    OutputString = 'TimeOn';
    PlottedParams = TimeOns;
    PlottedParamSEs = SETimeOns;
    ylab = 't_{on} (min)';
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = min(20, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(.1, GlobalPlotYmin);
    
elseif strcmpi(parameter, 'timeoffs') | strcmpi(parameter, 'timeoff')
    OutputString = 'TimeOff';
    PlottedParams = TimeOffs;
    PlottedParamSEs = SETimeOffs;
    ylab = 't_{off} (min)';
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = min(60, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(.1, GlobalPlotYmin);
elseif strcmpi(parameter, 'elongationtimes') | strcmpi(parameter, 'elongationtime')
    OutputString = 'ElongationTime';
    PlottedParams = ElongationTimes;
    PlottedParamSEs = SEElongationTimes;
    ylab = 'Elongation time (min)';
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = min(30, GlobalPlotYmax);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(1, GlobalPlotYmin);
elseif strcmpi(parameter, 'elongationrates') | strcmpi(parameter, 'elongationrate')
    OutputString = 'ElongationRate';
    PlottedParams = this.GeneLength/ElongationTimes;
    PlottedParamSEs = PlottedParams.*SEElongationTimes./ElongationTimes;
    ylab = 'Elongation rate (bp/min)';
    [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
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
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = min(5000, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(10, GlobalPlotYmin);
elseif strcmpi(parameter, 'transcriptionwindows') | strcmpi(parameter, 'transcriptionwindow')
    OutputString = 'TranscriptionWindow';
    PlottedParams = TimeOffs-TimeOns;
    PlottedParamSEs = sqrt(SETimeOffs.^2+SETimeOns.^2);
    ylab = 'Transcriotion Window Duration (min)';
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = min(60, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(.1, GlobalPlotYmin);
elseif strcmpi(parameter, 'maxfluos') | strcmpi(parameter, 'maxfluo')
    OutputString = 'MaxFluos';
    PlottedParams = MaxFluos;
    PlottedParamSEs = SEMaxFluos;
    ylab = 'Max Fluoresence (AU)';
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = GlobalPlotYmax*1.05;
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(10, GlobalPlotYmin);
    
elseif strcmpi(parameter, 'plateauheights') | strcmpi(parameter, 'plateauheight')
    OutputString = 'PlateauHeights';
    PlottedParams = InitiationRates.*ElongationTimes;
    PlottedParamSEs = sqrt((InitiationRates.^2).*(SEElongationTimes.^2)+...
        (ElongationTimes.^2).*(SEInitiationRates.^2) );
    ylab = 'Fitted max polymerases loaded (AU)';
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    
    GlobalPlotYmax = min(12000, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(10, GlobalPlotYmin);
elseif strcmpi(parameter, 'unloadingrates') | strcmpi(parameter, 'unloadingrate')
    OutputString = 'UnloadingRate';
    PlottedParams = -UnloadingRates;
    PlottedParamSEs = SEUnloadingRates;
    ylab = '- unloading Rate (AU/min)';
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
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
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    
    GlobalPlotYmax = min(100000, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(10, GlobalPlotYmin);
elseif strcmpi(parameter, 'unweightedmrnaproductions') | strcmpi(parameter, 'unweightedmrnaproduction') | strcmpi(parameter, 'unweightedtotalmrnaproductions') | strcmpi(parameter, 'unweightedtotalmrnaproduction')
    OutputString = 'UnweightedTotalmRNAProduction';
    
    PlottedParams = InitiationRates.*(TimeOffs-TimeOns);
    PlottedParamSEs = sqrt(((TimeOffs-TimeOns).^2).*(SEInitiationRates.^2)+...
        (InitiationRates.^2).*(SETimeOffs.^2+SETimeOns.^2));
    ylab = 'mean mRNA produced per cell (AU)';
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    
    GlobalPlotYmax = min(100000, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(10, GlobalPlotYmin);
else
    error('Invalid choice of parameter.')
end

