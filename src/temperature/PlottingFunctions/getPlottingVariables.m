function [PlottedParams, PlottedParamSEs,R2s, ylab,OutputString,GlobalPlotYmax,GlobalPlotYmin,LogPlotYmin] = ...
    getPlottingVariables(this, parameter,  TraceType, R2bound, UseRescaledFluo, UseRescaledParamTiming,...
    UsePerNucleusTraces, UseBinnedTraces, UseBinnedPerNucleusTraces)
if ~exist('UseRescaledFluo', 'var')
    UseRescaledFluo = false;
end

if ~exist('UseRescaledParamTiming', 'var')
    UseRescaledParamTiming = false;
end

if ~exist('UsePerNucleusTraces', 'var')
    UsePerNucleusTraces = false;
end
if ~exist('UseBinnedTraces', 'var')
    UseBinnedTraces = false;
end
if ~exist('UseBinnedPerNucleusTraces', 'var')
    UseBinnedPerNucleusTraces = false;
end
%% Load relevant parameters into memory
if UseBinnedPerNucleusTraces
    InitiationRates = getBinnedPerNucleusTrapezoidParameters(this, 'MeanInitiationRates', TraceType, false);
    SEInitiationRates = getBinnedPerNucleusTrapezoidParameters(this, 'MeanInitiationRates', TraceType, true);
    TimeOns = getBinnedPerNucleusTrapezoidParameters(this, 'TimeOns', TraceType, false);
    SETimeOns = getBinnedPerNucleusTrapezoidParameters(this, 'TimeOns', TraceType, true);
    TimeOffs = getBinnedPerNucleusTrapezoidParameters(this, 'TimeOffs', TraceType, false);
    SETimeOffs = getBinnedPerNucleusTrapezoidParameters(this, 'TimeOffs', TraceType, true);
    ElongationTimes = getBinnedPerNucleusTrapezoidParameters(this, 'ElongationTimes', TraceType, false);
    SEElongationTimes = getBinnedPerNucleusTrapezoidParameters(this, 'ElongationTimes', TraceType, true);
    UnloadingRates = getBinnedPerNucleusTrapezoidParameters(this, 'UnloadingRates', TraceType, false);
    SEUnloadingRates = getBinnedPerNucleusTrapezoidParameters(this, 'UnloadingRates', TraceType, true);
    R2s = getBinnedPerNucleusTrapezoidParameters(this, 'R2s', TraceType);
    [FractionOns, SchnitzCounts] = getCycleFractionOns(this, TraceType, UseBinnedTraces, UseBinnedPerNucleusTraces);
    [MaxFluos, SEMaxFluos] = getBinnedPerNucleusMaxFluoMatForPlotting(this, TraceType);
elseif UseBinnedTraces
    InitiationRates = getBinnedTrapezoidParameters(this, 'MeanInitiationRates', TraceType, false);
    SEInitiationRates = getBinnedTrapezoidParameters(this, 'MeanInitiationRates', TraceType, true);
    TimeOns = getBinnedTrapezoidParameters(this, 'TimeOns', TraceType, false);
    SETimeOns = getBinnedTrapezoidParameters(this, 'TimeOns', TraceType, true);
    TimeOffs = getBinnedTrapezoidParameters(this, 'TimeOffs', TraceType, false);
    SETimeOffs = getBinnedTrapezoidParameters(this, 'TimeOffs', TraceType, true);
    ElongationTimes = getBinnedTrapezoidParameters(this, 'ElongationTimes', TraceType, false);
    SEElongationTimes = getBinnedTrapezoidParameters(this, 'ElongationTimes', TraceType, true);
    UnloadingRates = getBinnedTrapezoidParameters(this, 'UnloadingRates', TraceType, false);
    SEUnloadingRates = getBinnedTrapezoidParameters(this, 'UnloadingRates', TraceType, true);
    R2s = getBinnedTrapezoidParameters(this, 'R2s', TraceType);
    [FractionOns, SchnitzCounts] = getCycleFractionOns(this, TraceType, UseBinnedTraces, UseBinnedPerNucleusTraces);
    [MaxFluos, SEMaxFluos] = getBinnedMaxFluoMatForPlotting(this, TraceType);
elseif UsePerNucleusTraces
        InitiationRates = getPerNucleusTrapezoidParameters(this, 'MeanInitiationRates', TraceType, false);
    SEInitiationRates = getPerNucleusTrapezoidParameters(this, 'MeanInitiationRates', TraceType, true);
    TimeOns = getPerNucleusTrapezoidParameters(this, 'TimeOns', TraceType, false);
    SETimeOns = getPerNucleusTrapezoidParameters(this, 'TimeOns', TraceType, true);
    TimeOffs = getPerNucleusTrapezoidParameters(this, 'TimeOffs', TraceType, false);
    SETimeOffs = getPerNucleusTrapezoidParameters(this, 'TimeOffs', TraceType, true);
    ElongationTimes = getPerNucleusTrapezoidParameters(this, 'ElongationTimes', TraceType, false);
    SEElongationTimes = getPerNucleusTrapezoidParameters(this, 'ElongationTimes', TraceType, true);
    UnloadingRates = getPerNucleusTrapezoidParameters(this, 'UnloadingRates', TraceType, false);
    SEUnloadingRates = getPerNucleusTrapezoidParameters(this, 'UnloadingRates', TraceType, true);
    R2s = getPerNucleusTrapezoidParameters(this, 'R2s', TraceType);
    [FractionOns, SchnitzCounts] = getCycleFractionOns(this, TraceType, UseBinnedTraces, UseBinnedPerNucleusTraces);
    [MaxFluos, SEMaxFluos] = getPerNucleusTrapezoidParameters(this, TraceType);
else
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
[FractionOns, SchnitzCounts] = getCycleFractionOns(this, TraceType, UseBinnedTraces, UseBinnedPerNucleusTraces);
[MaxFluos, SEMaxFluos] = getMaxFluoMatForPlotting(this, TraceType);
end

%%
if ~UseBinnedTraces & ~UseBinnedPerNucleusTraces
NumSets = length(this.ExperimentPrefixes);

if UseRescaledFluo
    temp_match_indices = NaN(1, NumSets);
    SetFluoCoeffs = NaN(1, NumSets);

    for i = NumSets
        temp_match_indices(i) = find(round(this.UniqueTemperatures, 2) == round(this.Temp_sps(i), 2));
        SetFluoCoeffs(i) = this.FluoCoeffs(temp_match_indices(i));
        if ~isempty(temp_match_indices(i))
            InitiationRates(i, :,:) = SetFluoCoeffs(i)*InitiationRates(i, :,:);
            SEInitiationRates(i, :,:) = SetFluoCoeffs(i)*SEInitiationRates(i, :,:);
            MeanSpotFluos(i,:,:) =  SetFluoCoeffs(i)*MeanSpotFluos(i, :,:);
            SEMeanSpotFluos(i,:,:) =  SetFluoCoeffs(i)*SEMeanSpotFluos(i, :,:);
            MaxFluos(i,:,:) =  SetFluoCoeffs(i)*MaxFluos(i, :,:);
            SEMaxFluos(i,:,:) =  SetFluoCoeffs(i)*SEMaxFluos(i, :,:);
        else
            InitiationRates(i, :,:) = NaN(size(InitiationRates(i, :,:) ));
            SEInitiationRates(i, :,:) = NaN(size(SEInitiationRates(i, :,:) ));
            MeanSpotFluos(i, :,:) = NaN(size(MeanSpotFluos(i, :,:) ));
            SEMeanSpotFluos(i, :,:) = NaN(size(SEMeanSpotFluos(i, :,:) ));
            MaxFluos(i, :,:) = NaN(size(MaxFluos(i, :,:) ));
            SEMaxFluos(i, :,:) = NaN(size(SEMaxFluos(i, :,:) ));
        end
        
    end
end




if UseRescaledParamTiming
    temp_match_indices = NaN(1, NumSets);
    TimingCoeffs = NaN(NumSets, 6);
    for i = 1:NumSets
        temp_match_indices(i) = find(round(this.UniqueTemperatures, 2) == round(this.Temp_sps(i), 2));
        if ~isempty(temp_match_indices(i))
            for NC = 9:14
                if NC == 14
                    if ~isnan(this.TimeScalingInfo.PropNCDivisionInfo(temp_match_indices(i), NC-9))
                        TimingCoeffs(i, NC-8) =  1/this.TimeScalingInfo.PropNCDivisionInfo(temp_match_indices(i), NC-9);
                    else
                        TimingCoeffs(i, NC-8) = this.TimingCoeffs(temp_match_indices(i));
                    end
                elseif ~isnan(this.TimeScalingInfo.PropNCDivisionInfo(temp_match_indices(i), NC-8))
                     TimingCoeffs(i, NC-8) = 1/this.TimeScalingInfo.PropNCDivisionInfo(temp_match_indices(i), NC-8);
                else
                    TimingCoeffs(i, NC-8) = this.TimingCoeffs(temp_match_indices(i));
                end
                
                InitiationRates(i, :,NC-8) = InitiationRates(i, :,NC-8)/TimingCoeffs(i, NC-8);
                SEInitiationRates(i, :,NC-8) = SEInitiationRates(i, :,NC-8)/TimingCoeffs(i, NC-8);
                TimeOns(i, :,NC-8) = TimeOns(i, :,NC-8)*TimingCoeffs(i, NC-8);
                SETimeOns(i, :,NC-8) = SETimeOns(i, :,NC-8)*TimingCoeffs(i, NC-8);
                TimeOffs(i, :,NC-8) = TimeOffs(i, :,NC-8)*TimingCoeffs(i, NC-8);
                SETimeOffs(i, :,NC-8) = SETimeOffs(i, :,NC-8)*TimingCoeffs(i, NC-8);
                ElongationTimes(i, :,NC-8) = ElongationTimes(i, :,NC-8)*TimingCoeffs(i, NC-8);
                SEElongationTimes(i, :,NC-8) = SEElongationTimes(i, :,NC-8)*TimingCoeffs(i, NC-8);
                UnloadingRates(i, :,NC-8) = UnloadingRates(i, :,NC-8)/TimingCoeffs(i, NC-8);
                SEUnloadingRates(i, :,NC-8) = SEUnloadingRates(i, :,NC-8)/TimingCoeffs(i, NC-8);

            end
        else
            InitiationRates(i, :,:) = NaN(size(InitiationRates(i, :,:) ));
            SEInitiationRates(i, :,:) = NaN(size(SEInitiationRates(i, :,:) ));
            TimeOns(i, :,:) = NaN(size(TimeOns(i, :,:) ));
            SETimeOns(i, :,:) = NaN(size(SETimeOns(i, :,:) ));
            TimeOffs(i, :,:) = NaN(size(TimeOffs(i, :,:) ));
            SETimeOffs(i, :,:) = NaN(size(SETimeOffs(i, :,:) ));
            ElongationTimes(i, :,:) = NaN(size(ElongationTimes(i, :,:) ));
            SEElongationTimes(i, :,:) = NaN(size(SEElongationTimes(i, :,:) ));
            UnloadingRates(i, :,:) = NaN(size(UnloadingRates(i, :,:) ));
            SEUnloadingRates(i, :,:) = NaN(size(SEUnloadingRates(i, :,:) ));
        end
    end
    
end

else
    NumTemperatures = length(this.UniqueTemperatures);

if UseRescaledFluo

    SetFluoCoeffs = this.FluoCoeffs;

    for i = NumTemperatures
            InitiationRates(i, :,:) = SetFluoCoeffs(i)*InitiationRates(i, :,:);
            SEInitiationRates(i, :,:) = SetFluoCoeffs(i)*SEInitiationRates(i, :,:);
            MeanSpotFluos(i,:,:) =  SetFluoCoeffs(i)*MeanSpotFluos(i, :,:);
            SEMeanSpotFluos(i,:,:) =  SetFluoCoeffs(i)*SEMeanSpotFluos(i, :,:);
            MaxFluos(i,:,:) =  SetFluoCoeffs(i)*MaxFluos(i, :,:);
            SEMaxFluos(i,:,:) =  SetFluoCoeffs(i)*SEMaxFluos(i, :,:);

        
    end
end




if UseRescaledParamTiming
    TimingCoeffs = NaN(NumTemperatures, 6);
    for i = 1:NumTemperatures
        for NC = 9:14
            if NC == 14
                if ~isnan(this.TimeScalingInfo.PropNCDivisionInfo(i, NC-9))
                    TimingCoeffs(i, NC-8) =  1/this.TimeScalingInfo.PropNCDivisionInfo(i, NC-9);
                else
                    TimingCoeffs(i, NC-8) = this.TimingCoeffs(i);
                end
            elseif ~isnan(this.TimeScalingInfo.PropNCDivisionInfo(i, NC-8))
                TimingCoeffs(i, NC-8) = 1/this.TimeScalingInfo.PropNCDivisionInfo(i, NC-8);
            else
                TimingCoeffs(i, NC-8) = this.TimingCoeffs(i);
            end
            
            InitiationRates(i, :,NC-8) = InitiationRates(i, :,NC-8)/TimingCoeffs(i, NC-8);
            SEInitiationRates(i, :,NC-8) = SEInitiationRates(i, :,NC-8)/TimingCoeffs(i, NC-8);
            TimeOns(i, :,NC-8) = TimeOns(i, :,NC-8)*TimingCoeffs(i, NC-8);
            SETimeOns(i, :,NC-8) = SETimeOns(i, :,NC-8)*TimingCoeffs(i, NC-8);
            TimeOffs(i, :,NC-8) = TimeOffs(i, :,NC-8)*TimingCoeffs(i, NC-8);
            SETimeOffs(i, :,NC-8) = SETimeOffs(i, :,NC-8)*TimingCoeffs(i, NC-8);
            ElongationTimes(i, :,NC-8) = ElongationTimes(i, :,NC-8)*TimingCoeffs(i, NC-8);
            SEElongationTimes(i, :,NC-8) = SEElongationTimes(i, :,NC-8)*TimingCoeffs(i, NC-8);
            UnloadingRates(i, :,NC-8) = UnloadingRates(i, :,NC-8)/TimingCoeffs(i, NC-8);
            SEUnloadingRates(i, :,NC-8) = SEUnloadingRates(i, :,NC-8)/TimingCoeffs(i, NC-8);
            
        end
        
    end
    
end
end


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
    if UseBinnedTraces | UseBinnedPerNucleusTraces | UsePerNucleusTraces
        error('Cannot plot mean spot values for PerNucleus or Binned Traces because calculation was not done.')
    end
    OutputString = 'MeanSpotFluos';
    PlottedParams = MeanSpotFluos;
    PlottedParamSEs = SEMeanSpotFluos;
    R2s = ones(size(PlottedParams));
    R2s(isnan(PlottedParams)) = NaN;
    R2s(SchnitzCounts < this.MinimumSchnitzCount)  = 0;
    if UseRescaledFluo
        ylab = 're-scaled mean spot fluo. (AU)';
    else
        ylab = 'mean spot fluorescence (AU)';
    end
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = GlobalPlotYmax*1.05;
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(10, GlobalPlotYmin);
elseif strcmpi(parameter, 'timeons') | strcmpi(parameter, 'timeon')
    OutputString = 'TimeOn';
    PlottedParams = TimeOns;
    PlottedParamSEs = SETimeOns;
    if UseRescaledParamTiming
        ylab = 'T-adjusted t_{on} (25ºC min)';
    else
        ylab = 't_{on} (min)';
    end
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = min(20, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(.1, GlobalPlotYmin);
    
elseif strcmpi(parameter, 'timeoffs') | strcmpi(parameter, 'timeoff')
    OutputString = 'TimeOff';
    PlottedParams = TimeOffs;
    PlottedParamSEs = SETimeOffs;
    if UseRescaledParamTiming
        ylab = 'T-adjusted t_{off} (25ºC min)';
    else
        ylab = 't_{off} (min)';
    end
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = min(60, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(.1, GlobalPlotYmin);
elseif strcmpi(parameter, 'elongationtimes') | strcmpi(parameter, 'elongationtime')
    OutputString = 'ElongationTime';
    PlottedParams = ElongationTimes;
    PlottedParamSEs = SEElongationTimes;
    if UseRescaledParamTiming
        ylab = 'T-adjusted elongation time (25ºC min)';
    else
        ylab = 'Elongation time (min)';
    end
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = min(30, GlobalPlotYmax);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(1, GlobalPlotYmin);
elseif strcmpi(parameter, 'elongationrates') | strcmpi(parameter, 'elongationrate')
    OutputString = 'ElongationRate';
    PlottedParams = this.GeneLength/ElongationTimes;
    PlottedParamSEs = PlottedParams.*SEElongationTimes./ElongationTimes;
    if UseRescaledParamTiming
        ylab = 'T-adjusted elongation rate (bp/25ºC min)';
    else
        ylab = 'Elongation rate (bp/min)';
    end
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
    if UseRescaledFluo

        if UseRescaledParamTiming
            ylab = 'T-adjusted re-scaled loading rate (AU/25ºC min)';
        else
            ylab = 're-scaled loading rate (AU/min)';
        end
    else
        if UseRescaledParamTiming
            ylab = 'T-adjusted loading rate (AU/25ºC min)';
        else
            ylab = 'loading rate (AU/min)';
        end
    end
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = min(5000, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(10, GlobalPlotYmin);
elseif strcmpi(parameter, 'transcriptionwindows') | strcmpi(parameter, 'transcriptionwindow')
    OutputString = 'TranscriptionWindow';
    PlottedParams = TimeOffs-TimeOns;
    PlottedParamSEs = sqrt(SETimeOffs.^2+SETimeOns.^2);
    if UseRescaledParamTiming
        ylab = 'T-adjusted txn window (25ºC min)';
    else
        ylab = 'Transcriotion Window Duration (min)';
    end
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = min(60, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(.1, GlobalPlotYmin);
elseif strcmpi(parameter, 'maxfluos') | strcmpi(parameter, 'maxfluo')
    OutputString = 'MaxFluos';
    PlottedParams = MaxFluos;
    PlottedParamSEs = SEMaxFluos;
    if UseRescaledFluo
        ylab = 're-scaled max fluoresence (AU)';
    else
        ylab = 'max fluoresence (AU)';
    end
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
    if UseRescaledFluo
        ylab = 're-scaled Fitted max polymerases loaded (AU)';
    else
        ylab = 'Fitted max polymerases loaded (AU)';
    end
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    
    GlobalPlotYmax = min(12000, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(10, GlobalPlotYmin);
elseif strcmpi(parameter, 'unloadingrates') | strcmpi(parameter, 'unloadingrate')
    OutputString = 'UnloadingRate';
    PlottedParams = -UnloadingRates;
    PlottedParamSEs = SEUnloadingRates;
    if UseRescaledFluo

        if UseRescaledParamTiming
            ylab = '- T-adjusted re-scaled unloading rate (AU/25ºC min)';
        else
            ylab = '- re-scaled unloading Rate (AU/min)';
        end
    else
        if UseRescaledParamTiming
            ylab = '- T-adjusted unloading rate (AU/25ºC min)';
        else
            ylab = '- unloading Rate (AU/min)';
        end
    end
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    GlobalPlotYmax = min(5000, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(10, GlobalPlotYmin);
elseif strcmpi(parameter, 'mrnaproductions') | strcmpi(parameter, 'mrnaproduction') | strcmpi(parameter, 'totalmrnaproductions') | strcmpi(parameter, 'totalmrnaproduction')
    OutputString = 'TotalmRNAProduction';
    if ~UseBinnedPerNucleusTraces & ~UsePerNucleusTraces
        FractionOnsTemp = FractionOns;
        FractionOnsTemp(SchnitzCounts < this.MinimumSchnitzCount) = NaN;
        PlottedParams = FractionOnsTemp.*InitiationRates.*(TimeOffs-TimeOns);
        PlottedParamSEs = sqrt(((FractionOnsTemp.^2).*(TimeOffs-TimeOns).^2).*(SEInitiationRates.^2)+...
            (FractionOnsTemp.^2).*(InitiationRates.^2).*(SETimeOffs.^2+SETimeOns.^2));
    else
        PlottedParams = InitiationRates.*(TimeOffs-TimeOns);
        PlottedParamSEs = sqrt(((TimeOffs-TimeOns).^2).*(SEInitiationRates.^2)+...
            (InitiationRates.^2).*(SETimeOffs.^2+SETimeOns.^2));
    end
    
    if UseRescaledFluo
        ylab = 're-scaled mean mRNA produced per cell (AU)';
    else
        ylab = 'mean mRNA produced per cell (AU)';
    end
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    
    GlobalPlotYmax = min(100000, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(10, GlobalPlotYmin);
elseif strcmpi(parameter, 'unweightedmrnaproductions') | strcmpi(parameter, 'unweightedmrnaproduction') | strcmpi(parameter, 'unweightedtotalmrnaproductions') | strcmpi(parameter, 'unweightedtotalmrnaproduction')
    OutputString = 'UnweightedTota lmRNAProduction';
    
    PlottedParams = InitiationRates.*(TimeOffs-TimeOns);
    PlottedParamSEs = sqrt(((TimeOffs-TimeOns).^2).*(SEInitiationRates.^2)+...
        (InitiationRates.^2).*(SETimeOffs.^2+SETimeOns.^2));
    if UseRescaledFluo
        ylab = 're-scaled mean mRNA produced per cell (AU)';
    else
        ylab = 'mean mRNA produced per cell (AU)';
    end
     [~, ~, GlobalPlotYmax, GlobalPlotYmin] =...
        getGlobalYlims(PlottedParams, PlottedParamSEs, R2s, R2bound);
    
    GlobalPlotYmax = min(100000, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(10, GlobalPlotYmin);
else
    error('Invalid choice of parameter.')
end

