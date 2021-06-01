function [MeanValues, StdErrors, TimeVector, Ymax, Ymin, YLabel, OutString, LogYmax, LogYmin] = ...
    GetPlottingMats(CompiledParameters, hmmVarString, UseRescaledTiming, UseRescaledFluo,...
    UseRescaledParamTiming, nrows, ReferenceTemperature)

if ~exist('ReferenceTemperature', 'var')
    ReferenceTemperature = 25;
end

if single(round(ReferenceTemperature, 0)) == single(ReferenceTemperature)
    TempString = sprintf('%.0f ', ReferenceTemperature);
else
    TempString = sprintf('%.1f ', ReferenceTemperature);
end

shortYlabel = false;



if nrows > 2
    shortYlabel = true;
end

if ~UseRescaledTiming & ~UseRescaledParamTiming
    TimeVector = CompiledParameters.TimeVector;
else
    TimeVector = CompiledParameters.ScaledTimeVector;
end


if strcmpi(hmmVarString, 'duration') | strcmpi(hmmVarString, 'durations') | ...
        strcmpi(hmmVarString, 'burstduration') | strcmpi(hmmVarString, 'burstdurations')
    if ~UseRescaledTiming & ~UseRescaledParamTiming
        MeanValues = CompiledParameters.Durations;
        StdErrors = CompiledParameters.DurationsStdErr;
    elseif UseRescaledParamTiming 
        MeanValues = CompiledParameters.ScaledDurationsV2;
        StdErrors = CompiledParameters.ScaledDurationsStdErrV2;
    else
        MeanValues = CompiledParameters.ScaledDurations;
        StdErrors = CompiledParameters.ScaledDurationsStdErr;
    end
    
    
        
    YminA = nanmin(nanmin(nanmin(MeanValues)));
    YmaxA = nanmax(nanmax(nanmax(MeanValues)));
    YminB = nanmin(nanmin(nanmin(MeanValues-StdErrors)));
    YmaxB = nanmax(nanmax(nanmax(MeanValues+StdErrors)));
    Ymax = min(YmaxA*1.2, YmaxB);
    Ymin = max(YminA*5/6, YminB);
    Ymin = max(floor(Ymin/2)*2, 0);
    Ymax = min(ceil(Ymax/2)*2, 10);
    LogYmin = Ymin;
    LogYmax = Ymax;
    
    if ~shortYlabel
        if ~UseRescaledParamTiming
            YLabel = 'burst duration (min)';
        else
            
            YLabel = ['burst duration (',TempString,'ºC min)'];
            
        end
    else
        if ~UseRescaledParamTiming
            YLabel = 'burst dur. (m)';
        else
            YLabel = ['burst dur. (',TempString,'ºC m)'];
        end
    end
    
    OutString = 'burst_durations';
    
elseif strcmpi(hmmVarString, 'frequency') | strcmpi(hmmVarString, 'frequencies') | ...
        strcmpi(hmmVarString, 'burstfrequency') | strcmpi(hmmVarString, 'burstfrequencies')
    if ~UseRescaledTiming & ~UseRescaledParamTiming
        MeanValues = CompiledParameters.Frequencies;
        StdErrors = CompiledParameters.FrequenciesStdErr;
    elseif UseRescaledParamTiming 
        MeanValues = CompiledParameters.ScaledFrequenciesV2;
        StdErrors = CompiledParameters.ScaledFrequenciesStdErrV2;
    else
        MeanValues = CompiledParameters.ScaledFrequencies;
        StdErrors = CompiledParameters.ScaledFrequenciesStdErr;
    end
    
    YminA = nanmin(nanmin(nanmin(MeanValues)));
    YmaxA = nanmax(nanmax(nanmax(MeanValues)));
    YminB = nanmin(nanmin(nanmin(MeanValues-StdErrors)));
    YmaxB = nanmax(nanmax(nanmax(MeanValues+StdErrors)));
    Ymax = min(YmaxA*1.2, YmaxB);
    Ymin = max(YminA*5/6, YminB);
    
    Ymin = max(floor(Ymin/0.1)*0.1, 0);
    Ymax = min(ceil(Ymax/0.1)*0.1, 2);
    LogYmin = Ymin;
    LogYmax = Ymax;
    
    if ~shortYlabel
        if ~UseRescaledParamTiming
            YLabel = 'burst frequency (1/min)';
        else
            YLabel = ['burst frequency (1/',TempString,'ºC min)'];
        end
    else
        if ~UseRescaledParamTiming
            YLabel = 'burst freq. (1/m)';
        else
            YLabel = ['burst freq. (1/',TempString,'ºC m)'];
        end
    end
    
    OutString = 'burst_frequencies';
    
elseif strcmpi(hmmVarString, 'initiationrate') | strcmpi(hmmVarString, 'initiationrates') | ...
        strcmpi(hmmVarString, 'burstinitiationrate') | strcmpi(hmmVarString, 'burstinitiationrates')
    
    if ~UseRescaledTiming & ~UseRescaledParamTiming
        MeanValues = CompiledParameters.InitiationRates;
        StdErrors = CompiledParameters.InitiationRatesStdErr;
    elseif UseRescaledParamTiming 
        MeanValues = CompiledParameters.ScaledInitiationRatesV2;
        StdErrors = CompiledParameters.ScaledInitiationRatesStdErrV2;
    else
        MeanValues = CompiledParameters.ScaledInitiationRates;
        StdErrors = CompiledParameters.ScaledInitiationRatesStdErr;
    end
    
    if UseRescaledFluo 
        for i = 1:size(MeanValues, 1)
            match_index = find(round(CompiledParameters.UniqueTemperatures, 2) == round(CompiledParameters.SetTemperatures(i), 2));
            if ~isempty(match_index)
                MeanValues(i, :,:) = CompiledParameters.FluoCoeffs(match_index)*MeanValues(i, :,:);
                StdErrors(i, :,:) = CompiledParameters.FluoCoeffs(match_index)*StdErrors(i, :,:);
            end
        end
    end
    
        
    
    YminA = nanmin(nanmin(nanmin(MeanValues)));
    YmaxA = nanmax(nanmax(nanmax(MeanValues)));
    YminB = nanmin(nanmin(nanmin(MeanValues-StdErrors)));
    YmaxB = nanmax(nanmax(nanmax(MeanValues+StdErrors)));
    Ymax = min(YmaxA*1.2, YmaxB);
    Ymin = max(YminA*5/6, YminB);
    
    Ymin = max(floor(Ymin/50)*50, 0);
    Ymax = ceil(Ymax/50)*50;
    LogYmin = Ymin;
    LogYmax = 10^(ceil(log10(Ymax)));
    
    if ~shortYlabel
        if ~UseRescaledFluo & ~UseRescaledParamTiming
            YLabel = 'initiation rate (au/min)';
        elseif ~UseRescaledFluo
            YLabel = ['initiation rate (au/',TempString,'ºC min)'];
        elseif ~UseRescaledParamTiming
            YLabel = 'scaled init. rate (au/min)';
        else
            YLabel = ['scaled init. rate (au/',TempString,'ºC min)'];
        end
    else
        if ~UseRescaledFluo & ~UseRescaledParamTiming
            YLabel = 'init. rate (au/m)';
        elseif ~UseRescaledParamTiming
            YLabel = 'scaled init. rate (au/m)';
        elseif ~UseRescaledFluo
            YLabel = ['init. rate (au/',TempString,'ºC  m)'];
        else
             YLabel = ['scaled init. rate (au/',TempString,'ºC m)'];
        end
    end
    
    OutString = 'burst_initiation_rates';
elseif strcmpi(hmmVarString, 'cycletime') | strcmpi(hmmVarString, 'cycletimes') | ...
        strcmpi(hmmVarString, 'burstcycletime') | strcmpi(hmmVarString, 'burstcycletimes')
    if ~UseRescaledTiming & ~UseRescaledParamTiming
        koff = 1./CompiledParameters.Durations;
        koffSE = CompiledParameters.DurationsStdErr./(CompiledParameters.Durations.^2);
        kon = CompiledParameters.Frequencies;
        konSE = CompiledParameters.FrequenciesStdErr;
    elseif UseRescaledParamTiming 
        koff = 1./CompiledParameters.ScaledDurationsV2;
        koffSE = CompiledParameters.ScaledDurationsStdErrV2./(CompiledParameters.ScaledDurationsV2.^2);
        kon = CompiledParameters.ScaledFrequenciesV2;
        konSE = CompiledParameters.ScaledFrequenciesStdErrV2;
    else
        koff = 1./CompiledParameters.ScaledDurations;
        koffSE = CompiledParameters.ScaledDurationsStdErr./(CompiledParameters.ScaledDurations.^2);
        kon = CompiledParameters.ScaledFrequencies;
        konSE = CompiledParameters.ScaledFrequenciesStdErr;
    end
    
    
    MeanValues = (kon+koff)./(kon.*koff);
    StdErrors = sqrt((konSE.^2)./(kon.^4) + (koffSE.^2)./(koff.^4));
    
        
    YminA = nanmin(nanmin(nanmin(MeanValues)));
    YmaxA = nanmax(nanmax(nanmax(MeanValues)));
    YminB = nanmin(nanmin(nanmin(MeanValues-StdErrors)));
    YmaxB = nanmax(nanmax(nanmax(MeanValues+StdErrors)));
    Ymax = min(YmaxA*1.2, YmaxB);
    Ymin = max(YminA*5/6, YminB);
    Ymin = max(floor(Ymin/2)*2, 0);
    Ymax = min(ceil(Ymax/2)*2, 10);
    LogYmin = Ymin;
    LogYmax = Ymax;
    
    if ~shortYlabel
        if ~UseRescaledParamTiming
            YLabel = 'burst cycle time (min)';
        else
            YLabel = ['burst cycle time (',TempString,'ºC min)'];
        end
    else
        if ~UseRescaledParamTiming
            YLabel = 'cycle time (m)';
        else
            YLabel = ['cycle time (',TempString,'ºC m)'];
        end
    end
    
    OutString = 'burst_cycle_times';
elseif strcmpi(hmmVarString, 'meaninitiationrate') | strcmpi(hmmVarString, 'meaninitiationrates') | ...
        strcmpi(hmmVarString, 'meanburstinitiationrate') | strcmpi(hmmVarString, 'meanburstinitiationrates')
    if ~UseRescaledTiming & ~UseRescaledParamTiming
        InitRates = CompiledParameters.InitiationRates;
        InitRateSEs = CompiledParameters.InitiationRatesStdErr;
        koff = 1./CompiledParameters.Durations;
        koffSE = CompiledParameters.DurationsStdErr./(CompiledParameters.Durations.^2);
        kon = CompiledParameters.Frequencies;
        konSE = CompiledParameters.FrequenciesStdErr;
    elseif UseRescaledParamTiming 
        InitRates = CompiledParameters.ScaledInitiationRatesV2;
        InitRateSEs = CompiledParameters.ScaledInitiationRatesStdErrV2;
        koff = 1./CompiledParameters.ScaledDurationsV2;
        koffSE = CompiledParameters.ScaledDurationsStdErrV2./(CompiledParameters.ScaledDurationsV2.^2);
        kon = CompiledParameters.ScaledFrequenciesV2;
        konSE = CompiledParameters.ScaledFrequenciesStdErrV2;
    else
        InitRates = CompiledParameters.ScaledInitiationRates;
        InitRateSEs = CompiledParameters.ScaledInitiationRatesStdErr;
        koff = 1./CompiledParameters.ScaledDurations;
        koffSE = CompiledParameters.ScaledDurationsStdErr./(CompiledParameters.ScaledDurations.^2);
        kon = CompiledParameters.ScaledFrequencies;
        konSE = CompiledParameters.ScaledFrequenciesStdErr;
    end
    
    if UseRescaledFluo 
        for i = 1:size(InitRates, 1)
            match_index = find(round(CompiledParameters.UniqueTemperatures, 2) == round(CompiledParameters.SetTemperatures(i), 2));
            if ~isempty(match_index)
                InitRates(i, :,:) = CompiledParameters.FluoCoeffs(match_index)*InitRates(i, :,:);
                InitRateSEs(i, :,:) = CompiledParameters.FluoCoeffs(match_index)*InitRateSEs(i, :,:);
            end
        end
    end 
    


    
    MeanValues = InitRates.*(kon)./(kon+koff);
    StdErrors = sqrt(((kon./(kon+koff)).^2).*(InitRateSEs.^2) + ...
        ((InitRates.*koff./((kon+koff).^2)).^2).*(konSE.^2) + ...
        ((-InitRates.*kon./((kon+koff).^2)).^2).*(koffSE.^2)  );
    
        
    YminA = nanmin(nanmin(nanmin(MeanValues)));
    YmaxA = nanmax(nanmax(nanmax(MeanValues)));
    YminB = nanmin(nanmin(nanmin(MeanValues-StdErrors)));
    YmaxB = nanmax(nanmax(nanmax(MeanValues+StdErrors)));
    Ymax = min(YmaxA*1.2, YmaxB);
    Ymin = max(YminA*5/6, YminB);
    
    Ymin = max(floor(Ymin/50)*50, 0);
    Ymax = ceil(Ymax/50)*50;
    LogYmin = Ymin;
    LogYmax = 10^(ceil(log10(Ymax)));
    
    if ~shortYlabel
        if ~UseRescaledFluo & ~UseRescaledParamTiming
            YLabel = 'mean loading rate (au/min)';
        elseif ~UseRescaledFluo
            YLabel = ['mean loading rate (au/',TempString,'ºC min)'];
        elseif ~UseRescaledParamTiming
            YLabel = 'scaled mean loading rate (au/min)';
        else
            YLabel = ['scaled mean loading rate (au/',TempString,'ºC min)'];
        end
    else
        if ~UseRescaledFluo & ~UseRescaledParamTiming
            YLabel = 'mean loading rate (au/m)';
        elseif ~UseRescaledParamTiming
            YLabel = 'scaled mean loading rate (au/m)';
        elseif ~UseRescaledFluo
            YLabel = ['mean loading rate (au/',TempString,'ºC  m)'];
        else
             YLabel = ['scaled mean loading rate (au/',TempString,'ºC m)'];
        end
    end
    
    OutString = 'mean_initiation_rates';
else
    error('hmmVarString does not match any of the hmm parameters stored in CompiledParameters.')
end
