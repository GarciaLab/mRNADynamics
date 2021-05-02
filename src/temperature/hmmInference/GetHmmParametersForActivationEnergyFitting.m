function [x, y, y_err]  = GetHmmParametersForActivationEnergyFitting(this, parameter, UseScaledFluo, UseScaledTiming, TimeIndex, APindex)
x = this.UniqueTemperatures;
R = 8.314*10^(-3); % kJ * K^(-1)*mol^(-1)

if strcmpi(parameter, 'initiationrates') | strcmpi(parameter, 'loadingrates') 
    if ~UseScaledTiming
        y =  this.InitiationRates(:, TimeIndex, APindex).';
        y_err =  this.InitiationRatesStdErr(:, TimeIndex, APindex).';
    else
        y =  this.ScaledInitiationRates(:, TimeIndex, APindex).';
        y_err =  this.ScaledInitiationRatesStdErr(:, TimeIndex, APindex).';
    end
    
    if UseScaledFluo
        for i = 1:length(y)
            match_index = find(round(this.UniqueTemperatures, 2) == round(this.SetTemperatures(i), 2));
            if ~isempty(match_index)
                y(i) = this.FluoCoeffs(match_index)*y(i);
                y_err(i) = this.FluoCoeffs(match_index)*y_err(i);
            end
        end
    end
    
        
elseif strcmpi(parameter, 'burstdurations') | strcmpi(parameter, 'durations') 
    if ~UseScaledTiming
        y =  this.Durations(:, TimeIndex, APindex).';
        y_err =  this.DurationsStdErr(:, TimeIndex, APindex).';
    else
        y =  this.ScaledDurations(:, TimeIndex, APindex).';
        y_err =  this.ScaledDurationsStdErr(:, TimeIndex, APindex).';
    end

elseif strcmpi(parameter, 'burstfrequencies') | strcmpi(parameter, 'frequencies') 
    if ~UseScaledTiming
        y =  this.Frequencies(:, TimeIndex, APindex).';
        y_err =  this.FrequenciesStdErr(:, TimeIndex, APindex).';
    else
        y =  this.ScaledFrequencies(:, TimeIndex, APindex).';
        y_err =  this.ScaledFrequenciesStdErr(:, TimeIndex, APindex).';
    end

else
   error(['Parameter: ', parameter, ' not supported.'])
end

if ~all(isnan(y_err))
    IncludedSets = ((y./y_err > 1) & ~isnan(y) & (y > 0));
else
    IncludedSets = (~isnan(y)  & (y > 0));
end

x = x(IncludedSets);
y = y(IncludedSets);
y_err = y_err(IncludedSets);

x = 1./(R*(x+273.15));
y = log(y);
y_err = abs(y_err./y);






