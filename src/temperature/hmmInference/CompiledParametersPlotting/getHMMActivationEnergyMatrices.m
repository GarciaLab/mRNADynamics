function [Ea, se_Ea, LogA, se_LogA, R2, PlotTitle] = ...
    getHMMActivationEnergyMatrices(this, parameter, UseRescaledFluo, UseRescaledTiming)
if ~exist('UseRescaledFluo', 'var')
    UseRescaledFluo = false;
end
if ~exist('UseRescaledTiming', 'var')
    UseRescaledTiming = false;
end

if UseRescaledFluo & UseRescaledTiming
        
    Ea = this.DoubleScaledActivationEnergyFits.(parameter).ActivationEnergies;
    se_Ea = this.DoubleScaledActivationEnergyFits.(parameter).SEActivationEnergies;
    LogA = this.DoubleScaledActivationEnergyFits.(parameter).LogScalingCoefficients;
    se_LogA = this.DoubleScaledActivationEnergyFits.(parameter).SELogScalingCoefficients;
    R2 = this.DoubleScaledActivationEnergyFits.(parameter).R2s;
    
    if strcmpi(parameter, 'initiationrates')
        PlotTitle = 'Fluo-Adjusted Initiation Rates (Time Bins Rescaled)';
    elseif strcmpi(parameter, 'durations')
        PlotTitle = 'Burst Durations (Time Bins Rescaled)';
    else
        PlotTitle = 'Burst Frequencies (Time Bins Rescaled)';
    end
elseif UseRescaledTiming
    
    Ea = this.TimeScaledActivationEnergyFits.(parameter).ActivationEnergies;
    se_Ea = this.TimeScaledActivationEnergyFits.(parameter).SEActivationEnergies;
    LogA = this.TimeScaledActivationEnergyFits.(parameter).LogScalingCoefficients;
    se_LogA = this.TimeScaledActivationEnergyFits.(parameter).SELogScalingCoefficients;
    R2 = this.TimeScaledActivationEnergyFits.(parameter).R2s;
    
    if strcmpi(parameter, 'initiationrates')
        PlotTitle = 'Fluo-Adjusted Initiation Rates';
    elseif strcmpi(parameter, 'durations')
        PlotTitle = 'Burst Durations';
    else
        PlotTitle = 'Burst Frequencies';
    end
elseif UseRescaledFluo
    
    Ea = this.FluoScaledActivationEnergyFits.(parameter).ActivationEnergies;
    se_Ea = this.FluoScaledActivationEnergyFits.(parameter).SEActivationEnergies;
    LogA = this.FluoScaledActivationEnergyFits.(parameter).LogScalingCoefficients;
    se_LogA = this.FluoScaledActivationEnergyFits.(parameter).SELogScalingCoefficients;
    R2 = this.FluoScaledActivationEnergyFits.(parameter).R2s;
    
    if strcmpi(parameter, 'initiationrates')
        PlotTitle = 'Initiation Rates (Time Bins Rescaled)';
    elseif strcmpi(parameter, 'durations')
        PlotTitle = 'Burst Durations (Time Bins Rescaled)';
    else
        PlotTitle = 'Burst Frequencies (Time Bins Rescaled)';
    end

else
    
    Ea = this.ActivationEnergyFits.(parameter).ActivationEnergies;
    se_Ea = this.ActivationEnergyFits.(parameter).SEActivationEnergies;
    LogA = this.ActivationEnergyFits.(parameter).LogScalingCoefficients;
    se_LogA = this.ActivationEnergyFits.(parameter).SELogScalingCoefficients;
    R2 = this.ActivationEnergyFits.(parameter).R2s;
    
    if strcmpi(parameter, 'initiationrates')
        PlotTitle = 'Initiation Rates';
    elseif strcmpi(parameter, 'durations')
        PlotTitle = 'Burst Durations';
    else
        PlotTitle = 'Burst Frequencies';
    end


end

