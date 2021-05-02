function [Ea, se_Ea, LogA, se_LogA, R2, PlotTitle] = ...
    getActivationEnergyMatrices(this, parameter, TraceType, UseFluoAdjusted)
if ~exist('UseFluoAdjusted', 'var')
    UseFluoAdjusted = false;
end
if ~UseFluoAdjusted
    if ~strcmpi(parameter, 'cycledurations')
        Ea = this.ActivationEnergyFits.(parameter).ActivationEnergies.(TraceType);
        se_Ea = this.ActivationEnergyFits.(parameter).ActivationEnergies.([TraceType, 'StdError']);
        LogA = this.ActivationEnergyFits.(parameter).LogScalingCoefficients.(TraceType);
        se_LogA = this.ActivationEnergyFits.(parameter).LogScalingCoefficients.([TraceType, 'StdError']);
        
        R2 = this.ActivationEnergyFits.(parameter).R2s.(TraceType);
        if strcmpi(parameter, 'TimeOns') | strcmpi(parameter, 'TimeOn')
            PlotTitle = 't_{on}';
        elseif strcmpi(parameter, 'TimeOffs') | strcmpi(parameter, 'TimeOff')
            PlotTitle = 't_{off}';
        elseif strcmpi(parameter, 'TranscriptionWindows')| strcmpi(parameter, 'TranscriptionWindow')
            PlotTitle = 'Transcription Window Duration';
        elseif strcmpi(parameter, 'ElongationTimes') | strcmpi(parameter, 'ElongationTime')
            PlotTitle = 'Elongation Time';
        elseif strcmpi(parameter, 'ElongationRates') | strcmpi(parameter, 'ElongationRate')
            PlotTitle = 'Elongation Rate';
        elseif strcmpi(parameter, 'LoadingRates') | strcmpi(parameter, 'LoadingRate')
            PlotTitle = 'Loading Rate';
        elseif strcmpi(parameter, 'PostTranscriptionDurations') | strcmpi(parameter, 'PostTranscriptionDuration')
            PlotTitle = 'Post-transcriptional Duration';
        elseif strcmpi(parameter, 'PlateauHeights') | strcmpi(parameter, 'PlateauHeight')
            PlotTitle = 'Fitted Plateau Height';
        elseif strcmpi(parameter, 'MaxFluos') | strcmpi(parameter, 'MaxFluo')
            PlotTitle = 'Maximum Fluorescence';
        end
        
    else
        Ea = this.ActivationEnergyFits.(parameter).ActivationEnergies.General;
        se_Ea = this.ActivationEnergyFits.(parameter).ActivationEnergies.GeneralStdError;
        LogA = this.ActivationEnergyFits.(parameter).LogScalingCoefficients.General;
        se_LogA = this.ActivationEnergyFits.(parameter).LogScalingCoefficients.GeneralStdError;
        
        R2 = this.ActivationEnergyFits.(parameter).R2s.General;
        PlotTitle = 'Division Cycle Durations';
        
        
    end
else
    Ea = this.FluoAdjustedActivationEnergyFits.(parameter).ActivationEnergies.(TraceType);
    se_Ea = this.FluoAdjustedActivationEnergyFits.(parameter).ActivationEnergies.([TraceType, 'StdError']);
    LogA = this.FluoAdjustedActivationEnergyFits.(parameter).LogScalingCoefficients.(TraceType);
    se_LogA = this.FluoAdjustedActivationEnergyFits.(parameter).LogScalingCoefficients.([TraceType, 'StdError']);
    
    R2 = this.FluoAdjustedActivationEnergyFits.(parameter).R2s.(TraceType);
    if strcmpi(parameter, 'LoadingRates') | strcmpi(parameter, 'LoadingRate')
        PlotTitle = 'Re-scaled Loading Rate';
    elseif strcmpi(parameter, 'PlateauHeights') | strcmpi(parameter, 'PlateauHeight')
        PlotTitle = 'Re-scaled Fitted Plateau Height';
    elseif strcmpi(parameter, 'MaxFluos') | strcmpi(parameter, 'MaxFluo')
        PlotTitle = 'Re-scaled Maximum Fluorescence';
    end   
end