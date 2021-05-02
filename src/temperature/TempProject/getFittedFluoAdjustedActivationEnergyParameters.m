function [Ea, se_Ea, LogA, se_LogA, R2, fitresult] = ...
    getFittedFluoAdjustedActivationEnergyParameters(this, parameter, NC, APindex, TraceType)

    Ea = this.FluoAdjustedActivationEnergyFits.(parameter).ActivationEnergies.(TraceType)(APindex, NC-8);
    se_Ea = this.FluoAdjustedActivationEnergyFits.(parameter).ActivationEnergies.([TraceType, 'StdError'])(APindex, NC-8); 
    LogA = this.FluoAdjustedActivationEnergyFits.(parameter).LogScalingCoefficients.(TraceType)(APindex, NC-8);
    se_LogA = this.FluoAdjustedActivationEnergyFits.(parameter).LogScalingCoefficients.([TraceType, 'StdError'])(APindex, NC-8);
    
    R2 = this.FluoAdjustedActivationEnergyFits.(parameter).R2s.(TraceType)(APindex, NC-8);
    
    fitresult = this.FluoAdjustedActivationEnergyFits.(parameter).Fits.(TraceType){APindex, NC-8};
end