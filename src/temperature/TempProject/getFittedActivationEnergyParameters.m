function [Ea, se_Ea, LogA, se_LogA, R2, fitresult] = ...
    getFittedActivationEnergyParameters(this, parameter, NC, APindex, TraceType)
if ~strcmpi(parameter, 'cycledurations')
    Ea = this.ActivationEnergyFits.(parameter).ActivationEnergies.(TraceType)(APindex, NC-8);
    se_Ea = this.ActivationEnergyFits.(parameter).ActivationEnergies.([TraceType, 'StdError'])(APindex, NC-8); 
    LogA = this.ActivationEnergyFits.(parameter).LogScalingCoefficients.(TraceType)(APindex, NC-8);
    se_LogA = this.ActivationEnergyFits.(parameter).LogScalingCoefficients.([TraceType, 'StdError'])(APindex, NC-8);
    
    R2 = this.ActivationEnergyFits.(parameter).R2s.(TraceType)(APindex, NC-8);
    
    fitresult = this.ActivationEnergyFits.(parameter).Fits.(TraceType){APindex, NC-8};
else
    Ea = this.ActivationEnergyFits.(parameter).ActivationEnergies.General(NC-8);
    se_Ea = this.ActivationEnergyFits.(parameter).ActivationEnergies.GeneralStdError(NC-8); 
    LogA = this.ActivationEnergyFits.(parameter).LogScalingCoefficients.General(NC-8);
    se_LogA = this.ActivationEnergyFits.(parameter).LogScalingCoefficients.GeneralStdError(NC-8); 
    
    R2 = this.ActivationEnergyFits.(parameter).R2s.General(NC-8);
    
    fitresult = this.ActivationEnergyFits.(parameter).Fits.General{NC-8};
end