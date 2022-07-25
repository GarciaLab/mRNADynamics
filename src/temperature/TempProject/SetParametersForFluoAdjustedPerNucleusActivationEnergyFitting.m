function this = SetParametersForFluoAdjustedPerNucleusActivationEnergyFitting(this, f, gof,...
    parameter, NC, APindex, TraceType)
try
    ci = confint(f);
catch
    ci = NaN(2,2);
end
t = tinv((1+this.alpha)/2, 1);

Ea = f.p1;
ci_Ea = ci(:,1).';
se_Ea = (ci_Ea(2)-ci_Ea(1)) ./ (2*t); % Standard Error

LogA = f.p2;
ci_LogA = ci(:,2).';
se_LogA = (ci_LogA(2)-ci_LogA(1)) ./ (2*t); % Standard Error

r2 = gof.rsquare;

if strcmpi(parameter, 'loadingrates') | strcmpi(parameter, 'plateauheights') | strcmpi(parameter, 'maxfluos')
    if strcmpi(parameter, 'loadingrates')
        pstring = 'LoadingRates';
    elseif strcmpi(parameter, 'plateauheights')
        pstring = 'PlateauHeights';
    elseif strcmpi(parameter, 'maxfluos')
        pstring = 'MaxFluos';
    end
    
    this.FluoAdjustedPerNucleusActivationEnergyFits.(pstring).ActivationEnergies.(TraceType)(APindex, NC-8) = -Ea;
    this.FluoAdjustedPerNucleusActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'StdError'])(APindex, NC-8) = se_Ea;
    this.FluoAdjustedPerNucleusActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'CIlow'])(APindex, NC-8) = -ci_Ea(1);
    this.FluoAdjustedPerNucleusActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'CIlow'])(APindex, NC-8) = -ci_Ea(2);
    
    this.FluoAdjustedPerNucleusActivationEnergyFits.(pstring).LogScalingCoefficients.(TraceType)(APindex, NC-8) = LogA;
    this.FluoAdjustedPerNucleusActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'StdError'])(APindex, NC-8) = se_LogA;
    this.FluoAdjustedPerNucleusActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'CIlow'])(APindex, NC-8) = ci_LogA(1);
    this.FluoAdjustedPerNucleusActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'CIlow'])(APindex, NC-8) = ci_LogA(2);
    
    this.FluoAdjustedPerNucleusActivationEnergyFits.(pstring).R2s.(TraceType)(APindex, NC-8) = r2;
    
    this.FluoAdjustedPerNucleusActivationEnergyFits.(pstring).Fits.(TraceType){APindex, NC-8} = f;
else
    error('Invalid Parameter Choice.');
end



end
