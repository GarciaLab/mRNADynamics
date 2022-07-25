function this = SetParametersForPerNucleusActivationEnergyFitting(this, f, gof,...
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

if strcmpi(parameter, 'cycledurations')
    this.PerNucleusActivationEnergyFits.('CycleDurations').ActivationEnergies.General(NC-8) = Ea;
    this.PerNucleusActivationEnergyFits.('CycleDurations').ActivationEnergies.GeneralStdError(NC-8) = se_Ea;
    this.PerNucleusActivationEnergyFits.('CycleDurations').ActivationEnergies.GeneralCIlow(NC-8) = ci_Ea(1);
    this.PerNucleusActivationEnergyFits.('CycleDurations').ActivationEnergies.GeneralCIhigh(NC-8) = ci_Ea(2);
    
    this.PerNucleusActivationEnergyFits.('CycleDurations').LogScalingCoefficients.General(NC-8) = LogA;
    this.PerNucleusActivationEnergyFits.('CycleDurations').LogScalingCoefficients.GeneralStdError(NC-8) = se_LogA;
    this.PerNucleusActivationEnergyFits.('CycleDurations').LogScalingCoefficients.GeneralCIlow(NC-8) = ci_LogA(1);
    this.PerNucleusActivationEnergyFits.('CycleDurations').LogScalingCoefficients.GeneralCIhigh(NC-8) = ci_LogA(2);
    
    this.PerNucleusActivationEnergyFits.('CycleDurations').R2s.General(NC-8) = r2;
    
    this.PerNucleusActivationEnergyFits.('CycleDurations').Fits.General{NC-8} = f;
    
    
elseif strcmpi(parameter, 'loadingrates') | strcmpi(parameter, 'elongationrates') | ...
        strcmpi(parameter, 'plateauheights') | strcmpi(parameter, 'maxfluos')
    if strcmpi(parameter, 'elongationrates')
        pstring = 'ElongationRates';
    elseif strcmpi(parameter, 'loadingrates')
        pstring = 'LoadingRates';
    elseif strcmpi(parameter, 'plateauheights')
        pstring = 'PlateauHeights';
    elseif strcmpi(parameter, 'maxfluos')
        pstring = 'MaxFluos';
    end
    
    this.PerNucleusActivationEnergyFits.(pstring).ActivationEnergies.(TraceType)(APindex, NC-8) = -Ea;
    this.PerNucleusActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'StdError'])(APindex, NC-8) = se_Ea;
    this.PerNucleusActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'CIlow'])(APindex, NC-8) = -ci_Ea(1);
    this.PerNucleusActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'CIlow'])(APindex, NC-8) = -ci_Ea(2);
    
    this.PerNucleusActivationEnergyFits.(pstring).LogScalingCoefficients.(TraceType)(APindex, NC-8) = LogA;
    this.PerNucleusActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'StdError'])(APindex, NC-8) = se_LogA;
    this.PerNucleusActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'CIlow'])(APindex, NC-8) = ci_LogA(1);
    this.PerNucleusActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'CIlow'])(APindex, NC-8) = ci_LogA(2);
    
    this.PerNucleusActivationEnergyFits.(pstring).R2s.(TraceType)(APindex, NC-8) = r2;
    
    this.PerNucleusActivationEnergyFits.(pstring).Fits.(TraceType){APindex, NC-8} = f;
else
    if strcmpi(parameter, 'timeons')
        pstring = 'TimeOns';
    elseif strcmpi(parameter, 'elongationtimes')
        pstring = 'ElongationTimes';
    elseif strcmpi(parameter, 'transcriptionwindows')
        pstring = 'TranscriptionWindows';
    elseif strcmpi(parameter, 'posttranscriptiondurations')
        pstring = 'PostTranscriptionDurations';
    end
    this.PerNucleusActivationEnergyFits.(pstring).ActivationEnergies.(TraceType)(APindex, NC-8) = Ea;
    this.PerNucleusActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'StdError'])(APindex, NC-8) = se_Ea;
    this.PerNucleusActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'CIlow'])(APindex, NC-8) = ci_Ea(1);
    this.PerNucleusActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'CIlow'])(APindex, NC-8) = ci_Ea(2);
    
    this.PerNucleusActivationEnergyFits.(pstring).LogScalingCoefficients.(TraceType)(APindex, NC-8) = LogA;
    this.PerNucleusActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'StdError'])(APindex, NC-8) = se_LogA;
    this.PerNucleusActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'CIlow'])(APindex, NC-8) = ci_LogA(1);
    this.PerNucleusActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'CIlow'])(APindex, NC-8) = ci_LogA(2);
    
    this.PerNucleusActivationEnergyFits.(pstring).R2s.(TraceType)(APindex, NC-8) = r2;
    
    this.PerNucleusActivationEnergyFits.(pstring).Fits.(TraceType){APindex, NC-8} = f;
end



end
