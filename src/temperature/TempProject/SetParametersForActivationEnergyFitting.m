function this = SetParametersForActivationEnergyFitting(this, f, gof,...
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
    this.ActivationEnergyFits.('CycleDurations').ActivationEnergies.General(NC-8) = Ea;
    this.ActivationEnergyFits.('CycleDurations').ActivationEnergies.GeneralStdError(NC-8) = se_Ea;
    this.ActivationEnergyFits.('CycleDurations').ActivationEnergies.GeneralCIlow(NC-8) = ci_Ea(1);
    this.ActivationEnergyFits.('CycleDurations').ActivationEnergies.GeneralCIhigh(NC-8) = ci_Ea(2);
    
    this.ActivationEnergyFits.('CycleDurations').LogScalingCoefficients.General(NC-8) = LogA;
    this.ActivationEnergyFits.('CycleDurations').LogScalingCoefficients.GeneralStdError(NC-8) = se_LogA;
    this.ActivationEnergyFits.('CycleDurations').LogScalingCoefficients.GeneralCIlow(NC-8) = ci_LogA(1);
    this.ActivationEnergyFits.('CycleDurations').LogScalingCoefficients.GeneralCIhigh(NC-8) = ci_LogA(2);
    
    this.ActivationEnergyFits.('CycleDurations').R2s.General(NC-8) = r2;
    
    this.ActivationEnergyFits.('CycleDurations').Fits.General{NC-8} = f;
    
    
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
    
    this.ActivationEnergyFits.(pstring).ActivationEnergies.(TraceType)(APindex, NC-8) = -Ea;
    this.ActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'StdError'])(APindex, NC-8) = se_Ea;
    this.ActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'CIlow'])(APindex, NC-8) = -ci_Ea(1);
    this.ActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'CIlow'])(APindex, NC-8) = -ci_Ea(2);
    
    this.ActivationEnergyFits.(pstring).LogScalingCoefficients.(TraceType)(APindex, NC-8) = LogA;
    this.ActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'StdError'])(APindex, NC-8) = se_LogA;
    this.ActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'CIlow'])(APindex, NC-8) = ci_LogA(1);
    this.ActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'CIlow'])(APindex, NC-8) = ci_LogA(2);
    
    this.ActivationEnergyFits.(pstring).R2s.(TraceType)(APindex, NC-8) = r2;
    
    this.ActivationEnergyFits.(pstring).Fits.(TraceType){APindex, NC-8} = f;
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
    this.ActivationEnergyFits.(pstring).ActivationEnergies.(TraceType)(APindex, NC-8) = Ea;
    this.ActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'StdError'])(APindex, NC-8) = se_Ea;
    this.ActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'CIlow'])(APindex, NC-8) = ci_Ea(1);
    this.ActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'CIlow'])(APindex, NC-8) = ci_Ea(2);
    
    this.ActivationEnergyFits.(pstring).LogScalingCoefficients.(TraceType)(APindex, NC-8) = LogA;
    this.ActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'StdError'])(APindex, NC-8) = se_LogA;
    this.ActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'CIlow'])(APindex, NC-8) = ci_LogA(1);
    this.ActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'CIlow'])(APindex, NC-8) = ci_LogA(2);
    
    this.ActivationEnergyFits.(pstring).R2s.(TraceType)(APindex, NC-8) = r2;
    
    this.ActivationEnergyFits.(pstring).Fits.(TraceType){APindex, NC-8} = f;
end



end
