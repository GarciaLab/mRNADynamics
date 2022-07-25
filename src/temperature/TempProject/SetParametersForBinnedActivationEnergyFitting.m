function this = SetParametersForBinnedActivationEnergyFitting(this, f, gof,...
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
    this.BinnedActivationEnergyFits.('CycleDurations').ActivationEnergies.General(NC-8) = Ea;
    this.BinnedActivationEnergyFits.('CycleDurations').ActivationEnergies.GeneralStdError(NC-8) = se_Ea;
    this.BinnedActivationEnergyFits.('CycleDurations').ActivationEnergies.GeneralCIlow(NC-8) = ci_Ea(1);
    this.BinnedActivationEnergyFits.('CycleDurations').ActivationEnergies.GeneralCIhigh(NC-8) = ci_Ea(2);
    
    this.BinnedActivationEnergyFits.('CycleDurations').LogScalingCoefficients.General(NC-8) = LogA;
    this.BinnedActivationEnergyFits.('CycleDurations').LogScalingCoefficients.GeneralStdError(NC-8) = se_LogA;
    this.BinnedActivationEnergyFits.('CycleDurations').LogScalingCoefficients.GeneralCIlow(NC-8) = ci_LogA(1);
    this.BinnedActivationEnergyFits.('CycleDurations').LogScalingCoefficients.GeneralCIhigh(NC-8) = ci_LogA(2);
    
    this.BinnedActivationEnergyFits.('CycleDurations').R2s.General(NC-8) = r2;
    
    this.BinnedActivationEnergyFits.('CycleDurations').Fits.General{NC-8} = f;
    
    
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
    
    this.BinnedActivationEnergyFits.(pstring).ActivationEnergies.(TraceType)(APindex, NC-8) = -Ea;
    this.BinnedActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'StdError'])(APindex, NC-8) = se_Ea;
    this.BinnedActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'CIlow'])(APindex, NC-8) = -ci_Ea(1);
    this.BinnedActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'CIlow'])(APindex, NC-8) = -ci_Ea(2);
    
    this.BinnedActivationEnergyFits.(pstring).LogScalingCoefficients.(TraceType)(APindex, NC-8) = LogA;
    this.BinnedActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'StdError'])(APindex, NC-8) = se_LogA;
    this.BinnedActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'CIlow'])(APindex, NC-8) = ci_LogA(1);
    this.BinnedActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'CIlow'])(APindex, NC-8) = ci_LogA(2);
    
    this.BinnedActivationEnergyFits.(pstring).R2s.(TraceType)(APindex, NC-8) = r2;
    
    this.BinnedActivationEnergyFits.(pstring).Fits.(TraceType){APindex, NC-8} = f;
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
    this.BinnedActivationEnergyFits.(pstring).ActivationEnergies.(TraceType)(APindex, NC-8) = Ea;
    this.BinnedActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'StdError'])(APindex, NC-8) = se_Ea;
    this.BinnedActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'CIlow'])(APindex, NC-8) = ci_Ea(1);
    this.BinnedActivationEnergyFits.(pstring).ActivationEnergies.([TraceType, 'CIlow'])(APindex, NC-8) = ci_Ea(2);
    
    this.BinnedActivationEnergyFits.(pstring).LogScalingCoefficients.(TraceType)(APindex, NC-8) = LogA;
    this.BinnedActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'StdError'])(APindex, NC-8) = se_LogA;
    this.BinnedActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'CIlow'])(APindex, NC-8) = ci_LogA(1);
    this.BinnedActivationEnergyFits.(pstring).LogScalingCoefficients.([TraceType, 'CIlow'])(APindex, NC-8) = ci_LogA(2);
    
    this.BinnedActivationEnergyFits.(pstring).R2s.(TraceType)(APindex, NC-8) = r2;
    
    this.BinnedActivationEnergyFits.(pstring).Fits.(TraceType){APindex, NC-8} = f;
end



end
