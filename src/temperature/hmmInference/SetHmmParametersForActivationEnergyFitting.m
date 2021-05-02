function this = SetHmmParametersForActivationEnergyFitting(this, f, gof,...
    parameter, UseScaledFluo, UseScaledTiming, TimeIndex, APindex)
try
    ci = confint(f);
catch
    ci = NaN(2,2);
end
alpha = 0.95;

t = tinv((1+alpha)/2, 1);

Ea = f.p1;
ci_Ea = ci(:,1).';
se_Ea = (ci_Ea(2)-ci_Ea(1)) ./ (2*t); % Standard Error

LogA = f.p2;
ci_LogA = ci(:,2).';
se_LogA = (ci_LogA(2)-ci_LogA(1)) ./ (2*t); % Standard Error

r2 = gof.rsquare;


Ea = -Ea;

if UseScaledFluo & UseScaledTiming
    this.DoubleScaledActivationEnergyFits.(parameter).ActivationEnergies(TimeIndex, APindex) = Ea;
    this.DoubleScaledActivationEnergyFits.(parameter).SEActivationEnergies(TimeIndex, APindex) = se_Ea;
    this.DoubleScaledActivationEnergyFits.(parameter).LogScalingCoefficients(TimeIndex, APindex) = LogA;
    this.DoubleScaledActivationEnergyFits.(parameter).SELogScalingCoefficients(TimeIndex, APindex) = se_LogA;
    this.DoubleScaledActivationEnergyFits.(parameter).R2s(TimeIndex, APindex) = r2;
    this.DoubleScaledActivationEnergyFits.(parameter).Fits{TimeIndex, APindex} = f;
elseif UseScaledFluo
    this.FluoScaledActivationEnergyFits.(parameter).ActivationEnergies(TimeIndex, APindex) = Ea;
    this.FluoScaledActivationEnergyFits.(parameter).SEActivationEnergies(TimeIndex, APindex) = se_Ea;
    this.FluoScaledActivationEnergyFits.(parameter).LogScalingCoefficients(TimeIndex, APindex) = LogA;
    this.FluoScaledActivationEnergyFits.(parameter).SELogScalingCoefficients(TimeIndex, APindex) = se_LogA;
    this.FluoScaledActivationEnergyFits.(parameter).R2s(TimeIndex, APindex) = r2;
    this.FluoScaledActivationEnergyFits.(parameter).Fits{TimeIndex, APindex} = f;
elseif UseScaledTiming
    this.TimeScaledActivationEnergyFits.(parameter).ActivationEnergies(TimeIndex, APindex) = Ea;
    this.TimeScaledActivationEnergyFits.(parameter).SEActivationEnergies(TimeIndex, APindex) = se_Ea;
    this.TimeScaledActivationEnergyFits.(parameter).LogScalingCoefficients(TimeIndex, APindex) = LogA;
    this.TimeScaledActivationEnergyFits.(parameter).SELogScalingCoefficients(TimeIndex, APindex) = se_LogA;
    this.TimeScaledActivationEnergyFits.(parameter).R2s(TimeIndex, APindex) = r2;
    this.TimeScaledActivationEnergyFits.(parameter).Fits{TimeIndex, APindex} = f;
else
    this.ActivationEnergyFits.(parameter).ActivationEnergies(TimeIndex, APindex) = Ea;
    this.ActivationEnergyFits.(parameter).SEActivationEnergies(TimeIndex, APindex) = se_Ea;
    this.ActivationEnergyFits.(parameter).LogScalingCoefficients(TimeIndex, APindex) = LogA;
    this.ActivationEnergyFits.(parameter).SELogScalingCoefficients(TimeIndex, APindex) = se_LogA;
    this.ActivationEnergyFits.(parameter).R2s(TimeIndex, APindex) = r2;
    this.ActivationEnergyFits.(parameter).Fits{TimeIndex, APindex} = f;
end


end
