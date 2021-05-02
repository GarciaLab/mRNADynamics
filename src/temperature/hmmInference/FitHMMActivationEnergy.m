function this = FitHMMActivationEnergy(this, parameter, UseScaledFluo, UseScaledTiming, TimeIndex, APindex)

[x, y, y_err]  = GetHmmParametersForActivationEnergyFitting(this, parameter,  UseScaledFluo, UseScaledTiming, TimeIndex, APindex);
alpha =0.95;


if length(x) >= 4
    [f, gof] = fit( x.', y.', 'poly1');
    
    try
        ci = confint(f);
    catch
        ci = NaN(2,2);
    end
    t = tinv((1+alpha)/2, 1);
    
    Ea = f.p1;
    ci_Ea = ci(:,1).';
    se_Ea = (ci_Ea(2)-ci_Ea(1)) ./ (2*t); % Standard Error
    
    LogA = f.p2;
    ci_LogA = ci(:,2).';
    se_LogA = (ci_LogA(2)-ci_LogA(1)) ./ (2*t); % Standard Error
    
    
    this = SetHmmParametersForActivationEnergyFitting(this, f, gof,...
        parameter, UseScaledFluo, UseScaledTiming, TimeIndex, APindex);

end