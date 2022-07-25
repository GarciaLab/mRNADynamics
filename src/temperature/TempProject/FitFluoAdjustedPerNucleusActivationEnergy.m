function this = FitFluoAdjustedPerNucleusActivationEnergy(this, parameter,  NC, APindex, TraceType)
if exist('APindex', 'var') & exist('TraceType', 'var')
    [t, t_setpoint, y, y_err]  = GetParametersForFluoAdjustedPerNucleusActivationEnergyFitting(this, parameter, NC, APindex, TraceType);
else
    [t, t_setpoint, y, y_err]  = GetParametersForFluoAdjustedPerNucleusActivationEnergyFitting(this, parameter, NC);
end




if length(unique(t_setpoint)) >= 3 &  length(t_setpoint) >= this.MinimumFittingPoints
    [f, gof] = fit( t.', y.', 'poly1');
    
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
    
    if exist('APindex', 'var') & exist('TraceType', 'var')
        this = SetParametersForFluoAdjustedPerNucleusActivationEnergyFitting(this, f, gof,...
            parameter, NC, APindex, TraceType);
    else
        this = SetParametersForFluoAdjustedPerNucleusActivationEnergyFitting(this, f, gof, parameter, NC);
    end
end