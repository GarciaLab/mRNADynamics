function [x, y, ci, Ea, se_Ea, LogA, se_LogA, R2] = ...
    getActivationEnergyFitTraces(this, parameter, NC, APindex, TraceType)
%% Load relevant parameters into memory
[Ea, se_Ea, LogA, se_LogA, R2, fitresult] = ...
    getFittedActivationEnergyParameters(this, parameter, NC, APindex, TraceType);

if all(~isnan([Ea, LogA]))
    t = 15:.01:30;
    x = 1./(this.R*(t+273));
    try
        [ci, y] = predint(fitresult, x, this.alpha, 'functional', 'off');
        y = exp(y);
        ci = exp(ci);
    catch 
        y = [];
        ci = [];
    end
else
    x = [];
    y = [];
    ci = [];
end

