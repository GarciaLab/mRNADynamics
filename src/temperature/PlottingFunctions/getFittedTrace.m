function [t_vector, fit_solution, ci, pos_slope, se_pos_slope, neg_slope,...
    se_neg_slope, time_on, se_time_on, time_off, se_time_off, time_peak, se_time_peak, R2] = ...
    getFittedTrace(this, SetIndex, APindex, NC, TraceType, MaxTime, UseBinnedParameters)
if ~exist('UseBinnedParameters', 'var')
    UseBinnedParameters = false;
end
[pos_slope, se_pos_slope, neg_slope, se_neg_slope, time_on, se_time_on,...
    time_off, se_time_off, time_peak, se_time_peak, R2, fitresult] = ...
    getFittedTrapezoidParameters(this, SetIndex, APindex, NC, TraceType, UseBinnedParameters);
if all(~isnan([pos_slope, neg_slope, time_on, time_off, time_peak]))
    plateau_height = pos_slope*time_peak+time_on;
    time_dark = (neg_slope*time_off-plateau_height)/neg_slope;
    pos_yintercept = -pos_slope*time_on;
    if NC < 14
        t_vector = time_on:0.1:min(time_dark, MaxTime);
    else
        t_vector = time_on:0.1:min(time_off, MaxTime);
    end
    fit_solution = trapezoidFitFunction(t_vector, pos_slope, pos_yintercept, neg_slope, time_peak, time_off);
    try
        ci = predint(fitresult, t_vector, this.alpha, 'functional', 'off');
    catch 
        ci = [];
    end
elseif all(~isnan([pos_slope, time_on, time_off, time_peak]))
    pos_yintercept = -pos_slope*time_on;
    t_vector = time_on:0.1:min(time_off, MaxTime);
    fit_solution = trapezoidFitFunction(t_vector, pos_slope, pos_yintercept, neg_slope, time_peak, time_off);
    try
        ci = predint(fitresult, t_vector, this.alpha, 'functional', 'off');
    catch 
        ci = [];
    end
elseif all(~isnan([pos_slope, time_on, time_peak]))
    pos_yintercept = -pos_slope*time_on;
    t_vector = time_on:0.1:min(time_peak, MaxTime);
    fit_solution = trapezoidFitFunction(t_vector, pos_slope, pos_yintercept, neg_slope, time_peak, time_off);
    try
        ci = predint(fitresult, t_vector, this.alpha, 'functional', 'off');
    catch 
        ci = [];
    end
else
    t_vector = [];
    fit_solution = [];
    ci = [];
end

