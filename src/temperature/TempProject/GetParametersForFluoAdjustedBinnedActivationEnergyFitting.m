function [x, x_setpoint, y, y_err]  = GetParametersForFluoAdjustedBinnedActivationEnergyFitting(this, parameter, NC, APindex, TraceType)
Temperatures = flip(unique(this.Temp_sps));
x = Temperatures;
x_setpoint = Temperatures;
NumTemperatures = length(x);


SetFluoCoeffs = this.FluoCoeffs;



if strcmpi(parameter, 'loadingrates')
    InitiationRates = getBinnedTrapezoidParameters(this, 'MeanInitiationRates', TraceType, false);
    SEInitiationRates = getBinnedTrapezoidParameters(this, 'MeanInitiationRates', TraceType, true);
    R2s = getBinnedTrapezoidParameters(this, 'R2s', TraceType);
    y = (InitiationRates(:,APindex, NC-8).').*SetFluoCoeffs;
    y_err = (SEInitiationRates(:,APindex,NC-8).').*SetFluoCoeffs;
    setR2s = R2s(:,APindex, NC-8).';
elseif strcmpi(parameter, 'maxfluos') | strcmpi(parameter, 'maxfluo') 
    [MaxFluos, SEMaxFluos] = getBinnedMaxFluoMatForPlotting(this, TraceType);
    R2s = ones(size(MaxFluos));
    y =  (MaxFluos(:,APindex, NC-8).').*SetFluoCoeffs;
    y_err = (SEMaxFluos(:, APindex,NC-8).').*SetFluoCoeffs;
    setR2s = R2s(:,APindex, NC-8).';
elseif strcmpi(parameter, 'plateauheights') | strcmpi(parameter, 'plateauheight')
    InitiationRates = getBinnedTrapezoidParameters(this, 'MeanInitiationRates', TraceType, false);
    SEInitiationRates = getBinnedTrapezoidParameters(this, 'MeanInitiationRates', TraceType, true);
    ElongationTimes = getBinnedTrapezoidParameters(this, 'ElongationTimes', TraceType, false);
    SEElongationTimes = getBinnedTrapezoidParameters(this, 'ElongationTimes', TraceType, true);
    PlateauHeights = InitiationRates.*ElongationTimes;
    SEPlateauHeights = sqrt((InitiationRates.^2).*(SEElongationTimes.^2)+...
        (ElongationTimes.^2).*(SEInitiationRates.^2) );
    R2s = getBinnedTrapezoidParameters(this, 'R2s', TraceType);
    y =  (PlateauHeights(:,APindex, NC-8).').*SetFluoCoeffs;
    y_err = (SEPlateauHeights(:, APindex,NC-8).').*SetFluoCoeffs;
    setR2s = R2s(:,APindex, NC-8).';
else
   error(['Parameter: ', parameter, ' not supported.'])
end

if ~all(isnan(y_err))
    IncludedSets = ((y./y_err > 1) & ~isnan(y) & (setR2s >= this.R2bound) & (y > 0));
else
    IncludedSets = (~isnan(y) & (setR2s >= this.R2bound) & (y > 0));
end

x = x(IncludedSets);
x_setpoint = x_setpoint(IncludedSets);
y = y(IncludedSets);
y_err = y_err(IncludedSets);

x = 1./(this.R*(x+273));
y = log(y);
y_err = abs(y_err./y);






