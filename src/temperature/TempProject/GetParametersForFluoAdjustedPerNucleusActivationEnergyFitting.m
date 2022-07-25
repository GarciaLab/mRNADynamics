function [x, x_setpoint, y, y_err]  = GetParametersForFluoAdjustedPerNucleusActivationEnergyFitting(this, parameter, NC, APindex, TraceType)
x = this.Temp_obs;
x_setpoint = this.Temp_sps;
NumSets = length(x);
UseSet = ismember(1:NumSets, this.ProcessedExperiments);

SetFluoCoeffs = NaN(1, NumSets);
temperatures = this.UniqueTemperatures;
for i = 1:NumSets
    temp_match_index = find(round(temperatures, 2) == round(this.Temp_sps(i), 2));
    SetFluoCoeffs(i) = this.FluoCoeffs(temp_match_index);
    
end


if strcmpi(parameter, 'loadingrates')
    InitiationRates = getPerNucleusTrapezoidParameters(this, 'MeanInitiationRates', TraceType, false);
    SEInitiationRates = getPerNucleusTrapezoidParameters(this, 'MeanInitiationRates', TraceType, true);
    R2s = getPerNucleusTrapezoidParameters(this, 'R2s', TraceType);
    y = (InitiationRates(:,APindex, NC-8).').*SetFluoCoeffs;
    y_err = (SEInitiationRates(:,APindex,NC-8).').*SetFluoCoeffs;
    setR2s = R2s(:,APindex, NC-8).';
elseif strcmpi(parameter, 'maxfluos') | strcmpi(parameter, 'maxfluo') 
    [MaxFluos, SEMaxFluos] = getPerNucleusMaxFluoMatForPlotting(this, TraceType);
    R2s = ones(size(MaxFluos));
    y =  (MaxFluos(:,APindex, NC-8).').*SetFluoCoeffs;
    y_err = (SEMaxFluos(:, APindex,NC-8).').*SetFluoCoeffs;
    setR2s = R2s(:,APindex, NC-8).';
elseif strcmpi(parameter, 'plateauheights') | strcmpi(parameter, 'plateauheight')
    InitiationRates = getPerNucleusTrapezoidParameters(this, 'MeanInitiationRates', TraceType, false);
    SEInitiationRates = getPerNucleusTrapezoidParameters(this, 'MeanInitiationRates', TraceType, true);
    ElongationTimes = getPerNucleusTrapezoidParameters(this, 'ElongationTimes', TraceType, false);
    SEElongationTimes = getPerNucleusTrapezoidParameters(this, 'ElongationTimes', TraceType, true);
    PlateauHeights = InitiationRates.*ElongationTimes;
    SEPlateauHeights = sqrt((InitiationRates.^2).*(SEElongationTimes.^2)+...
        (ElongationTimes.^2).*(SEInitiationRates.^2) );
    R2s = getPerNucleusTrapezoidParameters(this, 'R2s', TraceType);
    y =  (PlateauHeights(:,APindex, NC-8).').*SetFluoCoeffs;
    y_err = (SEPlateauHeights(:, APindex,NC-8).').*SetFluoCoeffs;
    setR2s = R2s(:,APindex, NC-8).';
else
   error(['Parameter: ', parameter, ' not supported.'])
end

if ~all(isnan(y_err))
    IncludedSets = (UseSet & (y./y_err > 1) & ~isnan(y) & (setR2s >= this.R2bound) & (y > 0));
else
    IncludedSets = (UseSet & ~isnan(y) & (setR2s >= this.R2bound) & (y > 0));
end

x = x(IncludedSets);
x_setpoint = x_setpoint(IncludedSets);
y = y(IncludedSets);
y_err = y_err(IncludedSets);

x = 1./(this.R*(x+273));
y = log(y);
y_err = abs(y_err./y);






