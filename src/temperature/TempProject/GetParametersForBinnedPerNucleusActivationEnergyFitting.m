function [x, x_setpoint, y, y_err]  = GetParametersForBinnedPerNucleusActivationEnergyFitting(this, parameter, NC, APindex, TraceType)
Temperatures = flip(unique(this.Temp_sps));
x =  Temperatures;
x_setpoint = Temperatures;
NumTemperatures = length(x);

if strcmpi(parameter, 'cycledurations')
    CycleDurations =  this.TimeScalingInfo.MeanNCDivisionInfo;
    SECycleDurations = this.TimeScalingInfo.StdNCDivisionInfo;
    y = CycleDurations(:,NC-8).';
    y_err = SECycleDurations(:,NC-8).';
    setR2s = ones(size(y));

elseif strcmpi(parameter, 'timeons')
    TimeOns = getBinnedPerNucleusTrapezoidParameters(this, 'TimeOns', TraceType, false);
    SETimeOns = getBinnedPerNucleusTrapezoidParameters(this, 'TimeOns', TraceType, true);
    R2s = getBinnedPerNucleusTrapezoidParameters(this, 'R2s', TraceType);
    y = TimeOns(:,APindex, NC-8).';
    y_err = SETimeOns(:,APindex,NC-8).';
    setR2s = R2s(:,APindex, NC-8).';

elseif strcmpi(parameter, 'elongationtimes')
    ElongationTimes = getBinnedPerNucleusTrapezoidParameters(this, 'ElongationTimes', TraceType, false);
    SEElongationTimes = getBinnedPerNucleusTrapezoidParameters(this, 'ElongationTimes', TraceType, true);
    R2s = getBinnedPerNucleusTrapezoidParameters(this, 'R2s', TraceType);
    y = ElongationTimes(:,APindex, NC-8).';
    y_err = SEElongationTimes(:,APindex, NC-8).';
    setR2s = R2s(:,APindex, NC-8).';

elseif strcmpi(parameter, 'elongationrates')
    ElongationTimes = getBinnedPerNucleusTrapezoidParameters(this, 'ElongationTimes', TraceType, false);
    SEElongationTimes = getBinnedPerNucleusTrapezoidParameters(this, 'ElongationTimes', TraceType, true);
    R2s = getBinnedPerNucleusTrapezoidParameters(this, 'R2s', TraceType);
    y = this.GeneLength./ElongationTimes(:,APindex, NC-8).';
    y_err = sqrt((((this.GeneLength./ElongationTimes(:,APindex, NC-8)).^2)./(ElongationTimes(:,APindex, NC-8).^2)).*SEElongationTimes(:,APindex, NC-8).^2).';
    setR2s = R2s(:,APindex, NC-8).';
    
elseif strcmpi(parameter, 'transcriptionwindows')
    TimeOns = getBinnedPerNucleusTrapezoidParameters(this, 'TimeOns', TraceType, false);
    SETimeOns = getBinnedPerNucleusTrapezoidParameters(this, 'TimeOns', TraceType, true);
    TimeOffs = getBinnedPerNucleusTrapezoidParameters(this, 'TimeOffs', TraceType, false);
    SETimeOffs = getBinnedPerNucleusTrapezoidParameters(this, 'TimeOffs', TraceType, true);
    R2s = getBinnedPerNucleusTrapezoidParameters(this, 'R2s', TraceType);
    y = TimeOffs(:,APindex, NC-8).'-TimeOns(:,APindex, NC-8).';
    y_err = sqrt(SETimeOns(:,APindex,NC-8).^2+SETimeOffs(:,APindex,NC-8).^2).';
    setR2s = R2s(:,APindex, NC-8).';

elseif strcmpi(parameter, 'posttranscriptiondurations')
    CycleDurations =  this.TimeScalingInfo.MeanNCDivisionInfo;
    SECycleDurations = this.TimeScalingInfo.StdNCDivisionInfo;
    
    TimeOffs = getBinnedPerNucleusTrapezoidParameters(this, 'TimeOffs', TraceType, false);
    SETimeOffs = getBinnedPerNucleusTrapezoidParameters(this, 'TimeOffs', TraceType, true);
    R2s = getBinnedPerNucleusTrapezoidParameters(this, 'R2s', TraceType);
    y = CycleDurations(:, NC-8).'-TimeOffs(:,APindex, NC-8).';
    y_err = sqrt(SECycleDurations(:,NC-8).^2+SETimeOffs(:,APindex,NC-8).^2).';
    setR2s = R2s(:,APindex, NC-8).'; 
    
elseif strcmpi(parameter, 'loadingrates')
    InitiationRates = getBinnedPerNucleusTrapezoidParameters(this, 'MeanInitiationRates', TraceType, false);
    SEInitiationRates = getBinnedPerNucleusTrapezoidParameters(this, 'MeanInitiationRates', TraceType, true);
    R2s = getBinnedPerNucleusTrapezoidParameters(this, 'R2s', TraceType);
    y = InitiationRates(:,APindex, NC-8).';
    y_err = SEInitiationRates(:,APindex,NC-8).';
    setR2s = R2s(:,APindex, NC-8).';
elseif strcmpi(parameter, 'maxfluos') | strcmpi(parameter, 'maxfluo') 
    [MaxFluos, SEMaxFluos] = getBinnedPerNucleusMaxFluoMatForPlotting(this, TraceType);
    R2s = ones(size(MaxFluos));
    y =  MaxFluos(:,APindex, NC-8).';
    y_err = SEMaxFluos(:, APindex,NC-8).';
    setR2s = R2s(:,APindex, NC-8).';
elseif strcmpi(parameter, 'plateauheights') | strcmpi(parameter, 'plateauheight')
    InitiationRates = getBinnedPerNucleusTrapezoidParameters(this, 'MeanInitiationRates', TraceType, false);
    SEInitiationRates = getBinnedPerNucleusTrapezoidParameters(this, 'MeanInitiationRates', TraceType, true);
    ElongationTimes = getBinnedPerNucleusTrapezoidParameters(this, 'ElongationTimes', TraceType, false);
    SEElongationTimes = getBinnedPerNucleusTrapezoidParameters(this, 'ElongationTimes', TraceType, true);
    PlateauHeights = InitiationRates.*ElongationTimes;
    SEPlateauHeights = sqrt((InitiationRates.^2).*(SEElongationTimes.^2)+...
        (ElongationTimes.^2).*(SEInitiationRates.^2) );
    R2s = getBinnedPerNucleusTrapezoidParameters(this, 'R2s', TraceType);
    y =  PlateauHeights(:,APindex, NC-8).';
    y_err = SEPlateauHeights(:, APindex,NC-8).';
    setR2s = R2s(:,APindex, NC-8).';
else
   error(['Parameter: ', parameter, ' not supported.'])
end

if ~all(isnan(y_err))
    IncludedSets = ((y./y_err > 1) & ~isnan(y) & (setR2s >= this.R2bound) & (y > 0));
else
    IncludedSets = ( ~isnan(y) & (setR2s >= this.R2bound) & (y > 0));
end

x = x(IncludedSets);
x_setpoint = x_setpoint(IncludedSets);
y = y(IncludedSets);
y_err = y_err(IncludedSets);

x = 1./(this.R*(x+273.15));
y = log(y);
y_err = abs(y_err./y);






