function [pos_slope, se_pos_slope, neg_slope, se_neg_slope, time_on, se_time_on,...
    time_off, se_time_off, time_peak, se_time_peak, R2, fitresult] = ...
    getFittedTrapezoidParameters(this, SetIndex, APindex, NC, TraceType, UseBinnedParameters)
if ~exist('UseBinnedParameters', 'var')
    UseBinnedParameters = false;
end
if ~UseBinnedParameters
if strcmp(lower(TraceType), 'anaphasealigned')
    pos_slope =  this.MeanInitiationRates.AnaphaseAligned(SetIndex, APindex, NC-8);
    se_pos_slope=  this.MeanInitiationRates.AnaphaseAlignedStdError(SetIndex, APindex, NC-8);
    neg_slope = this.UnloadingRates.AnaphaseAligned(SetIndex, APindex, NC-8);
    se_neg_slope=  this.UnloadingRates.AnaphaseAlignedStdError(SetIndex, APindex, NC-8);
    time_on = this.TimeOns.AnaphaseAligned(SetIndex, APindex, NC-8);
    se_time_on = this.TimeOns.AnaphaseAlignedStdError(SetIndex, APindex, NC-8);
    time_off = this.TimeOffs.AnaphaseAligned(SetIndex, APindex, NC-8);
    se_time_off = this.TimeOffs.AnaphaseAlignedStdError(SetIndex, APindex, NC-8);
    time_peak = this.ElongationTimes.AnaphaseAligned(SetIndex, APindex, NC-8)+time_on;
    se_time_peak = sqrt(this.ElongationTimes.AnaphaseAlignedStdError(SetIndex, APindex, NC-8)^2+se_time_on^2);
    R2 = this.MeanFitR2s.AnaphaseAligned(SetIndex, APindex, NC-8);
    fitresult = this.Fits.AnaphaseAligned{SetIndex, APindex, NC-8};
elseif strcmp(lower(TraceType), 'anaphasealigned3d')
    pos_slope =  this.MeanInitiationRates.AnaphaseAligned3D(SetIndex, APindex, NC-8);
    se_pos_slope=  this.MeanInitiationRates.AnaphaseAligned3DStdError(SetIndex, APindex, NC-8);
    neg_slope = this.UnloadingRates.AnaphaseAligned3D(SetIndex, APindex, NC-8);
    se_neg_slope=  this.UnloadingRates.AnaphaseAligned3DStdError(SetIndex, APindex, NC-8);
    time_on = this.TimeOns.AnaphaseAligned3D(SetIndex, APindex, NC-8);
    se_time_on = this.TimeOns.AnaphaseAligned3DStdError(SetIndex, APindex, NC-8);
    time_off = this.TimeOffs.AnaphaseAligned3D(SetIndex, APindex, NC-8);
    se_time_off = this.TimeOffs.AnaphaseAligned3DStdError(SetIndex, APindex, NC-8);
    time_peak = this.ElongationTimes.AnaphaseAligned3D(SetIndex, APindex, NC-8)+time_on;
    se_time_peak = sqrt(this.ElongationTimes.AnaphaseAligned3DStdError(SetIndex, APindex, NC-8)^2+se_time_on^2);
    R2 = this.MeanFitR2s.AnaphaseAligned3D(SetIndex, APindex, NC-8);
    fitresult = this.Fits.AnaphaseAligned3D{SetIndex, APindex, NC-8};
elseif strcmp(lower(TraceType), 'tbinned')
    pos_slope =  this.MeanInitiationRates.Tbinned(SetIndex, APindex, NC-8);
    se_pos_slope=  this.MeanInitiationRates.TbinnedStdError(SetIndex, APindex, NC-8);
    neg_slope = this.UnloadingRates.Tbinned(SetIndex, APindex, NC-8);
    se_neg_slope=  this.UnloadingRates.TbinnedStdError(SetIndex, APindex, NC-8);
    time_on = this.TimeOns.Tbinned(SetIndex, APindex, NC-8);
    se_time_on = this.TimeOns.TbinnedStdError(SetIndex, APindex, NC-8);
    time_off = this.TimeOffs.Tbinned(SetIndex, APindex, NC-8);
    se_time_off = this.TimeOffs.TbinnedStdError(SetIndex, APindex, NC-8);
    time_peak = this.ElongationTimes.Tbinned(SetIndex, APindex, NC-8)+time_on;
    se_time_peak = sqrt(this.ElongationTimes.TbinnedStdError(SetIndex, APindex, NC-8)^2+se_time_on^2);
    R2 = this.MeanFitR2s.Tbinned(SetIndex, APindex, NC-8);
    fitresult = this.Fits.Tbinned{SetIndex, APindex, NC-8};
elseif strcmp(lower(TraceType), 'tbinned3d')
    pos_slope =  this.MeanInitiationRates.Tbinned3D(SetIndex, APindex, NC-8);
    se_pos_slope=  this.MeanInitiationRates.Tbinned3DStdError(SetIndex, APindex, NC-8);
    neg_slope = this.UnloadingRates.Tbinned3D(SetIndex, APindex, NC-8);
    se_neg_slope=  this.UnloadingRates.Tbinned3DStdError(SetIndex, APindex, NC-8);
    time_on = this.TimeOns.Tbinned3D(SetIndex, APindex, NC-8);
    se_time_on = this.TimeOns.Tbinned3DStdError(SetIndex, APindex, NC-8);
    time_off = this.TimeOffs.Tbinned3D(SetIndex, APindex, NC-8);
    se_time_off = this.TimeOffs.Tbinned3DStdError(SetIndex, APindex, NC-8);
    time_peak = this.ElongationTimes.Tbinned3D(SetIndex, APindex, NC-8)+time_on;
    se_time_peak = sqrt(this.ElongationTimes.Tbinned3DStdError(SetIndex, APindex, NC-8)^2+se_time_on^2);
    R2 = this.MeanFitR2s.Tbinned3D(SetIndex, APindex, NC-8);
    fitresult = this.Fits.Tbinned3D{SetIndex, APindex, NC-8};
elseif strcmp(lower(TraceType), 'fluo')
    pos_slope =  this.MeanInitiationRates.Unaligned(SetIndex, APindex, NC-8);
    se_pos_slope=  this.MeanInitiationRates.UnalignedStdError(SetIndex, APindex, NC-8);
    neg_slope = this.UnloadingRates.Unaligned(SetIndex, APindex, NC-8);
    se_neg_slope=  this.UnloadingRates.UnalignedStdError(SetIndex, APindex, NC-8);
    time_on = this.TimeOns.Unaligned(SetIndex, APindex, NC-8);
    se_time_on = this.TimeOns.UnalignedStdError(SetIndex, APindex, NC-8);
    time_off = this.TimeOffs.Unaligned(SetIndex, APindex, NC-8);
    se_time_off = this.TimeOffs.UnalignedStdError(SetIndex, APindex, NC-8);
    time_peak = this.ElongationTimes.Unaligned(SetIndex, APindex, NC-8)+time_on;
    se_time_peak = sqrt(this.ElongationTimes.UnalignedStdError(SetIndex, APindex, NC-8)^2+se_time_on^2);
    R2 = this.MeanFitR2s.Unaligned(SetIndex, APindex, NC-8);
    fitresult = this.Fits.Unaligned{SetIndex, APindex, NC-8};
elseif strcmp(lower(TraceType), 'fluo3d')
    pos_slope =  this.MeanInitiationRates.Unaligned3D(SetIndex, APindex, NC-8);
    se_pos_slope=  this.MeanInitiationRates.Unaligned3DStdError(SetIndex, APindex, NC-8);
    neg_slope = this.UnloadingRates.Unaligned3D(SetIndex, APindex, NC-8);
    se_neg_slope=  this.UnloadingRates.Unaligned3DStdError(SetIndex, APindex, NC-8);
    time_on = this.TimeOns.Unaligned3D(SetIndex, APindex, NC-8);
    se_time_on = this.TimeOns.Unaligned3DStdError(SetIndex, APindex, NC-8);
    time_off = this.TimeOffs.Unaligned3D(SetIndex, APindex, NC-8);
    se_time_off = this.TimeOffs.Unaligned3DStdError(SetIndex, APindex, NC-8);
    time_peak = this.ElongationTimes.Unaligned3D(SetIndex, APindex, NC-8)+time_on;
    se_time_peak = sqrt(this.ElongationTimes.Unaligned3DStdError(SetIndex, APindex, NC-8)^2+se_time_on^2);
    R2 = this.MeanFitR2s.Unaligned3D(SetIndex, APindex, NC-8);
    fitresult = this.Fits.Unaligned3D{SetIndex, APindex, NC-8};
end
else
    if strcmp(lower(TraceType), 'anaphasealigned')
        traceName = 'AnaphaseAligned';
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        traceName = 'AnaphaseAligned3D';
    elseif strcmp(lower(TraceType), 'tbinned')
        traceName = 'Tbinned';
    elseif strcmp(lower(TraceType), 'tbinned3d')
        traceName = 'Tbinned3D';
    elseif strcmp(lower(TraceType), 'fluo')
        traceName = 'Unaligned';
    elseif strcmp(lower(TraceType), 'fluo3d')
        traceName = 'Unaligned3D';
    end
    pos_slope =  this.BinnedProfileParameters.MeanInitiationRates.(traceName)(SetIndex, APindex, NC-8);
    se_pos_slope=  this.BinnedProfileParameters.MeanInitiationRates.([traceName, 'StdError'])(SetIndex, APindex, NC-8);
    neg_slope = this.BinnedProfileParameters.UnloadingRates.(traceName)(SetIndex, APindex, NC-8);
    se_neg_slope=  this.BinnedProfileParameters.UnloadingRates.([traceName, 'StdError'])(SetIndex, APindex, NC-8);
    time_on = this.BinnedProfileParameters.TimeOns.(traceName)(SetIndex, APindex, NC-8);
    se_time_on = this.BinnedProfileParameters.TimeOns.([traceName, 'StdError'])(SetIndex, APindex, NC-8);
    time_off = this.BinnedProfileParameters.TimeOffs.(traceName)(SetIndex, APindex, NC-8);
    se_time_off = this.BinnedProfileParameters.TimeOffs.([traceName, 'StdError'])(SetIndex, APindex, NC-8);
    time_peak = this.BinnedProfileParameters.ElongationTimes.(traceName)(SetIndex, APindex, NC-8)+time_on;
    se_time_peak = sqrt(this.BinnedProfileParameters.ElongationTimes.([traceName, 'StdError'])(SetIndex, APindex, NC-8)^2+se_time_on^2);
    R2 = this.BinnedProfileParameters.MeanFitR2s.(traceName)(SetIndex, APindex, NC-8);
    fitresult = this.BinnedProfileParameters.Fits.(traceName){SetIndex, APindex, NC-8};
    end
%disp('Check point')
