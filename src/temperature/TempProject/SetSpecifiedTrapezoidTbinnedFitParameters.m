function this = SetSpecifiedTrapezoidTbinnedFitParameters(this, p, se, ci,...
    TempIndex, NC, APindex, TraceType, ParamType)
if strcmpi(lower(ParamType), 'fits')
    if strcmpi(TraceType, 'AnaphaseAligned')
        this.BinnedProfileParameters.Fits.AnaphaseAligned{TempIndex, APindex, NC-8} = p;
    elseif strcmpi(TraceType, 'AnaphaseAligned3D')
        this.BinnedProfileParameters.Fits.AnaphaseAligned3D{TempIndex, APindex, NC-8} = p;
    elseif strcmpi(TraceType, 'Tbinned')
        this.BinnedProfileParameters.Fits.Tbinned{TempIndex, APindex, NC-8} = p;
    elseif strcmpi(TraceType, 'Tbinned3D')
        this.BinnedProfileParameters.Fits.Tbinned3D{TempIndex, APindex, NC-8} = p;
    end
elseif strcmp(lower(ParamType), 'meaninitiationrates')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.BinnedProfileParameters.MeanInitiationRates.AnaphaseAligned(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.MeanInitiationRates.AnaphaseAlignedStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedProfileParameters.MeanInitiationRates.AnaphaseAlignedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedProfileParameters.MeanInitiationRates.AnaphaseAlignedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.BinnedProfileParameters.MeanInitiationRates.AnaphaseAligned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.MeanInitiationRates.AnaphaseAligned3DStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedProfileParameters.MeanInitiationRates.AnaphaseAligned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedProfileParameters.MeanInitiationRates.AnaphaseAligned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned')
        this.BinnedProfileParameters.MeanInitiationRates.Tbinned(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.MeanInitiationRates.TbinnedStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedProfileParameters.MeanInitiationRates.TbinnedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedProfileParameters.MeanInitiationRates.TbinnedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned3d')
        this.BinnedProfileParameters.MeanInitiationRates.Tbinned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.MeanInitiationRates.Tbinned3DStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedProfileParameters.MeanInitiationRates.Tbinned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedProfileParameters.MeanInitiationRates.Tbinned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    end
elseif strcmp(lower(ParamType), 'timeons')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.BinnedProfileParameters.TimeOns.AnaphaseAligned(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.TimeOns.AnaphaseAlignedStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedProfileParameters.TimeOns.AnaphaseAlignedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedProfileParameters.TimeOns.AnaphaseAlignedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.BinnedProfileParameters.TimeOns.AnaphaseAligned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.TimeOns.AnaphaseAligned3DStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedProfileParameters.TimeOns.AnaphaseAligned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedProfileParameters.TimeOns.AnaphaseAligned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned')
        this.BinnedProfileParameters.TimeOns.Tbinned(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.TimeOns.TbinnedStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedProfileParameters.TimeOns.TbinnedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedProfileParameters.TimeOns.TbinnedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned3d')
        this.BinnedProfileParameters.TimeOns.Tbinned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.TimeOns.Tbinned3DStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedProfileParameters.TimeOns.Tbinned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedProfileParameters.TimeOns.Tbinned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    end
elseif strcmp(lower(ParamType), 'timeoffs')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.BinnedProfileParameters.TimeOffs.AnaphaseAligned(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.TimeOffs.AnaphaseAlignedStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.BinnedProfileParameters.TimeOffs.AnaphaseAligned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.TimeOffs.AnaphaseAligned3DStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'tbinned')
        this.BinnedProfileParameters.TimeOffs.Tbinned(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.TimeOffs.TbinnedStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'tbinned3d')
        this.BinnedProfileParameters.TimeOffs.Tbinned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.TimeOffs.Tbinned3DStdError(TempIndex, APindex, NC-8) = se;
    end
elseif strcmp(lower(ParamType), 'elongationtimes')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.BinnedProfileParameters.ElongationTimes.AnaphaseAligned(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.ElongationTimes.AnaphaseAlignedStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.BinnedProfileParameters.ElongationTimes.AnaphaseAligned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.ElongationTimes.AnaphaseAligned3DStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'tbinned')
        this.BinnedProfileParameters.ElongationTimes.Tbinned(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.ElongationTimes.TbinnedStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'tbinned3d')
        this.BinnedProfileParameters.ElongationTimes.Tbinned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.ElongationTimes.Tbinned3DStdError(TempIndex, APindex, NC-8) = se;
    end
elseif strcmp(lower(ParamType), 'unloadingrates')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.BinnedProfileParameters.UnloadingRates.AnaphaseAligned(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.UnloadingRates.AnaphaseAlignedStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedProfileParameters.UnloadingRates.AnaphaseAlignedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedProfileParameters.UnloadingRates.AnaphaseAlignedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.BinnedProfileParameters.UnloadingRates.AnaphaseAligned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.UnloadingRates.AnaphaseAligned3DStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedProfileParameters.UnloadingRates.AnaphaseAligned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedProfileParameters.UnloadingRates.AnaphaseAligned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned')
        this.BinnedProfileParameters.UnloadingRates.Tbinned(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.UnloadingRates.TbinnedStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedProfileParameters.UnloadingRates.TbinnedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedProfileParameters.UnloadingRates.TbinnedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned3d')
        this.BinnedProfileParameters.UnloadingRates.Tbinned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedProfileParameters.UnloadingRates.Tbinned3DStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedProfileParameters.UnloadingRates.Tbinned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedProfileParameters.UnloadingRates.Tbinned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    end
elseif strcmp(lower(ParamType), 'r2s')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.BinnedProfileParameters.MeanFitR2s.AnaphaseAligned(TempIndex, APindex, NC-8) = p;
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.BinnedProfileParameters.MeanFitR2s.AnaphaseAligned3D(TempIndex, APindex, NC-8) = p;
    elseif strcmp(lower(TraceType), 'tbinned')
        this.BinnedProfileParameters.MeanFitR2s.Tbinned(TempIndex, APindex, NC-8) = p;
    elseif strcmp(lower(TraceType), 'tbinned3d')
        this.BinnedProfileParameters.MeanFitR2s.Tbinned3D(TempIndex, APindex, NC-8) = p;
    end
end
end