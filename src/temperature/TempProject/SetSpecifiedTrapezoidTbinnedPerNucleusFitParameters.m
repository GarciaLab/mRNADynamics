function this = SetSpecifiedTrapezoidTbinnedPerNucleusFitParameters(this, p, se, ci,...
    TempIndex, NC, APindex, TraceType, ParamType)
if strcmpi(lower(ParamType), 'fits')
    if strcmpi(TraceType, 'AnaphaseAligned')
        this.BinnedPerNucleusProfileParameters.Fits.AnaphaseAligned{TempIndex, APindex, NC-8} = p;
    elseif strcmpi(TraceType, 'AnaphaseAligned3D')
        this.BinnedPerNucleusProfileParameters.Fits.AnaphaseAligned3D{TempIndex, APindex, NC-8} = p;
    elseif strcmpi(TraceType, 'Tbinned')
        this.BinnedPerNucleusProfileParameters.Fits.Tbinned{TempIndex, APindex, NC-8} = p;
    elseif strcmpi(TraceType, 'Tbinned3D')
        this.BinnedPerNucleusProfileParameters.Fits.Tbinned3D{TempIndex, APindex, NC-8} = p;
    end
elseif strcmp(lower(ParamType), 'meaninitiationrates')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.BinnedPerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.MeanInitiationRates.AnaphaseAlignedStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedPerNucleusProfileParameters.MeanInitiationRates.AnaphaseAlignedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedPerNucleusProfileParameters.MeanInitiationRates.AnaphaseAlignedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.BinnedPerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned3DStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedPerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedPerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned')
        this.BinnedPerNucleusProfileParameters.MeanInitiationRates.Tbinned(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.MeanInitiationRates.TbinnedStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedPerNucleusProfileParameters.MeanInitiationRates.TbinnedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedPerNucleusProfileParameters.MeanInitiationRates.TbinnedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned3d')
        this.BinnedPerNucleusProfileParameters.MeanInitiationRates.Tbinned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.MeanInitiationRates.Tbinned3DStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedPerNucleusProfileParameters.MeanInitiationRates.Tbinned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedPerNucleusProfileParameters.MeanInitiationRates.Tbinned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    end
elseif strcmp(lower(ParamType), 'timeons')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.BinnedPerNucleusProfileParameters.TimeOns.AnaphaseAligned(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.TimeOns.AnaphaseAlignedStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedPerNucleusProfileParameters.TimeOns.AnaphaseAlignedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedPerNucleusProfileParameters.TimeOns.AnaphaseAlignedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.BinnedPerNucleusProfileParameters.TimeOns.AnaphaseAligned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.TimeOns.AnaphaseAligned3DStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedPerNucleusProfileParameters.TimeOns.AnaphaseAligned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedPerNucleusProfileParameters.TimeOns.AnaphaseAligned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned')
        this.BinnedPerNucleusProfileParameters.TimeOns.Tbinned(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.TimeOns.TbinnedStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedPerNucleusProfileParameters.TimeOns.TbinnedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedPerNucleusProfileParameters.TimeOns.TbinnedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned3d')
        this.BinnedPerNucleusProfileParameters.TimeOns.Tbinned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.TimeOns.Tbinned3DStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedPerNucleusProfileParameters.TimeOns.Tbinned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedPerNucleusProfileParameters.TimeOns.Tbinned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    end
elseif strcmp(lower(ParamType), 'timeoffs')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.BinnedPerNucleusProfileParameters.TimeOffs.AnaphaseAligned(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.TimeOffs.AnaphaseAlignedStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.BinnedPerNucleusProfileParameters.TimeOffs.AnaphaseAligned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.TimeOffs.AnaphaseAligned3DStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'tbinned')
        this.BinnedPerNucleusProfileParameters.TimeOffs.Tbinned(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.TimeOffs.TbinnedStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'tbinned3d')
        this.BinnedPerNucleusProfileParameters.TimeOffs.Tbinned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.TimeOffs.Tbinned3DStdError(TempIndex, APindex, NC-8) = se;
    end
elseif strcmp(lower(ParamType), 'elongationtimes')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.BinnedPerNucleusProfileParameters.ElongationTimes.AnaphaseAligned(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.ElongationTimes.AnaphaseAlignedStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.BinnedPerNucleusProfileParameters.ElongationTimes.AnaphaseAligned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.ElongationTimes.AnaphaseAligned3DStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'tbinned')
        this.BinnedPerNucleusProfileParameters.ElongationTimes.Tbinned(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.ElongationTimes.TbinnedStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'tbinned3d')
        this.BinnedPerNucleusProfileParameters.ElongationTimes.Tbinned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.ElongationTimes.Tbinned3DStdError(TempIndex, APindex, NC-8) = se;
    end
elseif strcmp(lower(ParamType), 'unloadingrates')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.BinnedPerNucleusProfileParameters.UnloadingRates.AnaphaseAligned(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.UnloadingRates.AnaphaseAlignedStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedPerNucleusProfileParameters.UnloadingRates.AnaphaseAlignedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedPerNucleusProfileParameters.UnloadingRates.AnaphaseAlignedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.BinnedPerNucleusProfileParameters.UnloadingRates.AnaphaseAligned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.UnloadingRates.AnaphaseAligned3DStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedPerNucleusProfileParameters.UnloadingRates.AnaphaseAligned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedPerNucleusProfileParameters.UnloadingRates.AnaphaseAligned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned')
        this.BinnedPerNucleusProfileParameters.UnloadingRates.Tbinned(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.UnloadingRates.TbinnedStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedPerNucleusProfileParameters.UnloadingRates.TbinnedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedPerNucleusProfileParameters.UnloadingRates.TbinnedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned3d')
        this.BinnedPerNucleusProfileParameters.UnloadingRates.Tbinned3D(TempIndex, APindex, NC-8) = p;
        this.BinnedPerNucleusProfileParameters.UnloadingRates.Tbinned3DStdError(TempIndex, APindex, NC-8) = se;
        this.BinnedPerNucleusProfileParameters.UnloadingRates.Tbinned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.BinnedPerNucleusProfileParameters.UnloadingRates.Tbinned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    end
elseif strcmp(lower(ParamType), 'r2s')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.BinnedPerNucleusProfileParameters.MeanFitR2s.AnaphaseAligned(TempIndex, APindex, NC-8) = p;
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.BinnedPerNucleusProfileParameters.MeanFitR2s.AnaphaseAligned3D(TempIndex, APindex, NC-8) = p;
    elseif strcmp(lower(TraceType), 'tbinned')
        this.BinnedPerNucleusProfileParameters.MeanFitR2s.Tbinned(TempIndex, APindex, NC-8) = p;
    elseif strcmp(lower(TraceType), 'tbinned3d')
        this.BinnedPerNucleusProfileParameters.MeanFitR2s.Tbinned3D(TempIndex, APindex, NC-8) = p;
    end
end
end