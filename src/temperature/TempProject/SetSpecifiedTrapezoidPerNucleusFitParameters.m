function this = SetSpecifiedTrapezoidPerNucleusFitParameters(this, p, se, ci,...
    TempIndex, NC, APindex, TraceType, ParamType)
if strcmpi(lower(ParamType), 'fits')
    if strcmpi(TraceType, 'AnaphaseAligned')
        this.PerNucleusProfileParameters.Fits.AnaphaseAligned{TempIndex, APindex, NC-8} = p;
    elseif strcmpi(TraceType, 'AnaphaseAligned3D')
        this.PerNucleusProfileParameters.Fits.AnaphaseAligned3D{TempIndex, APindex, NC-8} = p;
    elseif strcmpi(TraceType, 'Tbinned')
        this.PerNucleusProfileParameters.Fits.Tbinned{TempIndex, APindex, NC-8} = p;
    elseif strcmpi(TraceType, 'Tbinned3D')
        this.PerNucleusProfileParameters.Fits.Tbinned3D{TempIndex, APindex, NC-8} = p;
    elseif strcmpi(TraceType, 'fluo')
        this.PerNucleusProfileParameters.Fits.Unaligned{TempIndex, APindex, NC-8} = p;
    elseif strcmpi(TraceType, 'fluo3d')
        this.PerNucleusProfileParameters.Fits.Unaligned3D{TempIndex, APindex, NC-8} = p;
    end
elseif strcmp(lower(ParamType), 'meaninitiationrates')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.PerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.MeanInitiationRates.AnaphaseAlignedStdError(TempIndex, APindex, NC-8) = se;
        this.PerNucleusProfileParameters.MeanInitiationRates.AnaphaseAlignedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.PerNucleusProfileParameters.MeanInitiationRates.AnaphaseAlignedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.PerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned3D(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned3DStdError(TempIndex, APindex, NC-8) = se;
        this.PerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.PerNucleusProfileParameters.MeanInitiationRates.AnaphaseAligned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned')
        this.PerNucleusProfileParameters.MeanInitiationRates.Tbinned(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.MeanInitiationRates.TbinnedStdError(TempIndex, APindex, NC-8) = se;
        this.PerNucleusProfileParameters.MeanInitiationRates.TbinnedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.PerNucleusProfileParameters.MeanInitiationRates.TbinnedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned3d')
        this.PerNucleusProfileParameters.MeanInitiationRates.Tbinned3D(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.MeanInitiationRates.Tbinned3DStdError(TempIndex, APindex, NC-8) = se;
        this.PerNucleusProfileParameters.MeanInitiationRates.Tbinned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.PerNucleusProfileParameters.MeanInitiationRates.Tbinned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'fluo')
        this.PerNucleusProfileParameters.MeanInitiationRates.Unaligned(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.MeanInitiationRates.UnalignedStdError(TempIndex, APindex, NC-8) = se;
        this.PerNucleusProfileParameters.MeanInitiationRates.UnalignedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.PerNucleusProfileParameters.MeanInitiationRates.UnalignedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'fluo3d')
        this.PerNucleusProfileParameters.MeanInitiationRates.Unaligned3D(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.MeanInitiationRates.Unaligned3DStdError(TempIndex, APindex, NC-8) = se;
        this.PerNucleusProfileParameters.MeanInitiationRates.Unaligned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.PerNucleusProfileParameters.MeanInitiationRates.Unaligned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    end
elseif strcmp(lower(ParamType), 'timeons')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.PerNucleusProfileParameters.TimeOns.AnaphaseAligned(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.TimeOns.AnaphaseAlignedStdError(TempIndex, APindex, NC-8) = se;
        this.PerNucleusProfileParameters.TimeOns.AnaphaseAlignedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.PerNucleusProfileParameters.TimeOns.AnaphaseAlignedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.PerNucleusProfileParameters.TimeOns.AnaphaseAligned3D(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.TimeOns.AnaphaseAligned3DStdError(TempIndex, APindex, NC-8) = se;
        this.PerNucleusProfileParameters.TimeOns.AnaphaseAligned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.PerNucleusProfileParameters.TimeOns.AnaphaseAligned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned')
        this.PerNucleusProfileParameters.TimeOns.Tbinned(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.TimeOns.TbinnedStdError(TempIndex, APindex, NC-8) = se;
        this.PerNucleusProfileParameters.TimeOns.TbinnedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.PerNucleusProfileParameters.TimeOns.TbinnedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned3d')
        this.PerNucleusProfileParameters.TimeOns.Tbinned3D(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.TimeOns.Tbinned3DStdError(TempIndex, APindex, NC-8) = se;
        this.PerNucleusProfileParameters.TimeOns.Tbinned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.PerNucleusProfileParameters.TimeOns.Tbinned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'fluo')
        this.PerNucleusProfileParameters.TimeOns.Unaligned(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.TimeOns.UnalignedStdError(TempIndex, APindex, NC-8) = se;
        this.PerNucleusProfileParameters.TimeOns.UnalignedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.PerNucleusProfileParameters.TimeOns.UnalignedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'fluo3d')
        this.PerNucleusProfileParameters.TimeOns.Unaligned3D(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.TimeOns.Unaligned3DStdError(TempIndex, APindex, NC-8) = se;
        this.PerNucleusProfileParameters.TimeOns.Unaligned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.PerNucleusProfileParameters.TimeOns.Unaligned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    end
elseif strcmp(lower(ParamType), 'timeoffs')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.PerNucleusProfileParameters.TimeOffs.AnaphaseAligned(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.TimeOffs.AnaphaseAlignedStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.PerNucleusProfileParameters.TimeOffs.AnaphaseAligned3D(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.TimeOffs.AnaphaseAligned3DStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'tbinned')
        this.PerNucleusProfileParameters.TimeOffs.Tbinned(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.TimeOffs.TbinnedStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'tbinned3d')
        this.PerNucleusProfileParameters.TimeOffs.Tbinned3D(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.TimeOffs.Tbinned3DStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'fluo')
        this.PerNucleusProfileParameters.TimeOffs.Unaligned(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.TimeOffs.UnalignedStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'fluo3d')
        this.PerNucleusProfileParameters.TimeOffs.Unaligned3D(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.TimeOffs.Unaligned3DStdError(TempIndex, APindex, NC-8) = se;
    end
elseif strcmp(lower(ParamType), 'elongationtimes')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.PerNucleusProfileParameters.ElongationTimes.AnaphaseAligned(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.ElongationTimes.AnaphaseAlignedStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.PerNucleusProfileParameters.ElongationTimes.AnaphaseAligned3D(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.ElongationTimes.AnaphaseAligned3DStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'tbinned')
        this.PerNucleusProfileParameters.ElongationTimes.Tbinned(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.ElongationTimes.TbinnedStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'tbinned3d')
        this.PerNucleusProfileParameters.ElongationTimes.Tbinned3D(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.ElongationTimes.Tbinned3DStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'fluo')
        this.PerNucleusProfileParameters.ElongationTimes.Unaligned(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.ElongationTimes.UnalignedStdError(TempIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'fluo3d')
        this.PerNucleusProfileParameters.ElongationTimes.Unaligned3D(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.ElongationTimes.Unaligned3DStdError(TempIndex, APindex, NC-8) = se;
    end
elseif strcmp(lower(ParamType), 'unloadingrates')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.PerNucleusProfileParameters.UnloadingRates.AnaphaseAligned(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.UnloadingRates.AnaphaseAlignedStdError(TempIndex, APindex, NC-8) = se;
        this.PerNucleusProfileParameters.UnloadingRates.AnaphaseAlignedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.PerNucleusProfileParameters.UnloadingRates.AnaphaseAlignedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.PerNucleusProfileParameters.UnloadingRates.AnaphaseAligned3D(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.UnloadingRates.AnaphaseAligned3DStdError(TempIndex, APindex, NC-8) = se;
        this.PerNucleusProfileParameters.UnloadingRates.AnaphaseAligned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.PerNucleusProfileParameters.UnloadingRates.AnaphaseAligned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned')
        this.PerNucleusProfileParameters.UnloadingRates.Tbinned(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.UnloadingRates.TbinnedStdError(TempIndex, APindex, NC-8) = se;
        this.PerNucleusProfileParameters.UnloadingRates.TbinnedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.PerNucleusProfileParameters.UnloadingRates.TbinnedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned3d')
        this.PerNucleusProfileParameters.UnloadingRates.Tbinned3D(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.UnloadingRates.Tbinned3DStdError(TempIndex, APindex, NC-8) = se;
        this.PerNucleusProfileParameters.UnloadingRates.Tbinned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.PerNucleusProfileParameters.UnloadingRates.Tbinned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'fluo')
        this.PerNucleusProfileParameters.UnloadingRates.Unaligned(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.UnloadingRates.UnalignedStdError(TempIndex, APindex, NC-8) = se;
        this.PerNucleusProfileParameters.UnloadingRates.UnalignedCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.PerNucleusProfileParameters.UnloadingRates.UnalignedCIlow(TempIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'fluo3d')
        this.PerNucleusProfileParameters.UnloadingRates.Unaligned3D(TempIndex, APindex, NC-8) = p;
        this.PerNucleusProfileParameters.UnloadingRates.Unaligned3DStdError(TempIndex, APindex, NC-8) = se;
        this.PerNucleusProfileParameters.UnloadingRates.Unaligned3DCIhigh(TempIndex, APindex, NC-8) = ci(2);
        this.PerNucleusProfileParameters.UnloadingRates.Unaligned3DCIlow(TempIndex, APindex, NC-8) = ci(1);
    end
elseif strcmp(lower(ParamType), 'r2s')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.PerNucleusProfileParameters.MeanFitR2s.AnaphaseAligned(TempIndex, APindex, NC-8) = p;
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.PerNucleusProfileParameters.MeanFitR2s.AnaphaseAligned3D(TempIndex, APindex, NC-8) = p;
    elseif strcmp(lower(TraceType), 'tbinned')
        this.PerNucleusProfileParameters.MeanFitR2s.Tbinned(TempIndex, APindex, NC-8) = p;
    elseif strcmp(lower(TraceType), 'tbinned3d')
        this.PerNucleusProfileParameters.MeanFitR2s.Tbinned3D(TempIndex, APindex, NC-8) = p;
    elseif strcmp(lower(TraceType), 'fluo')
        this.PerNucleusProfileParameters.MeanFitR2s.Unaligned(TempIndex, APindex, NC-8) = p;
    elseif strcmp(lower(TraceType), 'fluo3d')
        this.PerNucleusProfileParameters.MeanFitR2s.Unaligned3D(TempIndex, APindex, NC-8) = p;
    end
end
end