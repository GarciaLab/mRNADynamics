function this = SetSpecifiedTrapezoidFitParameters(this, p, se, ci,...
    SetIndex, NC, APindex, TraceType, ParamType)
if strcmpi(lower(ParamType), 'fits')
    if strcmpi(TraceType, 'AnaphaseAligned')
        this.Fits.AnaphaseAligned{SetIndex, APindex, NC-8} = p;
    elseif strcmpi(TraceType, 'AnaphaseAligned3D')
        this.Fits.AnaphaseAligned3D{SetIndex, APindex, NC-8} = p;
    elseif strcmpi(TraceType, 'Unaligned') | strcmpi(TraceType, 'fluo')
        this.Fits.Unaligned{SetIndex, APindex, NC-8} = p;
    elseif strcmpi(TraceType, 'Unaligned3D') | strcmpi(TraceType, 'fluo3d')
        this.Fits.Unaligned3D{SetIndex, APindex, NC-8} = p;
    elseif strcmpi(TraceType, 'Tbinned')
        this.Fits.Tbinned{SetIndex, APindex, NC-8} = p;
    elseif strcmpi(TraceType, 'Tbinned3D')
        this.Fits.Tbinned3D{SetIndex, APindex, NC-8} = p;
    end
elseif strcmp(lower(ParamType), 'meaninitiationrates')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.MeanInitiationRates.AnaphaseAligned(SetIndex, APindex, NC-8) = p;
        this.MeanInitiationRates.AnaphaseAlignedStdError(SetIndex, APindex, NC-8) = se;
        this.MeanInitiationRates.AnaphaseAlignedCIhigh(SetIndex, APindex, NC-8) = ci(2);
        this.MeanInitiationRates.AnaphaseAlignedCIlow(SetIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.MeanInitiationRates.AnaphaseAligned3D(SetIndex, APindex, NC-8) = p;
        this.MeanInitiationRates.AnaphaseAligned3DStdError(SetIndex, APindex, NC-8) = se;
        this.MeanInitiationRates.AnaphaseAligned3DCIhigh(SetIndex, APindex, NC-8) = ci(2);
        this.MeanInitiationRates.AnaphaseAligned3DCIlow(SetIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned')
        this.MeanInitiationRates.Tbinned(SetIndex, APindex, NC-8) = p;
        this.MeanInitiationRates.TbinnedStdError(SetIndex, APindex, NC-8) = se;
        this.MeanInitiationRates.TbinnedCIhigh(SetIndex, APindex, NC-8) = ci(2);
        this.MeanInitiationRates.TbinnedCIlow(SetIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.MeanInitiationRates.Tbinned3D(SetIndex, APindex, NC-8) = p;
        this.MeanInitiationRates.Tbinned3DStdError(SetIndex, APindex, NC-8) = se;
        this.MeanInitiationRates.Tbinned3DCIhigh(SetIndex, APindex, NC-8) = ci(2);
        this.MeanInitiationRates.Tbinned3DCIlow(SetIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'unaligned') | strcmpi(TraceType, 'fluo')
        this.MeanInitiationRates.Unaligned(SetIndex, APindex, NC-8) = p;
        this.MeanInitiationRates.UnalignedStdError(SetIndex, APindex, NC-8) = se;
        this.MeanInitiationRates.UnalignedCIhigh(SetIndex, APindex, NC-8) = ci(2);
        this.MeanInitiationRates.UnalignedCIlow(SetIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'unaligned3d') | strcmpi(TraceType, 'fluo3d')
        this.MeanInitiationRates.Unaligned3D(SetIndex, APindex, NC-8) = p;
        this.MeanInitiationRates.Unaligned3DStdError(SetIndex, APindex, NC-8) = se;
        this.MeanInitiationRates.Unaligned3DCIhigh(SetIndex, APindex, NC-8) = ci(2);
        this.MeanInitiationRates.Unaligned7DCIlow(SetIndex, APindex, NC-8) = ci(1);
    end
elseif strcmp(lower(ParamType), 'timeons')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.TimeOns.AnaphaseAligned(SetIndex, APindex, NC-8) = p;
        this.TimeOns.AnaphaseAlignedStdError(SetIndex, APindex, NC-8) = se;
        this.TimeOns.AnaphaseAlignedCIhigh(SetIndex, APindex, NC-8) = ci(2);
        this.TimeOns.AnaphaseAlignedCIlow(SetIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.TimeOns.AnaphaseAligned3D(SetIndex, APindex, NC-8) = p;
        this.TimeOns.AnaphaseAligned3DStdError(SetIndex, APindex, NC-8) = se;
        this.TimeOns.AnaphaseAligned3DCIhigh(SetIndex, APindex, NC-8) = ci(2);
        this.TimeOns.AnaphaseAligned3DCIlow(SetIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned')
        this.TimeOns.Tbinned(SetIndex, APindex, NC-8) = p;
        this.TimeOns.TbinnedStdError(SetIndex, APindex, NC-8) = se;
        this.TimeOns.TbinnedCIhigh(SetIndex, APindex, NC-8) = ci(2);
        this.TimeOns.TbinnedCIlow(SetIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.TimeOns.Tbinned3D(SetIndex, APindex, NC-8) = p;
        this.TimeOns.Tbinned3DStdError(SetIndex, APindex, NC-8) = se;
        this.TimeOns.Tbinned3DCIhigh(SetIndex, APindex, NC-8) = ci(2);
        this.TimeOns.Tbinned3DCIlow(SetIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'unaligned') | strcmp(lower(TraceType), 'fluo') 
        this.TimeOns.Unaligned(SetIndex, APindex, NC-8) = p;
        this.TimeOns.UnalignedStdError(SetIndex, APindex, NC-8) = se;
        this.TimeOns.UnalignedCIhigh(SetIndex, APindex, NC-8) = ci(2);
        this.TimeOns.UnalignedCIlow(SetIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'unaligned3d')| strcmp(lower(TraceType), 'fluo3d') 
        this.TimeOns.Unaligned3D(SetIndex, APindex, NC-8) = p;
        this.TimeOns.Unaligned3DStdError(SetIndex, APindex, NC-8) = se;
        this.TimeOns.Unaligned3DCIhigh(SetIndex, APindex, NC-8) = ci(2);
        this.TimeOns.Unaligned7DCIlow(SetIndex, APindex, NC-8) = ci(1);
    end
elseif strcmp(lower(ParamType), 'timeoffs')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.TimeOffs.AnaphaseAligned(SetIndex, APindex, NC-8) = p;
        this.TimeOffs.AnaphaseAlignedStdError(SetIndex, APindex, NC-8) = se;
    
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.TimeOffs.AnaphaseAligned3D(SetIndex, APindex, NC-8) = p;
        this.TimeOffs.AnaphaseAligned3DStdError(SetIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'tbinned')
        this.TimeOffs.Tbinned(SetIndex, APindex, NC-8) = p;
        this.TimeOffs.TbinnedStdError(SetIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.TimeOffs.Tbinned3D(SetIndex, APindex, NC-8) = p;
        this.TimeOffs.Tbinned3DStdError(SetIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'unaligned') | strcmp(lower(TraceType), 'fluo') 
        this.TimeOffs.Unaligned(SetIndex, APindex, NC-8) = p;
        this.TimeOffs.UnalignedStdError(SetIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'unaligned3d') | strcmp(lower(TraceType), 'fluo3d') 
        this.TimeOffs.Unaligned3D(SetIndex, APindex, NC-8) = p;
        this.TimeOffs.Unaligned3DStdError(SetIndex, APindex, NC-8) = se;
    end
elseif strcmp(lower(ParamType), 'elongationtimes')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.ElongationTimes.AnaphaseAligned(SetIndex, APindex, NC-8) = p;
        this.ElongationTimes.AnaphaseAlignedStdError(SetIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.ElongationTimes.AnaphaseAligned3D(SetIndex, APindex, NC-8) = p;
        this.ElongationTimes.AnaphaseAligned3DStdError(SetIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'tbinned')
        this.ElongationTimes.Tbinned(SetIndex, APindex, NC-8) = p;
        this.ElongationTimes.TbinnedStdError(SetIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'tbinned3d')
        this.ElongationTimes.Tbinned3D(SetIndex, APindex, NC-8) = p;
        this.ElongationTimes.Tbinned3DStdError(SetIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'unaligned') | strcmp(lower(TraceType), 'fluo') 
        this.ElongationTimes.Unaligned(SetIndex, APindex, NC-8) = p;
        this.ElongationTimes.UnalignedStdError(SetIndex, APindex, NC-8) = se;
    elseif strcmp(lower(TraceType), 'unaligned3d') | strcmp(lower(TraceType), 'fluo3d') 
        this.ElongationTimes.Unaligned3D(SetIndex, APindex, NC-8) = p;
        this.ElongationTimes.Unaligned3DStdError(SetIndex, APindex, NC-8) = se;
    end
elseif strcmp(lower(ParamType), 'unloadingrates')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.UnloadingRates.AnaphaseAligned(SetIndex, APindex, NC-8) = p;
        this.UnloadingRates.AnaphaseAlignedStdError(SetIndex, APindex, NC-8) = se;
        this.UnloadingRates.AnaphaseAlignedCIhigh(SetIndex, APindex, NC-8) = ci(2);
        this.UnloadingRates.AnaphaseAlignedCIlow(SetIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.UnloadingRates.AnaphaseAligned3D(SetIndex, APindex, NC-8) = p;
        this.UnloadingRates.AnaphaseAligned3DStdError(SetIndex, APindex, NC-8) = se;
        this.UnloadingRates.AnaphaseAligned3DCIhigh(SetIndex, APindex, NC-8) = ci(2);
        this.UnloadingRates.AnaphaseAligned3DCIlow(SetIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'tbinned')
        this.UnloadingRates.Tbinned(SetIndex, APindex, NC-8) = p;
        this.UnloadingRates.TbinnedStdError(SetIndex, APindex, NC-8) = se;
        this.UnloadingRates.TbinnedCIhigh(SetIndex, APindex, NC-8) = ci(2);
        this.UnloadingRates.TbinnedCIlow(SetIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.UnloadingRates.Tbinned3D(SetIndex, APindex, NC-8) = p;
        this.UnloadingRates.Tbinned3DStdError(SetIndex, APindex, NC-8) = se;
        this.UnloadingRates.Tbinned3DCIhigh(SetIndex, APindex, NC-8) = ci(2);
        this.UnloadingRates.Tbinned3DCIlow(SetIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'unaligned') | strcmp(lower(TraceType), 'fluo') 
        this.UnloadingRates.Unaligned(SetIndex, APindex, NC-8) = p;
        this.UnloadingRates.UnalignedStdError(SetIndex, APindex, NC-8) = se;
        this.UnloadingRates.UnalignedCIhigh(SetIndex, APindex, NC-8) = ci(2);
        this.UnloadingRates.UnalignedCIlow(SetIndex, APindex, NC-8) = ci(1);
    elseif strcmp(lower(TraceType), 'unaligned3d') | strcmp(lower(TraceType), 'fluo3d') 
        this.UnloadingRates.Unaligned3D(SetIndex, APindex, NC-8) = p;
        this.UnloadingRates.Unaligned3DStdError(SetIndex, APindex, NC-8) = se;
        this.UnloadingRates.Unaligned3DCIhigh(SetIndex, APindex, NC-8) = ci(2);
        this.UnloadingRates.Unaligned7DCIlow(SetIndex, APindex, NC-8) = ci(1);
    end
elseif strcmp(lower(ParamType), 'r2s')
    if strcmp(lower(TraceType), 'anaphasealigned')
        this.MeanFitR2s.AnaphaseAligned(SetIndex, APindex, NC-8) = p;
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.MeanFitR2s.AnaphaseAligned3D(SetIndex, APindex, NC-8) = p;
    elseif strcmp(lower(TraceType), 'tbinned')
        this.MeanFitR2s.Tbinned(SetIndex, APindex, NC-8) = p;
    elseif strcmp(lower(TraceType), 'anaphasealigned3d')
        this.MeanFitR2s.Tbinned3D(SetIndex, APindex, NC-8) = p;
    elseif strcmp(lower(TraceType), 'unaligned') | strcmp(lower(TraceType), 'fluo') 
        this.MeanFitR2s.Unaligned(SetIndex, APindex, NC-8) = p;
    elseif strcmp(lower(TraceType), 'unaligned3d') | strcmp(lower(TraceType), 'fluo3d') 
        this.MeanFitR2s.Unaligned3D(SetIndex, APindex, NC-8) = p;
    end
end
end