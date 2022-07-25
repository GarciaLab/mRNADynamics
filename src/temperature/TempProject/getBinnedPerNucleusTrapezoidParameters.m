function params = getBinnedPerNucleusTrapezoidParameters(this, ParamType, TraceType, getSE)
if ~exist('getSE', 'var')
    getSE = false;
end
% First make sure ParamType is a valid argument
if strcmp(lower(ParamType), 'meaninitiationrates') 
    paramName = 'MeanInitiationRates';
elseif strcmp(lower(ParamType), 'timeons')
    paramName = 'TimeOns';
elseif strcmp(lower(ParamType), 'timeoffs')
    paramName = 'TimeOffs';
elseif strcmp(lower(ParamType), 'elongationtimes')
    paramName = 'ElongationTimes';
elseif strcmp(lower(ParamType), 'unloadingrates')
    paramName = 'UnloadingRates';
elseif strcmpi(ParamType, 'meanspotfluos')
    paramName = 'MeanSpotFluo';
elseif strcmp(lower(ParamType), 'r2s')
    paramName = 'MeanFitR2s';
    getSE = false;
else
    error(['Invalid choice of ParamType: ', ParamType])
end

% Next make sure TraceType is a valid argument
if strcmp(lower(TraceType), 'anaphasealigned')
    if ~getSE
        traceName = 'AnaphaseAligned';
    else
        traceName = 'AnaphaseAlignedStdError';
    end
elseif strcmp(lower(TraceType), 'anaphasealigned3d')
    if ~getSE
        traceName = 'AnaphaseAligned3D';
    else
        traceName = 'AnaphaseAligned3DStdError';
    end
elseif strcmp(lower(TraceType), 'tbinned')
    if ~getSE
        traceName = 'Tbinned';
    else
        traceName = 'TbinnedStdError';
    end
elseif strcmp(lower(TraceType), 'tbinned3d')
    if ~getSE
        traceName = 'Tbinned3D';
    else
        traceName = 'Tbinned3DStdError';
    end

else
    error(['Invalid choice of TraceType: ', TraceType]);
end


params = this.BinnedPerNucleusProfileParameters.(paramName).(traceName);
        