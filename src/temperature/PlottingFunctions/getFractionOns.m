function params = getFractionOns(this, TraceType, UsePerNucleusTraces, UseBinnedTraces, UseBinnedPerNucleusTraces)

% Next make sure TraceType is a valid argument
if strcmp(lower(TraceType), 'anaphasealigned')
    if ~getSE
        traceName = 'AnaphaseAligned';
    end
elseif strcmp(lower(TraceType), 'anaphasealigned3d')
    if ~getSE
        traceName = 'AnaphaseAligned3D';
    end
elseif strcmp(lower(TraceType), 'tbinned')
    if ~getSE
        traceName = 'Tbinned';
    end
elseif strcmp(lower(TraceType), 'tbinned3d')
    if ~getSE
        traceName = 'Tbinned3D';
    end
     
elseif strcmp(lower(TraceType), 'fluo') | strcmpi(lower(TraceType), 'unaligned') & ~UseBinnedTraces & ~UseBinnedPerNucleusTraces 
    if ~getSE
        traceName = 'Unaligned';
    end
elseif strcmp(lower(TraceType), 'fluo3d')| strcmpi(lower(TraceType), 'unaligned3d')& ~UseBinnedTraces & ~UseBinnedPerNucleusTraces
    if ~getSE
        traceName = 'Unaligned3D';
    end
else
    error(['Invalid choice of TraceType: ', TraceType]);
end
NumAPbins = 41;
params = NaN(length(this.ExperimentPrefixes), NumAPbins, 6);
for i = 1:length(this.ExperimentPrefixes)
    if ~ismember(i, this.IncludedExperiments)
        continue
    end
    if UseBinnedPerNucleusTraces | UseBinnedTraces
    elseif UsePerNucleusTraces
        AllNumNuclei = this.MeanProfiles{i}.([traceName, 'CycleFractionOn']);
    end
    for NC = 9:14
        for 1:NumAPbins
            
        end
    end
end

        