function [params, counts] = getCycleFractionOns(this, TraceType, UseBinnedTraces, UseBinnedPerNucleusTraces)

% Next make sure TraceType is a valid argument
if strcmp(lower(TraceType), 'anaphasealigned')
    traceName = 'AnaphaseAligned';
elseif strcmp(lower(TraceType), 'anaphasealigned3d')
    traceName = 'AnaphaseAligned3D';
elseif strcmp(lower(TraceType), 'tbinned')
    traceName = 'Tbinned';
elseif strcmp(lower(TraceType), 'tbinned3d')
    traceName = 'Tbinned3D';
elseif strcmp(lower(TraceType), 'fluo') | strcmpi(lower(TraceType), 'unaligned') & ~UseBinnedTraces & ~UseBinnedPerNucleusTraces
    traceName = 'Unaligned';
elseif strcmp(lower(TraceType), 'fluo3d')| strcmpi(lower(TraceType), 'unaligned3d')& ~UseBinnedTraces & ~UseBinnedPerNucleusTraces
    traceName = 'Unaligned3D';
else
    error(['Invalid choice of TraceType: ', TraceType]);
end
NumAPbins = 41;

if UseBinnedPerNucleusTraces | UseBinnedTraces
    params = NaN(length(this.UniqueTemperatures), NumAPbins, 6);
    counts = NaN(length(this.UniqueTemperatures), NumAPbins, 6);
    AllNumNuclei = this.BinnedMeanProfiles.([traceName, 'CycleTotalNuclei']);
    AllNumOffNuclei = this.BinnedMeanProfiles.([traceName, 'CycleTotalOffNuclei']);
    FractionOn = (AllNumNuclei-AllNumOffNuclei)./AllNumNuclei;
    FractionOn(AllNumNuclei <= this.MinimumSchnitzCount) = NaN;
    CycleFractionOn = squeeze(mean(FractionOn, 1, 'omitnan'));
    CycleSchnitzCount = squeeze(max(AllNumNuclei,[], 1, 'omitnan'));
    for i = 1:length(this.UniqueTemperatures)
        params(i,:,:) = CycleFractionOn(:,:,i);
        counts(i,:,:) = CycleSchnitzCount(:,:,i);
    end
    
    
else
    params = NaN(length(this.ExperimentPrefixes), NumAPbins, 6);
    counts = NaN(length(this.ExperimentPrefixes), NumAPbins, 6);
    for i = 1:length(this.ExperimentPrefixes)
        if ~ismember(i, this.IncludedExperiments)
            continue
        end
        AllNumNuclei = this.MeanProfiles{i}.([traceName, 'CycleNumNuclei']);
        AllNumOffNuclei = this.MeanProfiles{i}.([traceName, 'CycleNumOffNuclei']);
        FractionOn = (AllNumNuclei-AllNumOffNuclei)./AllNumNuclei;
        FractionOn(AllNumNuclei <= this.MinimumSchnitzCount) = NaN;
        
        CycleFractionOn = mean(FractionOn, 1, 'omitnan');
        params(i,:,:) = CycleFractionOn;
        
        CycleSchnitzCount = squeeze(max(AllNumNuclei,[], 1, 'omitnan'));
        counts(i,:,:) = CycleSchnitzCount;
        
    end
end
