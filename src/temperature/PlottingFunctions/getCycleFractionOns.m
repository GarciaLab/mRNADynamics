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
    TimeOffs = this.BinnedProfileParameters.TimeOffs.(traceName);
    TimeOns = this.BinnedProfileParameters.TimeOns.(traceName);
    
    FractionOn = (AllNumNuclei-AllNumOffNuclei)./AllNumNuclei;
    FractionOn(AllNumNuclei <= this.MinimumSchnitzCount) = NaN;
    
    
    MeanFractionOn = NaN(length(this.UniqueTemperatures), NumAPbins, 6);
    SEFractionOn = NaN(length(this.UniqueTemperatures), NumAPbins, 6);
    for APindex = 1:NumAPbins
        for nc_index = 1:6
            for temp_index = 1:5
                FractionOnVec = FractionOn(:,APindex, nc_index, temp_index).';
                TimeOnBin = TimeOns(temp_index, APindex, nc_index);
                TimeOffBin = TimeOffs(temp_index, APindex, nc_index);
                if ~isnan(TimeOnBin) & ~isnan(TimeOffBin)
                    RoundedTimeOnBin = ceil(TimeOnBin/(this.time_delta/60))+1;
                    RoundedTimeOffBin = floor(TimeOffBin/(this.time_delta/60))+1;
                    MeanFractionOn(temp_index, APindex, nc_index) = mean(FractionOnVec(RoundedTimeOnBin:RoundedTimeOffBin), 'omitnan');
                    SEFractionOn(temp_index, APindex, nc_index) = std(FractionOnVec(RoundedTimeOnBin:RoundedTimeOffBin), 'omitnan');
                    
                end
            end
        end
    end
    
    
    CycleSchnitzCount = squeeze(max(AllNumNuclei,[], 1, 'omitnan'));
    for i = 1:length(this.UniqueTemperatures)
        params(i,:,:) = MeanFractionOn(i,:,:);
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
        TimeOffs = squeeze(this.TimeOffs.(traceName)(i,:,:));
        TimeOns = squeeze(this.TimeOns.(traceName)(i,:,:));
        FractionOn = (AllNumNuclei-AllNumOffNuclei)./AllNumNuclei;
        FractionOn(AllNumNuclei <= this.MinimumSchnitzCount) = NaN;
        
        MeanFractionOn = NaN(NumAPbins, 6);
        SEFractionOn = NaN(NumAPbins, 6);
        for APindex = 1:NumAPbins
            for nc_index = 1:6
                for temp_index = 1:5
                    FractionOnVec = FractionOn(:,APindex, nc_index).';
                    TimeOnBin = TimeOns(APindex, nc_index);
                    TimeOffBin = TimeOffs(APindex, nc_index);
                    if ~isnan(TimeOnBin) & ~isnan(TimeOffBin)
                        RoundedTimeOnBin = ceil(TimeOnBin/(this.time_delta/60))+1;
                        RoundedTimeOffBin = floor(TimeOffBin/(this.time_delta/60))+1;
                        MeanFractionOn(APindex, nc_index) = mean(FractionOnVec(RoundedTimeOnBin:RoundedTimeOffBin), 'omitnan');
                        SEFractionOn( APindex, nc_index) = std(FractionOnVec(RoundedTimeOnBin:RoundedTimeOffBin), 'omitnan');
                        
                    end
                end
            end
        end
        
        params(i,:,:) = MeanFractionOn;
        
        CycleSchnitzCount = squeeze(max(AllNumNuclei,[], 1, 'omitnan'));
        counts(i,:,:) = CycleSchnitzCount;
        
    end
end
