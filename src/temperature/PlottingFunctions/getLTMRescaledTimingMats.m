function [MeanFluoMats, StdFluoMats, NumNucMats, NCTimes, MaximumNCTimes,...
    MaxFluos, NumFrames, MinAPs, MaxAPs] = getLTMRescaledTimingMats(this, traceName, NC, useRescaledFluo)
%%
NumSets = length(this.ExperimentPrefixes);
temperatures = this.UniqueTemperatures;
delta_t = this.time_delta/60;
Temp_sp = this.Temp_sps;

MeanFluoMats = cell(1, NumSets);
StdFluoMats = cell(1, NumSets);
NumNucMats = cell(1, NumSets);
NCTimes = cell(1, NumSets);
MaximumNCTimes = NaN(1, NumSets);
MaxFluos = NaN(1, NumSets);
NumFrames = NaN(1, NumSets);
MinAPs = NaN(1, NumSets);
MaxAPs = NaN(1, NumSets);

for idx = 1:NumSets
    if ~ismember(idx, this.ProcessedExperiments)
        continue
    end
    
    ExpNCTimes = this.MeanProfiles{idx}.([traceName, 'CycleFrameTimes']){NC-8};
    IncludedRows = 1:length(ExpNCTimes);
    ExpFluoMat = squeeze(this.MeanProfiles{idx}.([traceName, 'CycleMeanTraces'])(IncludedRows,:,NC-8));
    IncludedRows = find(sum(~isnan(ExpFluoMat),2).' > 0);
    temp_match = find(round(temperatures, 2) ==  round(Temp_sp(idx), 2));
    if NC == 14
        if ~isnan(this.TimeScalingInfo.PropNCDivisionInfo(temp_match, NC-9))
            NCTimes{idx} = ExpNCTimes(IncludedRows)/60/this.TimeScalingInfo.PropNCDivisionInfo(temp_match, NC-9);
        else
            NCTimes{idx} = this.TimingCoeffs(temp_match)*ExpNCTimes(IncludedRows)/60;
        end
    elseif ~isnan(this.TimeScalingInfo.PropNCDivisionInfo(temp_match, NC-8))
        NCTimes{idx} = ExpNCTimes(IncludedRows)/60/this.TimeScalingInfo.PropNCDivisionInfo(temp_match, NC-8);
    else
        NCTimes{idx} = this.TimingCoeffs(temp_match)*ExpNCTimes(IncludedRows)/60;
    end
    
    if isempty(NCTimes{idx})
        continue
    end
    MaximumNCTimes(idx) = max(NCTimes{idx});
    
end

SampledTimes = 0:delta_t:max(MaximumNCTimes);

%%

for idx=1:NumSets
    if ~ismember(idx, this.ProcessedExperiments)
        continue
    end
    ExpNCTimes = this.MeanProfiles{idx}.([traceName, 'CycleFrameTimes']){NC-8};
    IncludedRows = 1:length(ExpNCTimes);
    ExpFluoMat = squeeze(this.MeanProfiles{idx}.([traceName, 'CycleMeanTraces'])(IncludedRows,:,NC-8));
    ExpStdMat = squeeze(this.MeanProfiles{idx}.([traceName, 'CycleTraceStdErrors'])(IncludedRows,:,NC-8));
    ExpNumNucMat = squeeze(this.MeanProfiles{idx}.([traceName, 'CycleNumOnNuclei'])(IncludedRows,:,NC-8));
    IncludedRows = find(sum(~isnan(ExpFluoMat),2).' > 0);
    if isempty(IncludedRows)

        MeanFluoMats{idx} = [];
        StdFluoMats{idx} = [];
        NumNucMats{idx} = [];
        NCTimes{idx} = [];
    else
        
        ExpFluoMat = interp1(NCTimes{idx}, ExpFluoMat(IncludedRows,:), SampledTimes);
        ExpStdMat = InterpolateStdError(NCTimes{idx}, ExpStdMat(IncludedRows,:), SampledTimes);
        ExpNumNucMat = round(interp1(NCTimes{idx}, ExpNumNucMat(IncludedRows,:), SampledTimes));
        
        IncludedColumns = find(sum(~isnan(ExpFluoMat),1).' > 0);
        if ~isempty(IncludedColumns)
            MinAPs(idx) = min(IncludedColumns);
            MaxAPs(idx) = max(IncludedColumns);
        end
        
        IncludedRows = find(sum(~isnan(ExpFluoMat),2).' > 0);
        
        if isempty(IncludedRows)
            MeanFluoMats{idx} = [];
            StdFluoMats{idx} = [];
            NumNucMats{idx} = [];
            NCTimes{idx} = [];
        else
            MeanFluoMats{idx} = ExpFluoMat(IncludedRows,:);
            StdFluoMats{idx} = ExpStdMat(IncludedRows,:);
            NumNucMats{idx} = ExpNumNucMat(IncludedRows,:);
            
            if useRescaledFluo
                temp_match = find(round(temperatures, 2) ==  round(Temp_sp(idx), 2));
                MeanFluoMats{idx} = MeanFluoMats{idx}*this.FluoCoeffs(temp_match);
                StdFluoMats{idx} = StdFluoMats{idx}*this.FluoCoeffs(temp_match);
            end
            
            
            NCTimes{idx} = SampledTimes(IncludedRows);
            
            
            MaximumNCTimes(idx) = max(NCTimes{idx});
            MaxFluos(idx) = max(max(MeanFluoMats{idx}+StdFluoMats{idx}));
            NumFrames(idx) = length(NCTimes{idx});
        end
    end
end
