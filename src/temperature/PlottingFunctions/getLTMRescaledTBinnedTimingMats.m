function [MeanFluoMats, StdFluoMats, NumEmbryoMats, NCTimes, MaximumNCTimes,...
    MaxFluos,MinFluos, NumFrames, MinAPs, MaxAPs] = getLTMRescaledTBinnedTimingMats(this, traceName, NC, useRescaledFluo,...
    DownsamplingRate, UseFractionOns, AltFractionOns)
%%

temperatures = this.UniqueTemperatures;
NumTemperatures = length(temperatures);
delta_t = this.time_delta/60;

MeanFluoMats = cell(1, NumTemperatures);
StdFluoMats = cell(1, NumTemperatures);
NumEmbryoMats = cell(1, NumTemperatures);
NCTimes = cell(1, NumTemperatures);
MaximumNCTimes = NaN(1, NumTemperatures);
MaxFluos = NaN(1, NumTemperatures);
MinFluos = NaN(1, NumTemperatures);
NumFrames = NaN(1, NumTemperatures);
MinAPs = NaN(1, NumTemperatures);
MaxAPs = NaN(1, NumTemperatures);

for idx = 1:NumTemperatures
    ExpNCTimes = this.BinnedMeanProfiles.([traceName, 'CycleTimes'])/60;
    IncludedRows = 1:length(ExpNCTimes);
    ExpNCTimes = ExpNCTimes(IncludedRows);
    if ~UseFractionOns & ~AltFractionOns
        ExpFluoMat = this.BinnedMeanProfiles.([traceName, 'CycleMeanTraces'])(IncludedRows,:,NC-8,idx);
        ExpStdMat =  this.BinnedMeanProfiles.([traceName, 'CycleStdErrors'])(IncludedRows,:,NC-8,idx);
    elseif ~AltFractionOns
        ExpFluoMat = this.BinnedMeanProfiles.([traceName, 'CycleTotalOnNuclei'])(IncludedRows,:,NC-8,idx)./...
            this.BinnedMeanProfiles.([traceName, 'CycleTotalNuclei'])(IncludedRows,:,NC-8,idx);
        ExpStdMat = NaN(size(ExpFluoMat));
    else
        ExpFluoMat = this.BinnedMeanProfiles.([traceName, 'CycleFractionOn'])(IncludedRows,:,NC-8,idx);
        ExpStdMat = NaN(size(ExpFluoMat));
    end
    
    ExpNumEmbryoMat = this.BinnedMeanProfiles.([traceName, 'CycleNumEmbryos'])(IncludedRows,:,NC-8,idx);
    NCLength = this.TimeScalingInfo.MeanNCDivisionInfo(idx,NC-8);
    IncludedRows = find(sum(~isnan(ExpFluoMat),2).' > 0);
    if NC == 14
        if ~isnan(this.TimeScalingInfo.PropNCDivisionInfo(idx, NC-9))
            NCTimes{idx} = ExpNCTimes(IncludedRows)/this.TimeScalingInfo.PropNCDivisionInfo(idx, NC-9);
        else
            NCTimes{idx} = this.TimingCoeffs(idx)*ExpNCTimes(IncludedRows);
        end
    elseif ~isnan(this.TimeScalingInfo.PropNCDivisionInfo(idx, NC-8))
        NCTimes{idx} = ExpNCTimes(IncludedRows)/this.TimeScalingInfo.PropNCDivisionInfo(idx, NC-8);
    else
        NCTimes{idx} = this.TimingCoeffs(idx)*ExpNCTimes(IncludedRows)/60;
    end
    
     if isempty(NCTimes{idx})
        continue
    end
    MaximumNCTimes(idx) = max(NCTimes{idx});
    MeanFluoMats{idx} = ExpFluoMat;
    StdFluoMats{idx} = ExpStdMat;
    NumEmbryoMats{idx} = ExpNumEmbryoMat;
end

SampledTimes = 0:delta_t:max(MaximumNCTimes);
TimeVectorIndices = 1:DownsamplingRate:length(SampledTimes);
if ~ismember(length(SampledTimes), TimeVectorIndices)
    TimeVectorIndices(end+1) = length(SampledTimes);
end
SampledTimes = SampledTimes(TimeVectorIndices);
%%

for idx=1:NumTemperatures
    ExpFluoMat = MeanFluoMats{idx};
    ExpStdMat =  StdFluoMats{idx};
    ExpNumEmbryoMat = NumEmbryoMats{idx};
    
    ExpFluoMat(ExpNumEmbryoMat < this.MinimumEmbryos) = NaN;
    ExpStdMat(ExpNumEmbryoMat < this.MinimumEmbryos) = NaN;
    IncludedRows = find(sum(~isnan(ExpFluoMat),2).' > 0);
    if isempty(IncludedRows)

        MeanFluoMats{idx} = [];
        StdFluoMats{idx} = [];
        NumEmbryoMats{idx} = [];
        NCTimes{idx} = [];
    else
        
        ExpFluoMat = interp1(NCTimes{idx}(IncludedRows), ExpFluoMat(IncludedRows,:), SampledTimes);
        ExpStdMat = InterpolateStdError(NCTimes{idx}(IncludedRows), ExpStdMat(IncludedRows,:), SampledTimes);
        ExpNumEmbryoMat = round(interp1(NCTimes{idx}(IncludedRows), ExpNumEmbryoMat(IncludedRows,:), SampledTimes));
        
        IncludedColumns = find(sum(~isnan(ExpFluoMat),1).' > 0);
        if ~isempty(IncludedColumns)
            MinAPs(idx) = min(IncludedColumns);
            MaxAPs(idx) = max(IncludedColumns);
        end
        
        IncludedRows = find(sum(~isnan(ExpFluoMat),2).' > 0);
        
        if isempty(IncludedRows)
            MeanFluoMats{idx} = [];
            StdFluoMats{idx} = [];
            NumEmbryoMats{idx} = [];
            NCTimes{idx} = [];
        else
            MeanFluoMats{idx} = ExpFluoMat(IncludedRows,:);
            StdFluoMats{idx} = ExpStdMat(IncludedRows,:);
            NumEmbryoMats{idx} = ExpNumEmbryoMat(IncludedRows,:);
            
            if useRescaledFluo
                MeanFluoMats{idx} = MeanFluoMats{idx}*this.FluoCoeffs(idx);
                StdFluoMats{idx} = StdFluoMats{idx}*this.FluoCoeffs(idx);
            end
            
            
            NCTimes{idx} = SampledTimes(IncludedRows);
            
            
            MaximumNCTimes(idx) = max(NCTimes{idx});
            MaxFluos(idx) = max(max(MeanFluoMats{idx}+StdFluoMats{idx}));
            MinFluos(idx) =  min(min(MeanFluoMats{idx}+StdFluoMats{idx}));
            NumFrames(idx) = length(NCTimes{idx});
        end
    end
end

 


%%

