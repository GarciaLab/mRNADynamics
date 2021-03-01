function this = AddNBFluoOffsets(this)
%this.TempFluoOffsets
%this.SetFluoOffsets

NumSets = length(this.ExperimentPrefixes);
NumAPbins = uint16(1/this.Experiments{1}.APResolution)+1;

TraceTypes = {'AnaphaseAligned', 'Tbinned'};
temperatures = flip(unique(this.Temp_sps(this.ProcessedExperiments)));
all_temperatures =  flip(unique(this.Temp_sps));

IncludeSets = ismember(1:NumSets, this.ProcessedExperiments);


MinOffsetAPbin = uint16(this.MinOffsetAP/this.Experiments{1}.APResolution)+1;
MaxOffsetAPbin = uint16(this.MaxOffsetAP/this.Experiments{1}.APResolution)+1;
OffsetBins = MinOffsetAPbin:MaxOffsetAPbin;
NumOffsetBins = length(OffsetBins);

this.TempFluoOffsets= {};
this.SetFluoOffsets = {};


for tr_index = 1:length(TraceTypes)
    TraceType = TraceTypes{tr_index};
    AllSetNumTimePoints = NaN(1, NumSets);
    for i = 1:NumSets
        if ismember(i, this.ProcessedExperiments)
            AllSetNumTimePoints(i) = uint16(GetMaxTimePointsForAllNC(this, i, TraceType));
        end
    end
    
    NumTimePoints = nanmax(AllSetNumTimePoints);
    this.TempFluoOffsets.(TraceType) = {};
    this.TempFluoOffsets.(TraceType).mean = NaN(length(all_temperatures), NumTimePoints, 6);
    this.TempFluoOffsets.(TraceType).std = NaN(length(all_temperatures), NumTimePoints, 6);
    this.TempFluoOffsets.(TraceType).count = NaN(length(all_temperatures), NumTimePoints, 6);
    
    this.SetFluoOffsets.(TraceType) = {};
    this.SetFluoOffsets.(TraceType).mean = NaN(NumSets, NumTimePoints, 6);
    this.SetFluoOffsets.(TraceType).std = NaN(NumSets, NumTimePoints, 6);
    this.SetFluoOffsets.(TraceType).count = NaN(NumSets, NumTimePoints, 6);
    this.SetFluoOffsets.(TraceType).isself = zeros(NumSets, NumTimePoints, 6, 'logical');
    
    AllSetTraces = NaN(NumSets, NumTimePoints, NumOffsetBins, 6);
    AllSetStdErrors = NaN(NumSets, NumTimePoints, NumOffsetBins, 6);
    AllSetTimes = NaN(NumSets, NumTimePoints, NumOffsetBins, 6);
    
    for SetIndex = 1:NumSets
        if ~ismember(SetIndex, this.ProcessedExperiments)
            continue
        end
        
        for NC=9:14
            for APindex = 1:NumOffsetBins
                APbin = OffsetBins(APindex);
                [FluoTrace, TraceError, TimeTrace]= GetEmbryoNBProfile(this, SetIndex, NC,APbin, TraceType);
                if ~isempty(TimeTrace)
                    TimeIndices = uint16(TimeTrace/(this.time_delta/60)+1);
                    AllSetTraces(SetIndex,TimeIndices,APindex, NC-8) = FluoTrace;
                    AllSetStdErrors(SetIndex,TimeIndices,APindex, NC-8) = TraceError;
                    AllSetTimes(SetIndex,TimeIndices,APindex, NC-8) = TimeTrace;
                end
            end
        end
    end
    
    for t_index = 1:length(temperatures)
        CurrentTemperature = temperatures(t_index);
        MatchingSetIndices = find((this.Temp_sps == CurrentTemperature) & IncludeSets);
        for NC =9:14
            nc_idx = NC-8;
            for time_index = 1:NumTimePoints
                OffsetSamples = squeeze(AllSetTraces(MatchingSetIndices, time_index,:, NC-8));
                OffsetStdErrorSamples = squeeze(AllSetStdErrors(MatchingSetIndices, time_index,:, NC-8));
                OffsetTimes = squeeze(AllSetTimes(MatchingSetIndices, time_index,:, NC-8));
                TimeVector = reshape(OffsetTimes, 1, size(OffsetTimes, 1)*size(OffsetTimes, 2));
                TimeVector = TimeVector(~isnan(TimeVector));
                if ~isempty(TimeVector)
                    if all(TimeVector == TimeVector(1))
                        OffsetVector =  reshape(OffsetSamples, 1, size(OffsetSamples, 1)*size(OffsetSamples, 2));
                        ErrorVector =  reshape(OffsetStdErrorSamples, 1, size(OffsetStdErrorSamples, 1)*size(OffsetSamples, 2));
                        ErrorVector = ErrorVector(~isnan(OffsetVector));
                        OffsetVector = OffsetVector(~isnan(OffsetVector));
                        this.TempFluoOffsets.(TraceType).mean(t_index, time_index, NC-8) = mean(OffsetVector);
                        this.TempFluoOffsets.(TraceType).std(t_index, time_index, NC-8) = sqrt(var(OffsetVector));
                        this.TempFluoOffsets.(TraceType).count(t_index, time_index, NC-8) = length(OffsetVector);
                    end
                    for MatchSetIndex = 1:length(MatchingSetIndices)
                        SetIndex = MatchingSetIndices(MatchSetIndex);
                        if SetIndex == 27
                            disp('pause point')
                        end
                            
                        SetOffsetVector = squeeze(AllSetTraces(SetIndex, time_index,:, NC-8));
                        SetErrorVector = squeeze(AllSetStdErrors(SetIndex, time_index,:, NC-8));
                        SetTimeVector = squeeze(AllSetTimes(SetIndex, time_index,:, NC-8));
                        SetTimeVector = SetTimeVector(~isnan(SetOffsetVector));
                        SetErrorVector = SetErrorVector(~isnan(SetOffsetVector));
                        SetOffsetVector = SetOffsetVector(~isnan(SetOffsetVector));
                        if ~isempty(SetTimeVector)
                            if all(SetTimeVector == SetTimeVector(1)) &  (length(SetOffsetVector(~isnan(SetOffsetVector))) >= 2)
                                this.SetFluoOffsets.(TraceType).mean(SetIndex, time_index, NC-8) = mean(SetOffsetVector);
                                this.SetFluoOffsets.(TraceType).std(SetIndex, time_index, NC-8) = sqrt(var(SetOffsetVector));
                                this.SetFluoOffsets.(TraceType).count(SetIndex, time_index, NC-8) = length(SetOffsetVector);
                                this.SetFluoOffsets.(TraceType).isself(SetIndex, time_index, NC-8) = true;
                            else
                                this.SetFluoOffsets.(TraceType).mean(SetIndex, time_index, NC-8) = this.TempFluoOffsets.(TraceType).mean(t_index, time_index, NC-8);
                                this.SetFluoOffsets.(TraceType).std(SetIndex, time_index, NC-8) = this.TempFluoOffsets.(TraceType).std(t_index, time_index, NC-8) ;
                                this.SetFluoOffsets.(TraceType).count(SetIndex, time_index, NC-8) = this.TempFluoOffsets.(TraceType).count(t_index, time_index, NC-8);
                                this.SetFluoOffsets.(TraceType).isself(SetIndex, time_index, NC-8) = false;
                            end
                        else
                            this.SetFluoOffsets.(TraceType).mean(SetIndex, time_index, NC-8) = this.TempFluoOffsets.(TraceType).mean(t_index, time_index, NC-8);
                            this.SetFluoOffsets.(TraceType).std(SetIndex, time_index, NC-8) = this.TempFluoOffsets.(TraceType).std(t_index, time_index, NC-8) ;
                            this.SetFluoOffsets.(TraceType).count(SetIndex, time_index, NC-8) = this.TempFluoOffsets.(TraceType).count(t_index, time_index, NC-8);
                            this.SetFluoOffsets.(TraceType).isself(SetIndex, time_index, NC-8) = false;
                        end
                    end
                end
                
                
                
            end
        end
    end
end


