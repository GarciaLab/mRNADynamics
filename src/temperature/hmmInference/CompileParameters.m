function CompiledParameters = CompileParameters(ResultsPaths,includeRescaling, PrescaledBins, TempParamPath, FluoParamPath)
%%
if ~exist('includeRescaling', 'var')
    includeRescaling = true;
end
if ~exist('TempParamPath', 'var')
    TempParamPath = 'S:/Gabriella\Dropbox\TemperatureParameters\hbBAC-MS2\';
    TempParamPath = 'S:/Gabriella\Dropbox\TemperatureParameters\hbBAC-MS2-V6\';
    warning('Using general temperature parameters for hbBAC-MS2.');
end

if ~exist('FluoParamPath', 'var')
    FluoParamPath = 'S:/Gabriella\Dropbox\TemperatureParameters\GFP\';
    warning('Using fluo temperature parameters for GFP.');
end
if ~exist('PrescaledBins', 'var')
    PrescaledBins = false;
end

%%

if ~PrescaledBins
    
    
    APResolution = 0.025;
    APbins = 0:APResolution:1;
    NumAPbins = length(APbins);
    
    TimeBinWidth = 15;
    TimeBinSep = 5;
    RefTemperature = 25.0;
    ScalingDeltaT = 2.5; % min
    
    NumSets = length(ResultsPaths);
    
    
    CompiledParameters = {};
    CompiledParameters.ResultsPaths = ResultsPaths;
    CompiledParameters.FigurePaths = ResultsPaths;
    CompiledParameters.ReporterLabels = cell(1, NumSets);
    
    CompiledParameters.APVector = APbins;
    
    
    
    EnrichmentPath = getProcessedEnrichmentFolder();
    
    
    
    %%
    
    prefix_regex = '(?<projectDirInfo>[A-Za-z0-9\-_\\]+)\\cpHMM_results\\compiledResults_w(?<nSteps>[0-9]+)_K(?<nStates>[0-9]+)_p0_ap(?<nAPbins>[0-9]+)_t(?<nTimebins>[0-9]+)_f(?<FluoDim>[0-9]+)D(?<TimeRes>[a-z0-9._]+)';
    cycle_regex = '(?<projectInfo>[A-Za-z\-0-9_]+)\\NC(?<cycle>[0-9]+)';
    timeres_regex = '_dt(?<dt>[0-9]+).mat';
    projectInfo_regex = 'hbBAC-MS2-(?<temperature>[0-9_]+)C';
    
    
    %%
    for ResultIndex=1:NumSets
        results_path_split = regexp(ResultsPaths{ResultIndex}, prefix_regex, 'names');
        if ~isempty(results_path_split)
            [fp1, fp2, fp3] = fileparts(ResultsPaths{ResultIndex});
            figpath = [EnrichmentPath,  filesep, fp1, filesep, 'figures', filesep, fp2, filesep];
            if ~exist(figpath, 'dir')
                mkdir(figpath);
            end
            CompiledParameters.FigurePaths{ResultIndex} = figpath;
            
            cycle_split = regexp(results_path_split.projectDirInfo, cycle_regex, 'names');
            if ~isempty(cycle_split)
                CompiledParameters.NC(ResultIndex) = str2num(cycle_split.cycle);
                proj_split_v1 = regexp(cycle_split.projectInfo, projectInfo_regex, 'names');
                PotentialProjName = cycle_split.projectInfo;
            else
                CompiledParameters.NC(ResultIndex) = 14;
                proj_split_v1 = regexp(results_path_split.projectDirInfo, projectInfo_regex, 'names');
                PotentialProjName = results_path_split.projectDirInfo;
            end
            
            if ~isempty(proj_split_v1)
                CompiledParameters.ReporterLabels{ResultIndex} = 'hbBAC-MS2';
                CompiledParameters.SetTemperatures(ResultIndex) = str2num(strrep(proj_split_v1.temperature, '_', '.'));
            else
                CompiledParameters.ReporterLabels{ResultIndex} = PotentialProjName;
                CompiledParameters.SetTemperatures(ResultIndex) = 25;
            end
            
            
            
            CompiledParameters.nSteps(ResultIndex) = str2num(results_path_split.nSteps);
            CompiledParameters.nStates(ResultIndex) = str2num(results_path_split.nStates);
            CompiledParameters.nAPbins(ResultIndex) = str2num(results_path_split.nAPbins);
            CompiledParameters.nTimebins(ResultIndex) = str2num(results_path_split.nTimebins);
            CompiledParameters.nTimebins(ResultIndex) = str2num(results_path_split.nTimebins);
            CompiledParameters.FluoDims(ResultIndex) = str2num(results_path_split.FluoDim);
            timeres_split = regexp(results_path_split.TimeRes, timeres_regex, 'names');
            if ~isempty(timeres_split)
                CompiledParameters.dt(ResultIndex) = str2num(timeres_split.dt);
            else
                datadir = strsplit(ResultsPaths{ResultIndex}, '\\cpHMM_results');
                load([EnrichmentPath datadir{1} filesep 'spot_struct.mat']);
                for i = 1:length(spot_struct)
                    if length(spot_struct(i).timeInterp) > 1
                        if all(~isnan(spot_struct(i).timeInterp))
                            CompiledParameters.dt(ResultIndex) = spot_struct(i).timeInterp(2)-spot_struct(i).timeInterp(1);
                            break
                        end
                    end
                    
                end
            end
            
            CompiledParameters.ElongationTimes(ResultIndex) = CompiledParameters.nSteps(ResultIndex)*CompiledParameters.dt(ResultIndex)/60;
            
            
        else
            error('Inference Info path has unexpected format.')
        end
    end
    %%
    NumTimeBins = max(CompiledParameters.nTimebins);
    CompiledParameters.TimeLimits = NaN(NumTimeBins, 2);
    CompiledParameters.TimeVector = NaN(1, NumTimeBins);
    for i = 1:NumTimeBins
        CompiledParameters.TimeLimits(i, 1) = (i-1)*5;
        CompiledParameters.TimeLimits(i, 2) = (i-1)*5+15;
        CompiledParameters.TimeVector(i) = mean(CompiledParameters.TimeLimits(i, :));
    end
    
    CompiledParameters.APRange =NaN(1, 2);
    CompiledParameters.SetTemperatures = NaN(1,NumSets);
    CompiledParameters.NC = NaN(1,NumSets);
    CompiledParameters.nSteps = NaN(1,NumSets);
    CompiledParameters.nStates = NaN(1,NumSets);
    CompiledParameters.nAPbins = NaN(1,NumSets);
    CompiledParameters.nTimebins = NaN(1,NumSets);
    CompiledParameters.FluoDims = NaN(1, NumSets);
    CompiledParameters.dt = NaN(1,NumSets);
    CompiledParameters.TimeAxes = cell(1, NumSets);
    CompiledParameters.ElongationTimes = NaN(1,NumSets);
    CompiledParameters.Durations = NaN(NumSets, NumTimeBins, NumAPbins);
    CompiledParameters.Frequencies = NaN(NumSets, NumTimeBins, NumAPbins);
    CompiledParameters.InitiationRates = NaN(NumSets, NumTimeBins, NumAPbins);
    CompiledParameters.DurationsStdErr = NaN(NumSets, NumTimeBins, NumAPbins);
    CompiledParameters.FrequenciesStdErr = NaN(NumSets, NumTimeBins, NumAPbins);
    CompiledParameters.InitiationRatesStdErr = NaN(NumSets, NumTimeBins, NumAPbins);
    
    %%
    
    for ResultIndex=1:NumSets
        results_path_split = regexp(ResultsPaths{ResultIndex}, prefix_regex, 'names');
        if ~isempty(results_path_split)
            [fp1, fp2, fp3] = fileparts(ResultsPaths{ResultIndex});
            figpath = [EnrichmentPath,  filesep, fp1, filesep, 'figures', filesep, fp2, filesep];
            if ~exist(figpath, 'dir')
                mkdir(figpath);
            end
            CompiledParameters.FigurePaths{ResultIndex} = figpath;
            
            cycle_split = regexp(results_path_split.projectDirInfo, cycle_regex, 'names');
            if ~isempty(cycle_split)
                CompiledParameters.NC(ResultIndex) = str2num(cycle_split.cycle);
                proj_split_v1 = regexp(cycle_split.projectInfo, projectInfo_regex, 'names');
                PotentialProjName = cycle_split.projectInfo;
            else
                CompiledParameters.NC(ResultIndex) = 14;
                proj_split_v1 = regexp(results_path_split.projectDirInfo, projectInfo_regex, 'names');
                PotentialProjName = results_path_split.projectDirInfo;
            end
            
            if ~isempty(proj_split_v1)
                CompiledParameters.ReporterLabels{ResultIndex} = 'hbBAC-MS2';
                CompiledParameters.SetTemperatures(ResultIndex) = str2num(strrep(proj_split_v1.temperature, '_', '.'));
            else
                CompiledParameters.ReporterLabels{ResultIndex} = PotentialProjName;
                CompiledParameters.SetTemperatures(ResultIndex) = 25;
            end
            
            
            
            CompiledParameters.nSteps(ResultIndex) = str2num(results_path_split.nSteps);
            CompiledParameters.nStates(ResultIndex) = str2num(results_path_split.nStates);
            CompiledParameters.nAPbins(ResultIndex) = str2num(results_path_split.nAPbins);
            CompiledParameters.nTimebins(ResultIndex) = str2num(results_path_split.nTimebins);
            CompiledParameters.nTimebins(ResultIndex) = str2num(results_path_split.nTimebins);
            CompiledParameters.FluoDims(ResultIndex) = str2num(results_path_split.FluoDim);
            timeres_split = regexp(results_path_split.TimeRes, timeres_regex, 'names');
            if ~isempty(timeres_split)
                CompiledParameters.dt(ResultIndex) = str2num(timeres_split.dt);
            else
                datadir = strsplit(ResultsPaths{ResultIndex}, '\\cpHMM_results');
                load([EnrichmentPath datadir{1} filesep 'spot_struct.mat']);
                for i = 1:length(spot_struct)
                    if length(spot_struct(i).timeInterp) > 1
                        if all(~isnan(spot_struct(i).timeInterp))
                            CompiledParameters.dt(ResultIndex) = spot_struct(i).timeInterp(2)-spot_struct(i).timeInterp(1);
                            break
                        end
                    end
                    
                end
            end
            
            CompiledParameters.ElongationTimes(ResultIndex) = CompiledParameters.nSteps(ResultIndex)*CompiledParameters.dt(ResultIndex)/60;
            
            
        else
            error('Inference Info path has unexpected format.')
        end
        
        load([EnrichmentPath, ResultsPaths{ResultIndex}]);
        
        time_group_vec = compiledResults.timeGroupVec;
        time_group_index = unique(time_group_vec);
        ap_group_vec = compiledResults.apGroupVec;
        ap_group_index = unique(ap_group_vec);
        ap_axis = compiledResults.apBins(1:end-1) + diff(compiledResults.apBins);
        time_axis = [];
        for t = time_group_index
            time_axis(end+1) = mean(compiledResults.timeBins{t})/60;
        end
        CompiledParameters.TimeAxes{ResultIndex} = time_axis;
        
        for t = time_group_index
            time_filter = time_group_vec==t;
            ap_ids = ap_group_vec(time_filter);
            ap_bins = ap_axis(ap_ids)/100;
            [ap_bins, sort_order] = sort(ap_bins);
            AP_index = ismember(round(APbins, 3), round(ap_bins, 3));
            
            dur_vec_mean = compiledResults.dur_vec_mean(time_filter);
            dur_vec_mean = dur_vec_mean(sort_order);
            dur_vec_ste = compiledResults.dur_vec_ste(time_filter);
            dur_vec_ste = dur_vec_ste(sort_order);
            
            CompiledParameters.Durations(ResultIndex, t, AP_index) = dur_vec_mean;
            CompiledParameters.DurationsStdErr(ResultIndex, t, AP_index) = dur_vec_ste;
            
            freq_vec_mean = compiledResults.freq_vec_mean(time_filter);
            freq_vec_mean = freq_vec_mean(sort_order);
            freq_vec_ste = compiledResults.freq_vec_ste(time_filter);
            freq_vec_ste = freq_vec_ste(sort_order);
            
            CompiledParameters.Frequencies(ResultIndex, t, AP_index) = freq_vec_mean;
            CompiledParameters.FrequenciesStdErr(ResultIndex, t, AP_index) = freq_vec_ste;
            
            
            init_vec_mean = compiledResults.init_vec_mean(time_filter);
            init_vec_mean = init_vec_mean(sort_order);
            init_vec_ste = compiledResults.init_vec_ste(time_filter);
            init_vec_ste = init_vec_ste(sort_order);
            
            
            CompiledParameters.InitiationRates(ResultIndex, t, AP_index) = init_vec_mean;
            CompiledParameters.InitiationRatesStdErr(ResultIndex, t, AP_index) = init_vec_ste;
        end
        
    end
    
    
    %%
    if includeRescaling
        CompiledParameters.UniqueTemperatures = fliplr(unique(CompiledParameters.SetTemperatures));
        
        if ~isempty(TempParamPath)
            if isfile([TempParamPath filesep 'autocorrElongationTimes.mat'])
                load([TempParamPath filesep 'autocorrElongationTimes.mat']);
                ElongTemps =  ElongationTimeInfo.Temperatures;
                CompiledParameters.AutocorrElongationTimes = NaN(1, length(CompiledParameters.UniqueTemperatures));
                for i = 1:length(ElongTemps)
                    match_index = find(round(CompiledParameters.UniqueTemperatures, 1) == round(ElongTemps(i), 1));
                    if ~isempty(match_index)
                        CompiledParameters.AutocorrElongationTimes(match_index) = ElongationTimeInfo.ElongationTimes(i);
                    end
                end
            else
                CompiledParameters.AutocorrElongationTimes = NaN(1, length(CompiledParameters.UniqueTemperatures));
                 CompiledParameters.AutocorrElongationTimes= [175 210 245 315 405];
            end
            
            if isfile([TempParamPath filesep 'DevTimeCoefficients.mat'])
                load([TempParamPath filesep 'DevTimeCoefficients.mat']);
                DevTimeTemps =  TimingInfo.Temperatures;
                CompiledParameters.TimingCoeffs = NaN(1, length(CompiledParameters.UniqueTemperatures));
                for i = 1:length(DevTimeTemps)
                    match_index = find(round(CompiledParameters.UniqueTemperatures, 1) == round(DevTimeTemps(i), 1));
                    if ~isempty(match_index)
                        CompiledParameters.TimingCoeffs(match_index) =TimingInfo.MeanNCDivisionInfo(2, 5)/TimingInfo.MeanNCDivisionInfo(i, 5);
                    end
                end
            else
                CompiledParameters.TimingCoeffs = NaN(1, length(CompiledParameters.UniqueTemperatures));
            end
            
            
            
            
        else
            CompiledParameters.AutocorrElongationTimes = NaN(1, length(CompiledParameters.UniqueTemperatures));
            CompiledParameters.AutocorrElongationTimes= [175 210 245 315 405];
            CompiledParameters.TimingCoeffs = NaN(1, length(CompiledParameters.UniqueTemperatures));
            
        end
        
        
        if ~isempty(FluoParamPath)
            if isfile([FluoParamPath filesep 'FluoCoefficients.mat'])
                load([FluoParamPath filesep 'FluoCoefficients.mat']);
                FluoTemps =  FluoInfo.Temperatures;
                CompiledParameters.FluoCoeffs = NaN(1, length(CompiledParameters.UniqueTemperatures));
                for i = 1:length(FluoTemps)
                    match_index = find(round(CompiledParameters.UniqueTemperatures, 1) == round(FluoTemps(i), 1));
                    if ~isempty(match_index)
                        CompiledParameters.FluoCoeffs(match_index) = FluoInfo.FluoCoeffs(i);
                    end
                end
            else
                CompiledParameters.FluoCoeffs = NaN(1, length(CompiledParameters.UniqueTemperatures));
            end
            
        else
            CompiledParameters.FluoCoeffs = NaN(1, length(CompiledParameters.UniqueTemperatures));
            
        end
        
        %%
        
        RefIndex = find(round(CompiledParameters.UniqueTemperatures, 1) == round(RefTemperature, 1));
        
        AllObservedTimes = [];
        dts = [];
        min_time_bins = [];
        max_time_bins = [];
        for i = 1:length(CompiledParameters.TimeAxes)
            t_index = find(round(CompiledParameters.UniqueTemperatures, 1) == round(CompiledParameters.SetTemperatures(i), 1));
            t_coeff = CompiledParameters.TimingCoeffs(t_index);
            if CompiledParameters.TimeAxes{i} ~= Inf
                AllObservedTimes = [AllObservedTimes CompiledParameters.TimeAxes{i}*t_coeff];
                dts(end+1) = AllObservedTimes(end)-AllObservedTimes(end-1);
                min_time_bins(end+1) = min(CompiledParameters.TimeAxes{i}*t_coeff);
                max_time_bins(end+1) = max(CompiledParameters.TimeAxes{i}*t_coeff);
            end
            
        end
        
        if ~isempty(AllObservedTimes)
            
            Shared_min_value = min(ceil(min_time_bins/ScalingDeltaT));
            Shared_max_value = max(floor(max_time_bins/ScalingDeltaT));
            
            ReferenceTimeVector = (Shared_min_value:Shared_max_value)*ScalingDeltaT;
            if round(min_time_bins(RefIndex), 3) == round(Shared_min_value*ScalingDeltaT,3)
                HalfRefTimeVector = (Shared_min_value:2:Shared_max_value)*ScalingDeltaT;
            else
                HalfRefTimeVector = ((Shared_min_value+1):2:Shared_max_value)*ScalingDeltaT;
            end
            
            if HalfRefTimeVector(end) > ReferenceTimeVector(end)
                HalfRefTimeVector = HalfRefTimeVector(1:end-1);
            end
            
            
            IncludedHalfSampledBins = ismember(ReferenceTimeVector, HalfRefTimeVector);
            NumScaledBins = length(ReferenceTimeVector);
            BinnedTimes = true;
        else
            NumScaledBins = 1;
            ReferenceTimeVector = Inf;
            BinnedTimes = false;
        end
        
        %%
        CompiledParameters.ScaledTimeVector = ReferenceTimeVector;
        CompiledParameters.ScaledDurations = NaN(NumSets, NumScaledBins, NumAPbins);
        CompiledParameters.ScaledFrequencies = NaN(NumSets, NumScaledBins, NumAPbins);
        CompiledParameters.ScaledInitiationRates = NaN(NumSets, NumScaledBins, NumAPbins);
        CompiledParameters.ScaledDurationsStdErr = NaN(NumSets, NumScaledBins, NumAPbins);
        CompiledParameters.ScaledFrequenciesStdErr = NaN(NumSets, NumScaledBins, NumAPbins);
        CompiledParameters.ScaledInitiationRatesStdErr = NaN(NumSets, NumScaledBins, NumAPbins);
        
        CompiledParameters.ScaledDurationsV2 = NaN(NumSets, NumScaledBins, NumAPbins);
        CompiledParameters.ScaledFrequenciesV2 = NaN(NumSets, NumScaledBins, NumAPbins);
        CompiledParameters.ScaledInitiationRatesV2 = NaN(NumSets, NumScaledBins, NumAPbins);
        CompiledParameters.ScaledDurationsStdErrV2 = NaN(NumSets, NumScaledBins, NumAPbins);
        CompiledParameters.ScaledFrequenciesStdErrV2 = NaN(NumSets, NumScaledBins, NumAPbins);
        CompiledParameters.ScaledInitiationRatesStdErrV2 = NaN(NumSets, NumScaledBins, NumAPbins);
        
        
        %%
        if BinnedTimes
            for SetIndex = 1:NumSets
                t_index = find(round(CompiledParameters.UniqueTemperatures, 1) == round(CompiledParameters.SetTemperatures(SetIndex), 1));
                t_coeff = CompiledParameters.TimingCoeffs(t_index);
                SampleTimeVector = CompiledParameters.TimeVector*t_coeff;
                for APbin = 1:NumAPbins
                    if all(isnan(CompiledParameters.InitiationRates(SetIndex,:,APbin)))
                        continue
                    end
                    
                    SampleVector = interp1(SampleTimeVector, squeeze(CompiledParameters.InitiationRates(SetIndex,:,APbin)), ReferenceTimeVector);
                    %if dts(SetIndex)
                    %         if round(dts(SetIndex),3) >= round( 2*ScalingDeltaT,3)
                    %             CompiledParameters.ScaledInitiationRates(SetIndex,IncludedHalfSampledBins, APbin) = ...
                    %                 SampleVector(IncludedHalfSampledBins);
                    %         else
                    %             CompiledParameters.ScaledInitiationRates(SetIndex,:, APbin) = SampleVector;
                    %         end
                    
                    
                    
                    CompiledParameters.ScaledInitiationRates(SetIndex,:, APbin) = SampleVector;
                    CompiledParameters.ScaledInitiationRatesStdErr(SetIndex,:, APbin) = ...
                        InterpolateStdError(SampleTimeVector, squeeze(CompiledParameters.InitiationRatesStdErr(SetIndex,:,APbin)), ReferenceTimeVector);
                    
                    CompiledParameters.ScaledInitiationRatesV2(SetIndex,:, APbin) = ...
                        CompiledParameters.ScaledInitiationRates(SetIndex,:, APbin)/t_coeff;
                    CompiledParameters.ScaledInitiationRatesStdErrV2(SetIndex,:, APbin) = ...
                        CompiledParameters.ScaledInitiationRatesStdErr(SetIndex,:, APbin)/t_coeff;
                    
                end
                
                for APbin = 1:NumAPbins
                    if all(isnan(CompiledParameters.Durations(SetIndex,:,APbin)))
                        continue
                    end
                    
                    SampleVector = interp1(SampleTimeVector, squeeze(CompiledParameters.Durations(SetIndex,:,APbin)), ReferenceTimeVector);
                    %if dts(SetIndex)
                    %         if round(dts(SetIndex),3) >= round( 2*ScalingDeltaT,3)
                    %             CompiledParameters.ScaledDurations(SetIndex,IncludedHalfSampledBins, APbin) = ...
                    %                 SampleVector(IncludedHalfSampledBins);
                    %         else
                    %             CompiledParameters.ScaledDurations(SetIndex,:, APbin) = SampleVector;
                    %         end
                    
                    CompiledParameters.ScaledDurations(SetIndex,:, APbin) = SampleVector;
                    CompiledParameters.ScaledDurationsStdErr(SetIndex,:, APbin) = ...
                        InterpolateStdError(SampleTimeVector, squeeze(CompiledParameters.DurationsStdErr(SetIndex,:,APbin)), ReferenceTimeVector);
                    CompiledParameters.ScaledDurationsV2(SetIndex,:, APbin) = ...
                        CompiledParameters.ScaledDurations(SetIndex,:, APbin)*t_coeff;
                    CompiledParameters.ScaledDurationsStdErrV2(SetIndex,:, APbin) = ...
                        CompiledParameters.ScaledDurationsStdErr(SetIndex,:, APbin)*t_coeff;
                end
                
                for APbin = 1:NumAPbins
                    if all(isnan(CompiledParameters.Frequencies(SetIndex,:,APbin)))
                        continue
                    end
                    
                    SampleVector = interp1(SampleTimeVector, squeeze(CompiledParameters.Frequencies(SetIndex,:,APbin)), ReferenceTimeVector);
                    %if dts(SetIndex)
                    %         if round(dts(SetIndex),3) >= round( 2*ScalingDeltaT,3)
                    %             CompiledParameters.ScaledFrequencies(SetIndex,IncludedHalfSampledBins, APbin) = ...
                    %                 SampleVector(IncludedHalfSampledBins);
                    %         else
                    %             CompiledParameters.ScaledFrequencies(SetIndex,:, APbin) = SampleVector;
                    %         end
                    %
                    CompiledParameters.ScaledFrequencies(SetIndex,:, APbin) = SampleVector;
                    CompiledParameters.ScaledFrequenciesStdErr(SetIndex,:, APbin) = ...
                        InterpolateStdError(SampleTimeVector, squeeze(CompiledParameters.FrequenciesStdErr(SetIndex,:,APbin)), ReferenceTimeVector);
                    CompiledParameters.ScaledFrequenciesV2(SetIndex,:, APbin) = ...
                        CompiledParameters.ScaledFrequencies(SetIndex,:, APbin)/t_coeff;
                    CompiledParameters.ScaledFrequenciesStdErrV2(SetIndex,:, APbin) = ...
                        CompiledParameters.ScaledFrequenciesStdErr(SetIndex,:, APbin)/t_coeff;
                end
                
            end
        else
            for SetIndex = 1:NumSets
                t_index = find(round(CompiledParameters.UniqueTemperatures, 1) == round(CompiledParameters.SetTemperatures(SetIndex), 1));
                t_coeff = CompiledParameters.TimingCoeffs(t_index);
                for APbin = 1:NumAPbins
                    if all(isnan(CompiledParameters.InitiationRates(SetIndex,:,APbin)))
                        continue
                    end
                    
                    
                    
                    
                    CompiledParameters.ScaledInitiationRatesV2(SetIndex,:, APbin) = ...
                        CompiledParameters.InitiationRates(SetIndex,:, APbin)/t_coeff;
                    CompiledParameters.ScaledInitiationRatesStdErrV2(SetIndex,:, APbin) = ...
                        CompiledParameters.InitiationRatesStdErr(SetIndex,:, APbin)/t_coeff;
                    CompiledParameters.ScaledDurationsV2(SetIndex,:, APbin) = ...
                        CompiledParameters.Durations(SetIndex,:, APbin)*t_coeff;
                    CompiledParameters.ScaledDurationsStdErrV2(SetIndex,:, APbin) = ...
                        CompiledParameters.DurationsStdErr(SetIndex,:, APbin)*t_coeff;
                    CompiledParameters.ScaledFrequenciesV2(SetIndex,:, APbin) = ...
                        CompiledParameters.Frequencies(SetIndex,:, APbin)/t_coeff;
                    CompiledParameters.ScaledFrequenciesStdErrV2(SetIndex,:, APbin) = ...
                        CompiledParameters.FrequenciesStdErr(SetIndex,:, APbin)/t_coeff;
                    
                end
                
                
            end
        end
        
    end
    
else
    APResolution = 0.025;
    APbins = 0:APResolution:1;
    NumAPbins = length(APbins);
    
    EnrichmentPath = getProcessedEnrichmentFolder();
    
    NumSets = length(ResultsPaths);
    
    
    CompiledParameters = {};
    CompiledParameters.ResultsPaths = ResultsPaths;
    CompiledParameters.FigurePaths = ResultsPaths;
    CompiledParameters.ReporterLabels = cell(1, NumSets);
    
    CompiledParameters.APVector = APbins;
    
    ScalingDeltaT = 5; % min
    
    
    
    
    
    
    
    
    
    %%
    
    prefix_regex = '(?<projectDirInfo>[A-Za-z0-9\-_\\]+)\\cpHMM_results\\compiledResults_w(?<nSteps>[0-9]+)_K(?<nStates>[0-9]+)_p0_ap(?<nAPbins>[0-9]+)_t(?<nTimebins>[0-9]+)_f(?<FluoDim>[0-9]+)D(?<TimeRes>[a-z0-9._]+)';
    cycle_regex = '(?<projectInfo>[A-Za-z\-0-9_]+)\\NC(?<cycle>[0-9]+)';
    timeres_regex = '_dt(?<dt>[0-9]+).mat';
    projectInfo_regex = 'hbBAC-MS2-(?<temperature>[0-9_]+)C';
    
    
    %%
    for ResultIndex=1:NumSets
        results_path_split = regexp(ResultsPaths{ResultIndex}, prefix_regex, 'names');
        if ~isempty(results_path_split)
            [fp1, fp2, fp3] = fileparts(ResultsPaths{ResultIndex});
            figpath = [EnrichmentPath,  filesep, fp1, filesep, 'figures', filesep, fp2, filesep];
            if ~exist(figpath, 'dir')
                mkdir(figpath);
            end
            CompiledParameters.FigurePaths{ResultIndex} = figpath;
            
            cycle_split = regexp(results_path_split.projectDirInfo, cycle_regex, 'names');
            if ~isempty(cycle_split)
                CompiledParameters.NC(ResultIndex) = str2num(cycle_split.cycle);
                proj_split_v1 = regexp(cycle_split.projectInfo, projectInfo_regex, 'names');
                PotentialProjName = cycle_split.projectInfo;
            else
                CompiledParameters.NC(ResultIndex) = 14;
                proj_split_v1 = regexp(results_path_split.projectDirInfo, projectInfo_regex, 'names');
                PotentialProjName = results_path_split.projectDirInfo;
            end
            
            if ~isempty(proj_split_v1)
                CompiledParameters.ReporterLabels{ResultIndex} = 'hbBAC-MS2';
                CompiledParameters.SetTemperatures(ResultIndex) = str2num(strrep(proj_split_v1.temperature, '_', '.'));
            else
                CompiledParameters.ReporterLabels{ResultIndex} = PotentialProjName;
                CompiledParameters.SetTemperatures(ResultIndex) = 25;
            end
            
            
            
            CompiledParameters.nSteps(ResultIndex) = str2num(results_path_split.nSteps);
            CompiledParameters.nStates(ResultIndex) = str2num(results_path_split.nStates);
            CompiledParameters.nAPbins(ResultIndex) = str2num(results_path_split.nAPbins);
            CompiledParameters.nTimebins(ResultIndex) = str2num(results_path_split.nTimebins);
            CompiledParameters.nTimebins(ResultIndex) = str2num(results_path_split.nTimebins);
            CompiledParameters.FluoDims(ResultIndex) = str2num(results_path_split.FluoDim);
            timeres_split = regexp(results_path_split.TimeRes, timeres_regex, 'names');
            if ~isempty(timeres_split)
                CompiledParameters.dt(ResultIndex) = str2num(timeres_split.dt);
            else
                datadir = strsplit(ResultsPaths{ResultIndex}, '\\cpHMM_results');
                load([EnrichmentPath datadir{1} filesep 'spot_struct.mat']);
                for i = 1:length(spot_struct)
                    if length(spot_struct(i).timeInterp) > 1
                        if all(~isnan(spot_struct(i).timeInterp))
                            CompiledParameters.dt(ResultIndex) = spot_struct(i).timeInterp(2)-spot_struct(i).timeInterp(1);
                            break
                        end
                    end
                    
                end
            end
            
            CompiledParameters.ElongationTimes(ResultIndex) = CompiledParameters.nSteps(ResultIndex)*CompiledParameters.dt(ResultIndex)/60;
            
            
        else
            error('Inference Info path has unexpected format.')
        end
    end
    %%
    [RefTemperature, RefIndex] = max(CompiledParameters.SetTemperatures);
    NumTimeBins = max(CompiledParameters.nTimebins);
    CompiledParameters.TimeLimits = NaN(NumTimeBins, 2);
    CompiledParameters.TimeVector = NaN(1, NumTimeBins);
    load([EnrichmentPath, ResultsPaths{RefIndex}]);
    for i = 1:length(compiledResults.timeBins)
        CompiledParameters.TimeLimits(i, 1) = compiledResults.timeBins{i}(1)/60;
        CompiledParameters.TimeLimits(i, 2) = compiledResults.timeBins{i}(2)/60;
        CompiledParameters.TimeVector(i) = mean(compiledResults.timeBins{i})/60;
    end
    
    
    CompiledParameters.TimeBinWidths = NaN(1, NumSets);
    CompiledParameters.TimeBinSeps = NaN(1, NumSets);
    
    
    CompiledParameters.APRange =NaN(1, 2);
    CompiledParameters.SetTemperatures = NaN(1,NumSets);
    CompiledParameters.NC = NaN(1,NumSets);
    CompiledParameters.nSteps = NaN(1,NumSets);
    CompiledParameters.nStates = NaN(1,NumSets);
    CompiledParameters.nAPbins = NaN(1,NumSets);
    CompiledParameters.nTimebins = NaN(1,NumSets);
    CompiledParameters.FluoDims = NaN(1, NumSets);
    CompiledParameters.dt = NaN(1,NumSets);
    CompiledParameters.TimeAxes = cell(1, NumSets);
    CompiledParameters.ElongationTimes = NaN(1,NumSets);
    CompiledParameters.Durations = NaN(NumSets, NumTimeBins, NumAPbins);
    CompiledParameters.Frequencies = NaN(NumSets, NumTimeBins, NumAPbins);
    CompiledParameters.InitiationRates = NaN(NumSets, NumTimeBins, NumAPbins);
    CompiledParameters.DurationsStdErr = NaN(NumSets, NumTimeBins, NumAPbins);
    CompiledParameters.FrequenciesStdErr = NaN(NumSets, NumTimeBins, NumAPbins);
    CompiledParameters.InitiationRatesStdErr = NaN(NumSets, NumTimeBins, NumAPbins);
    
    %%
    
    for ResultIndex=1:NumSets
        results_path_split = regexp(ResultsPaths{ResultIndex}, prefix_regex, 'names');
        if ~isempty(results_path_split)
            [fp1, fp2, fp3] = fileparts(ResultsPaths{ResultIndex});
            figpath = [EnrichmentPath,  filesep, fp1, filesep, 'figures', filesep, fp2, filesep];
            if ~exist(figpath, 'dir')
                mkdir(figpath);
            end
            CompiledParameters.FigurePaths{ResultIndex} = figpath;
            
            cycle_split = regexp(results_path_split.projectDirInfo, cycle_regex, 'names');
            if ~isempty(cycle_split)
                CompiledParameters.NC(ResultIndex) = str2num(cycle_split.cycle);
                proj_split_v1 = regexp(cycle_split.projectInfo, projectInfo_regex, 'names');
                PotentialProjName = cycle_split.projectInfo;
            else
                CompiledParameters.NC(ResultIndex) = 14;
                proj_split_v1 = regexp(results_path_split.projectDirInfo, projectInfo_regex, 'names');
                PotentialProjName = results_path_split.projectDirInfo;
            end
            
            if ~isempty(proj_split_v1)
                CompiledParameters.ReporterLabels{ResultIndex} = 'hbBAC-MS2';
                CompiledParameters.SetTemperatures(ResultIndex) = str2num(strrep(proj_split_v1.temperature, '_', '.'));
            else
                CompiledParameters.ReporterLabels{ResultIndex} = PotentialProjName;
                CompiledParameters.SetTemperatures(ResultIndex) = 25;
            end
            
            
            
            CompiledParameters.nSteps(ResultIndex) = str2num(results_path_split.nSteps);
            CompiledParameters.nStates(ResultIndex) = str2num(results_path_split.nStates);
            CompiledParameters.nAPbins(ResultIndex) = str2num(results_path_split.nAPbins);
            CompiledParameters.nTimebins(ResultIndex) = str2num(results_path_split.nTimebins);
            CompiledParameters.nTimebins(ResultIndex) = str2num(results_path_split.nTimebins);
            CompiledParameters.FluoDims(ResultIndex) = str2num(results_path_split.FluoDim);
            timeres_split = regexp(results_path_split.TimeRes, timeres_regex, 'names');
            if ~isempty(timeres_split)
                CompiledParameters.dt(ResultIndex) = str2num(timeres_split.dt);
            else
                datadir = strsplit(ResultsPaths{ResultIndex}, '\\cpHMM_results');
                load([EnrichmentPath datadir{1} filesep 'spot_struct.mat']);
                for i = 1:length(spot_struct)
                    if length(spot_struct(i).timeInterp) > 1
                        if all(~isnan(spot_struct(i).timeInterp))
                            CompiledParameters.dt(ResultIndex) = spot_struct(i).timeInterp(2)-spot_struct(i).timeInterp(1);
                            break
                        end
                    end
                    
                end
            end
            
            CompiledParameters.ElongationTimes(ResultIndex) = CompiledParameters.nSteps(ResultIndex)*CompiledParameters.dt(ResultIndex)/60;
            
            
        else
            error('Inference Info path has unexpected format.')
        end
        
        load([EnrichmentPath, ResultsPaths{ResultIndex}]);
        
        time_group_vec = compiledResults.timeGroupVec;
        time_group_index = unique(time_group_vec);
        ap_group_vec = compiledResults.apGroupVec;
        ap_group_index = unique(ap_group_vec);
        ap_axis = compiledResults.apBins(1:end-1) + diff(compiledResults.apBins);
        time_axis = [];
        for t = time_group_index
            time_axis(end+1) = mean(compiledResults.timeBins{t})/60;
        end
        CompiledParameters.TimeAxes{ResultIndex} = time_axis;
        
        
        
        
        CompiledParameters.TimeBinWidths(ResultIndex) = CompiledParameters.TimeAxes{ResultIndex}(1)*2;
        CompiledParameters.TimeBinSeps(ResultIndex) = CompiledParameters.TimeAxes{ResultIndex}(2)-CompiledParameters.TimeAxes{ResultIndex}(1);
        
        for t = time_group_index
            time_filter = time_group_vec==t;
            ap_ids = ap_group_vec(time_filter);
            ap_bins = ap_axis(ap_ids)/100;
            [ap_bins, sort_order] = sort(ap_bins);
            AP_index = ismember(round(APbins, 3), round(ap_bins, 3));
            
            dur_vec_mean = compiledResults.dur_vec_mean(time_filter);
            dur_vec_mean = dur_vec_mean(sort_order);
            dur_vec_ste = compiledResults.dur_vec_ste(time_filter);
            dur_vec_ste = dur_vec_ste(sort_order);
            
            CompiledParameters.Durations(ResultIndex, t, AP_index) = dur_vec_mean;
            CompiledParameters.DurationsStdErr(ResultIndex, t, AP_index) = dur_vec_ste;
            
            freq_vec_mean = compiledResults.freq_vec_mean(time_filter);
            freq_vec_mean = freq_vec_mean(sort_order);
            freq_vec_ste = compiledResults.freq_vec_ste(time_filter);
            freq_vec_ste = freq_vec_ste(sort_order);
            
            CompiledParameters.Frequencies(ResultIndex, t, AP_index) = freq_vec_mean;
            CompiledParameters.FrequenciesStdErr(ResultIndex, t, AP_index) = freq_vec_ste;
            
            
            init_vec_mean = compiledResults.init_vec_mean(time_filter);
            init_vec_mean = init_vec_mean(sort_order);
            init_vec_ste = compiledResults.init_vec_ste(time_filter);
            init_vec_ste = init_vec_ste(sort_order);
            
            
            CompiledParameters.InitiationRates(ResultIndex, t, AP_index) = init_vec_mean;
            CompiledParameters.InitiationRatesStdErr(ResultIndex, t, AP_index) = init_vec_ste;
        end
        
    end
    
    
    %%
    if includeRescaling
        CompiledParameters.UniqueTemperatures = fliplr(unique(CompiledParameters.SetTemperatures));
        
        if ~isempty(TempParamPath)
            if isfile([TempParamPath filesep 'autocorrElongationTimes.mat'])
                load([TempParamPath filesep 'autocorrElongationTimes.mat']);
                ElongTemps =  ElongationTimeInfo.Temperatures;
                CompiledParameters.AutocorrElongationTimes = NaN(1, length(CompiledParameters.UniqueTemperatures));
                for i = 1:length(ElongTemps)
                    match_index = find(round(CompiledParameters.UniqueTemperatures, 1) == round(ElongTemps(i), 1));
                    if ~isempty(match_index)
                        CompiledParameters.AutocorrElongationTimes(match_index) = ElongationTimeInfo.ElongationTimes(i);
                    end
                end
            else
                CompiledParameters.AutocorrElongationTimes = NaN(1, length(CompiledParameters.UniqueTemperatures));
            end
            
            if isfile([TempParamPath filesep 'DevTimeCoefficients.mat'])
                load([TempParamPath filesep 'DevTimeCoefficients.mat']);
                DevTimeTemps =  TimingInfo.Temperatures;
                CompiledParameters.TimingCoeffs = NaN(1, length(CompiledParameters.UniqueTemperatures));
                for i = 1:length(DevTimeTemps)
                    match_index = find(round(CompiledParameters.UniqueTemperatures, 1) == round(DevTimeTemps(i), 1));
                    if ~isempty(match_index)
                        CompiledParameters.TimingCoeffs(match_index) = TimingInfo.TimeCoeffs(i);
                    end
                end
            else
                CompiledParameters.TimingCoeffs = NaN(1, length(CompiledParameters.UniqueTemperatures));
            end
            
            
            
            
        else
            CompiledParameters.AutocorrElongationTimes = NaN(1, length(CompiledParameters.UniqueTemperatures));
            CompiledParameters.TimingCoeffs = NaN(1, length(CompiledParameters.UniqueTemperatures));
            
        end
        
        
        if ~isempty(FluoParamPath)
            if isfile([FluoParamPath filesep 'FluoCoefficients.mat'])
                load([FluoParamPath filesep 'FluoCoefficients.mat']);
                FluoTemps =  FluoInfo.Temperatures;
                CompiledParameters.FluoCoeffs = NaN(1, length(CompiledParameters.UniqueTemperatures));
                for i = 1:length(FluoTemps)
                    match_index = find(round(CompiledParameters.UniqueTemperatures, 1) == round(FluoTemps(i), 1));
                    if ~isempty(match_index)
                        CompiledParameters.FluoCoeffs(match_index) = FluoInfo.FluoCoeffs(i);
                    end
                end
            else
                CompiledParameters.FluoCoeffs = NaN(1, length(CompiledParameters.UniqueTemperatures));
            end
            
        else
            CompiledParameters.FluoCoeffs = NaN(1, length(CompiledParameters.UniqueTemperatures));
            
        end
        
        %%
        
        ReferenceTimeVector = CompiledParameters.TimeVector;
        
        
        
        %%
        CompiledParameters.ScaledTimeVector = ReferenceTimeVector;
        CompiledParameters.ScaledDurations = CompiledParameters.Durations;
        CompiledParameters.ScaledFrequencies = CompiledParameters.Frequencies;
        CompiledParameters.ScaledInitiationRates = CompiledParameters.InitiationRates;
        CompiledParameters.ScaledDurationsStdErr = CompiledParameters.DurationsStdErr;
        CompiledParameters.ScaledFrequenciesStdErr = CompiledParameters.FrequenciesStdErr;
        CompiledParameters.ScaledInitiationRatesStdErr = CompiledParameters.InitiationRatesStdErr;
        
        CompiledParameters.ScaledDurationsV2 = NaN(size(CompiledParameters.ScaledDurations));
        CompiledParameters.ScaledFrequenciesV2 = NaN(size(CompiledParameters.ScaledDurations));
        CompiledParameters.ScaledInitiationRatesV2 = NaN(size(CompiledParameters.ScaledDurations));
        CompiledParameters.ScaledDurationsStdErrV2 = NaN(size(CompiledParameters.ScaledDurations));
        CompiledParameters.ScaledFrequenciesStdErrV2 = NaN(size(CompiledParameters.ScaledDurations));
        CompiledParameters.ScaledInitiationRatesStdErrV2 = NaN(size(CompiledParameters.ScaledDurations));
        
        
        %%
        if NumTimeBins > 1
            for SetIndex = 1:NumSets
                t_index = find(round(CompiledParameters.UniqueTemperatures, 1) == round(CompiledParameters.SetTemperatures(SetIndex), 1));
                t_coeff = CompiledParameters.TimingCoeffs(t_index);
                SampleTimeVector = CompiledParameters.TimeVector*t_coeff;
                for APbin = 1:NumAPbins
                    if all(isnan(CompiledParameters.InitiationRates(SetIndex,:,APbin)))
                        continue
                    end
                    
                    
                    
                    CompiledParameters.ScaledInitiationRatesV2(SetIndex,:, APbin) = ...
                        CompiledParameters.ScaledInitiationRates(SetIndex,:, APbin)/t_coeff;
                    CompiledParameters.ScaledInitiationRatesStdErrV2(SetIndex,:, APbin) = ...
                        CompiledParameters.ScaledInitiationRatesStdErr(SetIndex,:, APbin)/t_coeff;
                    
                end
                
                for APbin = 1:NumAPbins
                    if all(isnan(CompiledParameters.Durations(SetIndex,:,APbin)))
                        continue
                    end
                    
                    CompiledParameters.ScaledDurationsV2(SetIndex,:, APbin) = ...
                        CompiledParameters.ScaledDurations(SetIndex,:, APbin)*t_coeff;
                    CompiledParameters.ScaledDurationsStdErrV2(SetIndex,:, APbin) = ...
                        CompiledParameters.ScaledDurationsStdErr(SetIndex,:, APbin)*t_coeff;
                end
                
                for APbin = 1:NumAPbins
                    if all(isnan(CompiledParameters.Frequencies(SetIndex,:,APbin)))
                        continue
                    end
                    
                    
                    CompiledParameters.ScaledFrequenciesV2(SetIndex,:, APbin) = ...
                        CompiledParameters.ScaledFrequencies(SetIndex,:, APbin)/t_coeff;
                    CompiledParameters.ScaledFrequenciesStdErrV2(SetIndex,:, APbin) = ...
                        CompiledParameters.ScaledFrequenciesStdErr(SetIndex,:, APbin)/t_coeff;
                end
                
            end
        else
            for SetIndex = 1:NumSets
                t_index = find(round(CompiledParameters.UniqueTemperatures, 1) == round(CompiledParameters.SetTemperatures(SetIndex), 1));
                t_coeff = CompiledParameters.TimingCoeffs(t_index);
                for APbin = 1:NumAPbins
                    if all(isnan(CompiledParameters.InitiationRates(SetIndex,:,APbin)))
                        continue
                    end
                    
                    
                    
                    
                    CompiledParameters.ScaledInitiationRatesV2(SetIndex,:, APbin) = ...
                        CompiledParameters.InitiationRates(SetIndex,:, APbin)/t_coeff;
                    CompiledParameters.ScaledInitiationRatesStdErrV2(SetIndex,:, APbin) = ...
                        CompiledParameters.InitiationRatesStdErr(SetIndex,:, APbin)/t_coeff;
                    CompiledParameters.ScaledDurationsV2(SetIndex,:, APbin) = ...
                        CompiledParameters.Durations(SetIndex,:, APbin)*t_coeff;
                    CompiledParameters.ScaledDurationsStdErrV2(SetIndex,:, APbin) = ...
                        CompiledParameters.DurationsStdErr(SetIndex,:, APbin)*t_coeff;
                    CompiledParameters.ScaledFrequenciesV2(SetIndex,:, APbin) = ...
                        CompiledParameters.Frequencies(SetIndex,:, APbin)/t_coeff;
                    CompiledParameters.ScaledFrequenciesStdErrV2(SetIndex,:, APbin) = ...
                        CompiledParameters.FrequenciesStdErr(SetIndex,:, APbin)/t_coeff;
                    
                end
                
                
            end
        end
        
        
        
    end
    
    
    
    CompiledParameters.Durations = NaN(size(CompiledParameters.ScaledDurations));
    CompiledParameters.Frequencies = NaN(size(CompiledParameters.ScaledDurations));
    CompiledParameters.InitiationRates = NaN(size(CompiledParameters.ScaledDurations));
    CompiledParameters.DurationsStdErr = NaN(size(CompiledParameters.ScaledDurations));
    CompiledParameters.FrequenciesStdErr = NaN(size(CompiledParameters.ScaledDurations));
    CompiledParameters.InitiationRatesStdErr = NaN(size(CompiledParameters.ScaledDurations));
    
    
    
end
CompiledParameters.RefTemperature = RefTemperature;
try
CompiledParameters = AddHmmActivationEnergies(CompiledParameters, includeRescaling);


end












