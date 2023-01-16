%function IntegratedTraceFluos = CalculateEmbryoIntegratedTraceFluos(ltm, IncludeFractionOn, UseRescaledFluo)
% Hardcoded in anaphase aligned trace type
if ~exist('IncludeFractionOn', 'var')
    IncludeFractionOn = false;
end
if ~exist('UseRescaledFluo', 'var')
    UseRescaledFluo = false;
end

Temperatures = fliplr(unique(ltm.Temp_sps));
NumTemperatures = length(Temperatures);
EmbryoTemperatures = ltm.Temp_sps;
NumEmbryos = length(ltm.Temp_sps);
IncludeExperimentTF = ismember(1:NumEmbryos, ltm.IncludedExperiments);


NumAPbins = size(ltm.MeanProfiles{1}.AnaphaseAlignedCycleMeanTraces, 2);
EmbryoIndices = find(IncludeExperimentTF);
IntegratedFluos = NaN(NumEmbryos, NumAPbins, 6); % AU*min
for EmbryoIndex = 1:NumEmbryos
    if ismember(EmbryoIndex, EmbryoIndices)
        for NC = 9:14
            nc_index = NC-8;
            FrameTimes = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleFrameTimes{nc_index}/60;
            AllAPProfiles = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleMeanPerNucleusTraces(:,:,nc_index);
            if IncludeFractionOn
                 AllFractionOns = ltm.MeanProfiles{EmbryoIndex}.AnaphaseAlignedCycleFractionOn(:,:,nc_index);
            end
            for APindex = 1:NumAPbins
                FullAPProfile = AllAPProfiles(:,APindex).';
                UsableTimes = ~isnan(FullAPProfile);
               
                APFrameTimes =  FrameTimes(UsableTimes);
                APProfile = FullAPProfile(UsableTimes);
                if sum(UsableTimes) == 0
                    continue
                end
                if IncludeFractionOn
                    FullFractionOnProfile = AllFractionOns(:,APindex).';
                    FractionOnProfile = FullFractionOnProfile(UsableTimes);
                    APProfile = APProfile.*FractionOnProfile;
                end
                TimeBinAveragedProfile = (APProfile(1:end-1)+APProfile(2:end))/2;
                TimeBinWidths = diff(APFrameTimes);
                IntegratedFluos(EmbryoIndex, APindex, NC-8) = sum(TimeBinAveragedProfile.*TimeBinWidths);
            end
        end
    end
end

[PlottedParams, PlottedParamSEs,R2s, ylab,OutputString,GlobalPlotYmax,GlobalPlotYmin,LogPlotYmin] = ...
    getPlottingVariables(ltm, parameter,  TraceType, R2bound, UseRescaledFluo, UseRescaledTiming,...
    UsePerNucleusTraces, UseBinnedTraces, UseBinnedPerNucleusTraces);

%% First make average profiles binning anaphase aligned traces for each temperature 
