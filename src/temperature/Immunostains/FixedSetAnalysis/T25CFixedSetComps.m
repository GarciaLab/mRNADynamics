clear all
SizeDataPath = 'S:/Gabriella/Dropbox/EmbryoSizeMeasurements/EmbryoSizeData.mat';
load(SizeDataPath,'APLengths','DVLengths','VentralDistances','DorsalDistances','NGoodEmbryos',...
    'MeanAPLength','StdAPLength','MeanDVLength','StdDVLength','AspectRatios');
MeanAspectRatio = mean(AspectRatios);
StdAspectRatio = std(AspectRatios);

T25CRep1SetPrefixes = {'2022-04-03-yw104-180m25C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide3',...
    '2022-04-04-yw104-180m25C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide10',...
    '2022-04-07-yw104-180m25C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide9',...
    '2022-04-09-yw104-180m25C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep1Slide11'};
T25CRep1SetLabel = 'FixedT25CReplicate1';

T25CRep2SetPrefixes = {'2022-05-07-yw124-180m25C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep2Slide4',...
    '2022-05-11-yw124-180m25C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep2Slide6',...
    '2022-05-12-yw124-180m25C_Bcd488_Hb546_Hoechst_yw174-190m25C_Kni647_Rep2Slide5'};
T25CRep2SetLabel = 'FixedT25CReplicate2';


T25CFlippedSetPrefixes = {'2022-05-02-yw164-180m25C_Kni647_yw174-190m25C_Bcd488_Hb546_Hoechst_Flip1Slide1',...
    '2022-05-03-yw164-180m25C_Kni647_yw174-190m25C_Bcd488_Hb546_Hoechst_Flip1Slide2',...
    '2022-05-05-yw164-180m25C_Kni647_yw174-190m25C_Bcd488_Hb546_Hoechst_Flip1Slide4'};
T25CFlippedSetLabel = 'FixedT25CFlipped';

%%
T25CRep1_liveExperiments = cell(1, length(T25CRep1SetPrefixes));
for i = 1:length(T25CRep1SetPrefixes)
    T25CRep1_liveExperiments{i} = LiveExperiment(T25CRep1SetPrefixes{i});
end
FixedPixelSize_um = T25CRep1_liveExperiments{1}.pixelSize_um;
T25CRep1_CompiledEmbryos = CombineCompiledEmbryos(T25CRep1SetPrefixes);
T25CRep1_ApprovedEmbryos = T25CRep1_CompiledEmbryos.Approved;
T25CRep1_NEmbryos = length(T25CRep1_ApprovedEmbryos);

%%

T25CRep1_FixedAPLengths = NaN(1, T25CRep1_NEmbryos);
T25CRep1_FixedAPSlopes =  NaN(1, T25CRep1_NEmbryos);
T25CRep1_FixedAPIntercepts =  NaN(1, T25CRep1_NEmbryos);
T25CRep1_FixedDorsalDistances =  NaN(1, T25CRep1_NEmbryos);
T25CRep1_FixedVentralDistances =  NaN(1, T25CRep1_NEmbryos);
T25CRep1_FixedDVLengths =  NaN(1, T25CRep1_NEmbryos);

for i = 1:T25CRep1_NEmbryos
    if T25CRep1_CompiledEmbryos.Approved(i) & T25CRep1_CompiledEmbryos.Flags(i) == 0
        T25CRep1_FixedAPLengths(i) = FixedPixelSize_um*sqrt((T25CRep1_CompiledEmbryos.CoordAs(i,1)-T25CRep1_CompiledEmbryos.CoordPs(i,1))^2+...
            (T25CRep1_CompiledEmbryos.CoordAs(i,2)-T25CRep1_CompiledEmbryos.CoordPs(i,2))^2);
        T25CRep1_FixedAPSlopes(i) = (T25CRep1_CompiledEmbryos.CoordAs(i,2)-T25CRep1_CompiledEmbryos.CoordPs(i,2))/(T25CRep1_CompiledEmbryos.CoordAs(i,1)-T25CRep1_CompiledEmbryos.CoordPs(i,1));
        T25CRep1_FixedAPIntercepts(i) = T25CRep1_CompiledEmbryos.CoordAs(i,2)-T25CRep1_FixedAPSlopes(i)*T25CRep1_CompiledEmbryos.CoordAs(i,1);
        T25CRep1_FixedDorsalDistances(i) = FixedPixelSize_um*abs(T25CRep1_FixedAPIntercepts(i)+T25CRep1_FixedAPSlopes(i)*T25CRep1_CompiledEmbryos.CoordDs(i,1)-T25CRep1_CompiledEmbryos.CoordDs(i,2))/sqrt(T25CRep1_FixedAPSlopes(i)^2+1);
        T25CRep1_FixedVentralDistances(i) = FixedPixelSize_um*abs(T25CRep1_FixedAPIntercepts(i)+T25CRep1_FixedAPSlopes(i)*T25CRep1_CompiledEmbryos.CoordVs(i,1)-T25CRep1_CompiledEmbryos.CoordVs(i,2))/sqrt(T25CRep1_FixedAPSlopes(i)^2+1);
        T25CRep1_FixedDVLengths(i) = T25CRep1_FixedDorsalDistances(i)+T25CRep1_FixedVentralDistances(i);
    end
end

T25CRep1_CompiledEmbryos.APLengths = T25CRep1_FixedAPLengths;
T25CRep1_CompiledEmbryos.APSlopes = T25CRep1_FixedAPSlopes;
T25CRep1_CompiledEmbryos.APIntercepts = T25CRep1_FixedAPIntercepts;
T25CRep1_CompiledEmbryos.DorsalDistances = T25CRep1_FixedDorsalDistances;
T25CRep1_CompiledEmbryos.VentralDistances = T25CRep1_FixedVentralDistances;
T25CRep1_CompiledEmbryos.DVLengths = T25CRep1_FixedDVLengths;

T25CRep1_FixedAspectRatios = T25CRep1_FixedDVLengths./T25CRep1_FixedAPLengths;
T25CRep1_CompiledEmbryos.AspectRatios = T25CRep1_FixedAspectRatios;
T25CRep1_MeanFixedAspectRatio = mean(T25CRep1_FixedAspectRatios, 'omitnan');
T25CRep1_StdFixedAspectRatio = std(T25CRep1_FixedAspectRatios, 'omitnan');

T25CRep1_MeanFixedAPLength = mean(T25CRep1_FixedAPLengths, 'omitnan');
T25CRep1_StdFixedAPLength = std(T25CRep1_FixedAPLengths, 'omitnan');
T25CRep1_MeanFixedDVLength = mean(T25CRep1_FixedDVLengths, 'omitnan');
T25CRep1_StdFixedDVLength = std(T25CRep1_FixedDVLengths, 'omitnan');


T25CRep1_FixedAPLengthFit = fitdist(max(T25CRep1_FixedAPLengths)-T25CRep1_FixedAPLengths.', 'logistic');
T25CRep1_FixedAPDistCenter = -1*T25CRep1_FixedAPLengthFit.ParameterValues(1)+max(T25CRep1_FixedAPLengths);
T25CRep1_FixedAPLimits = [T25CRep1_FixedAPDistCenter+(min(APLengths)-MeanAPLength)*T25CRep1_FixedAPDistCenter/MeanAPLength,...
    T25CRep1_FixedAPDistCenter+(max(APLengths)-MeanAPLength)*T25CRep1_FixedAPDistCenter/MeanAPLength];
%
ARLimits = [min(AspectRatios,[], 'omitnan'), max(AspectRatios,[], 'omitnan')];
T25CRep1_ARKeepEmbryos = T25CRep1_FixedAspectRatios >= ARLimits(1) & T25CRep1_FixedAspectRatios <= ARLimits(2);

T25CRep1_APKeepEmbryos = T25CRep1_FixedAPLengths >= T25CRep1_FixedAPLimits(1) & T25CRep1_FixedAPLengths <= T25CRep1_FixedAPLimits(2);

T25CRep1_CompiledEmbryos.EmbryoDimensionsGood = T25CRep1_APKeepEmbryos & T25CRep1_ARKeepEmbryos;


T25CRep1_NC13Indices = T25CRep1_CompiledEmbryos.EmbryoDimensionsGood & T25CRep1_CompiledEmbryos.Flags == 0 &...
    T25CRep1_CompiledEmbryos.Approved & T25CRep1_CompiledEmbryos.nc == 13;
T25CRep1_NC14Indices = T25CRep1_CompiledEmbryos.EmbryoDimensionsGood & T25CRep1_CompiledEmbryos.Flags == 0 &...
    T25CRep1_CompiledEmbryos.Approved & T25CRep1_CompiledEmbryos.nc == 14;
%%
%T25CRep1_CompiledEmbryos.
%%
T25CRep2_liveExperiments = cell(1, length(T25CRep2SetPrefixes));
for i = 1:length(T25CRep2SetPrefixes)
    T25CRep2_liveExperiments{i} = LiveExperiment(T25CRep2SetPrefixes{i});
end
FixedPixelSize_um = T25CRep2_liveExperiments{1}.pixelSize_um;
T25CRep2_CompiledEmbryos = CombineCompiledEmbryos(T25CRep2SetPrefixes);
T25CRep2_ApprovedEmbryos = T25CRep2_CompiledEmbryos.Approved;
T25CRep2_NEmbryos = length(T25CRep2_ApprovedEmbryos);

%%

T25CRep2_FixedAPLengths = NaN(1, T25CRep2_NEmbryos);
T25CRep2_FixedAPSlopes =  NaN(1, T25CRep2_NEmbryos);
T25CRep2_FixedAPIntercepts =  NaN(1, T25CRep2_NEmbryos);
T25CRep2_FixedDorsalDistances =  NaN(1, T25CRep2_NEmbryos);
T25CRep2_FixedVentralDistances =  NaN(1, T25CRep2_NEmbryos);
T25CRep2_FixedDVLengths =  NaN(1, T25CRep2_NEmbryos);

for i = 1:T25CRep2_NEmbryos
    if T25CRep2_CompiledEmbryos.Approved(i) & T25CRep2_CompiledEmbryos.Flags(i) == 0
        T25CRep2_FixedAPLengths(i) = FixedPixelSize_um*sqrt((T25CRep2_CompiledEmbryos.CoordAs(i,1)-T25CRep2_CompiledEmbryos.CoordPs(i,1))^2+...
            (T25CRep2_CompiledEmbryos.CoordAs(i,2)-T25CRep2_CompiledEmbryos.CoordPs(i,2))^2);
        T25CRep2_FixedAPSlopes(i) = (T25CRep2_CompiledEmbryos.CoordAs(i,2)-T25CRep2_CompiledEmbryos.CoordPs(i,2))/(T25CRep2_CompiledEmbryos.CoordAs(i,1)-T25CRep2_CompiledEmbryos.CoordPs(i,1));
        T25CRep2_FixedAPIntercepts(i) = T25CRep2_CompiledEmbryos.CoordAs(i,2)-T25CRep2_FixedAPSlopes(i)*T25CRep2_CompiledEmbryos.CoordAs(i,1);
        T25CRep2_FixedDorsalDistances(i) = FixedPixelSize_um*abs(T25CRep2_FixedAPIntercepts(i)+T25CRep2_FixedAPSlopes(i)*T25CRep2_CompiledEmbryos.CoordDs(i,1)-T25CRep2_CompiledEmbryos.CoordDs(i,2))/sqrt(T25CRep2_FixedAPSlopes(i)^2+1);
        T25CRep2_FixedVentralDistances(i) = FixedPixelSize_um*abs(T25CRep2_FixedAPIntercepts(i)+T25CRep2_FixedAPSlopes(i)*T25CRep2_CompiledEmbryos.CoordVs(i,1)-T25CRep2_CompiledEmbryos.CoordVs(i,2))/sqrt(T25CRep2_FixedAPSlopes(i)^2+1);
        T25CRep2_FixedDVLengths(i) = T25CRep2_FixedDorsalDistances(i)+T25CRep2_FixedVentralDistances(i);
    end
end

T25CRep2_CompiledEmbryos.APLengths = T25CRep2_FixedAPLengths;
T25CRep2_CompiledEmbryos.APSlopes = T25CRep2_FixedAPSlopes;
T25CRep2_CompiledEmbryos.APIntercepts = T25CRep2_FixedAPIntercepts;
T25CRep2_CompiledEmbryos.DorsalDistances = T25CRep2_FixedDorsalDistances;
T25CRep2_CompiledEmbryos.VentralDistances = T25CRep2_FixedVentralDistances;
T25CRep2_CompiledEmbryos.DVLengths = T25CRep2_FixedDVLengths;

T25CRep2_FixedAspectRatios = T25CRep2_FixedDVLengths./T25CRep2_FixedAPLengths;
T25CRep2_CompiledEmbryos.AspectRatios = T25CRep2_FixedAspectRatios;

T25CRep2_MeanFixedAspectRatio = mean(T25CRep2_FixedAspectRatios, 'omitnan');
T25CRep2_StdFixedAspectRatio = std(T25CRep2_FixedAspectRatios, 'omitnan');

T25CRep2_MeanFixedAPLength = mean(T25CRep2_FixedAPLengths, 'omitnan');
T25CRep2_StdFixedAPLength = std(T25CRep2_FixedAPLengths, 'omitnan');
T25CRep2_MeanFixedDVLength = mean(T25CRep2_FixedDVLengths, 'omitnan');
T25CRep2_StdFixedDVLength = std(T25CRep2_FixedDVLengths, 'omitnan');


T25CRep2_FixedAPLengthFit = fitdist(max(T25CRep2_FixedAPLengths)-T25CRep2_FixedAPLengths.', 'logistic');
T25CRep2_FixedAPDistCenter = -1*T25CRep2_FixedAPLengthFit.ParameterValues(1)+max(T25CRep2_FixedAPLengths);
T25CRep2_FixedAPLimits = [T25CRep2_FixedAPDistCenter+(min(APLengths)-MeanAPLength)*T25CRep2_FixedAPDistCenter/MeanAPLength,...
    T25CRep2_FixedAPDistCenter+(max(APLengths)-MeanAPLength)*T25CRep2_FixedAPDistCenter/MeanAPLength];
%
ARLimits = [min(AspectRatios,[], 'omitnan'), max(AspectRatios,[], 'omitnan')];
T25CRep2_ARKeepEmbryos = T25CRep2_FixedAspectRatios >= ARLimits(1) & T25CRep2_FixedAspectRatios <= ARLimits(2);

T25CRep2_APKeepEmbryos = T25CRep2_FixedAPLengths >= T25CRep2_FixedAPLimits(1) & T25CRep2_FixedAPLengths <= T25CRep2_FixedAPLimits(2);

T25CRep2_CompiledEmbryos.EmbryoDimensionsGood = T25CRep2_APKeepEmbryos & T25CRep2_ARKeepEmbryos;


T25CRep2_NC13Indices = T25CRep2_CompiledEmbryos.EmbryoDimensionsGood & T25CRep2_CompiledEmbryos.Flags == 0 &...
    T25CRep2_CompiledEmbryos.Approved & T25CRep2_CompiledEmbryos.nc == 13;
T25CRep2_NC14Indices = T25CRep2_CompiledEmbryos.EmbryoDimensionsGood & T25CRep2_CompiledEmbryos.Flags == 0 &...
    T25CRep2_CompiledEmbryos.Approved & T25CRep2_CompiledEmbryos.nc == 14;


%%
T25CFlipped_liveExperiments = cell(1, length(T25CFlippedSetPrefixes));
for i = 1:length(T25CFlippedSetPrefixes)
    T25CFlipped_liveExperiments{i} = LiveExperiment(T25CFlippedSetPrefixes{i});
end
FixedPixelSize_um = T25CFlipped_liveExperiments{1}.pixelSize_um;
T25CFlipped_CompiledEmbryos = CombineCompiledEmbryos(T25CFlippedSetPrefixes);
T25CFlipped_ApprovedEmbryos = T25CFlipped_CompiledEmbryos.Approved;
T25CFlipped_NEmbryos = length(T25CFlipped_ApprovedEmbryos);

%%

T25CFlipped_FixedAPLengths = NaN(1, T25CFlipped_NEmbryos);
T25CFlipped_FixedAPSlopes =  NaN(1, T25CFlipped_NEmbryos);
T25CFlipped_FixedAPIntercepts =  NaN(1, T25CFlipped_NEmbryos);
T25CFlipped_FixedDorsalDistances =  NaN(1, T25CFlipped_NEmbryos);
T25CFlipped_FixedVentralDistances =  NaN(1, T25CFlipped_NEmbryos);
T25CFlipped_FixedDVLengths =  NaN(1, T25CFlipped_NEmbryos);

for i = 1:T25CFlipped_NEmbryos
    if T25CFlipped_CompiledEmbryos.Approved(i) & T25CFlipped_CompiledEmbryos.Flags(i) == 0
        T25CFlipped_FixedAPLengths(i) = FixedPixelSize_um*sqrt((T25CFlipped_CompiledEmbryos.CoordAs(i,1)-T25CFlipped_CompiledEmbryos.CoordPs(i,1))^2+...
            (T25CFlipped_CompiledEmbryos.CoordAs(i,2)-T25CFlipped_CompiledEmbryos.CoordPs(i,2))^2);
        T25CFlipped_FixedAPSlopes(i) = (T25CFlipped_CompiledEmbryos.CoordAs(i,2)-T25CFlipped_CompiledEmbryos.CoordPs(i,2))/(T25CFlipped_CompiledEmbryos.CoordAs(i,1)-T25CFlipped_CompiledEmbryos.CoordPs(i,1));
        T25CFlipped_FixedAPIntercepts(i) = T25CFlipped_CompiledEmbryos.CoordAs(i,2)-T25CFlipped_FixedAPSlopes(i)*T25CFlipped_CompiledEmbryos.CoordAs(i,1);
        T25CFlipped_FixedDorsalDistances(i) = FixedPixelSize_um*abs(T25CFlipped_FixedAPIntercepts(i)+T25CFlipped_FixedAPSlopes(i)*T25CFlipped_CompiledEmbryos.CoordDs(i,1)-T25CFlipped_CompiledEmbryos.CoordDs(i,2))/sqrt(T25CFlipped_FixedAPSlopes(i)^2+1);
        T25CFlipped_FixedVentralDistances(i) = FixedPixelSize_um*abs(T25CFlipped_FixedAPIntercepts(i)+T25CFlipped_FixedAPSlopes(i)*T25CFlipped_CompiledEmbryos.CoordVs(i,1)-T25CFlipped_CompiledEmbryos.CoordVs(i,2))/sqrt(T25CFlipped_FixedAPSlopes(i)^2+1);
        T25CFlipped_FixedDVLengths(i) = T25CFlipped_FixedDorsalDistances(i)+T25CFlipped_FixedVentralDistances(i);
    end
end

T25CFlipped_CompiledEmbryos.APLengths = T25CFlipped_FixedAPLengths;
T25CFlipped_CompiledEmbryos.APSlopes = T25CFlipped_FixedAPSlopes;
T25CFlipped_CompiledEmbryos.APIntercepts = T25CFlipped_FixedAPIntercepts;
T25CFlipped_CompiledEmbryos.DorsalDistances = T25CFlipped_FixedDorsalDistances;
T25CFlipped_CompiledEmbryos.VentralDistances = T25CFlipped_FixedVentralDistances;
T25CFlipped_CompiledEmbryos.DVLengths = T25CFlipped_FixedDVLengths;

T25CFlipped_FixedAspectRatios = T25CFlipped_FixedDVLengths./T25CFlipped_FixedAPLengths;
T25CFlipped_CompiledEmbryos.AspectRatios = T25CFlipped_FixedAspectRatios;

T25CFlipped_MeanFixedAspectRatio = mean(T25CFlipped_FixedAspectRatios, 'omitnan');
T25CFlipped_StdFixedAspectRatio = std(T25CFlipped_FixedAspectRatios, 'omitnan');

T25CFlipped_MeanFixedAPLength = mean(T25CFlipped_FixedAPLengths, 'omitnan');
T25CFlipped_StdFixedAPLength = std(T25CFlipped_FixedAPLengths, 'omitnan');
T25CFlipped_MeanFixedDVLength = mean(T25CFlipped_FixedDVLengths, 'omitnan');
T25CFlipped_StdFixedDVLength = std(T25CFlipped_FixedDVLengths, 'omitnan');

T25CFlipped_FixedAPLengthFit = fitdist(max(T25CFlipped_FixedAPLengths)-T25CFlipped_FixedAPLengths.', 'logistic');
T25CFlipped_FixedAPDistCenter = -1*T25CFlipped_FixedAPLengthFit.ParameterValues(1)+max(T25CFlipped_FixedAPLengths);
T25CFlipped_FixedAPLimits = [T25CFlipped_FixedAPDistCenter+(min(APLengths)-MeanAPLength)*T25CFlipped_FixedAPDistCenter/MeanAPLength,...
    T25CFlipped_FixedAPDistCenter+(max(APLengths)-MeanAPLength)*T25CFlipped_FixedAPDistCenter/MeanAPLength];
%
ARLimits = [min(AspectRatios,[], 'omitnan'), max(AspectRatios,[], 'omitnan')];
T25CFlipped_ARKeepEmbryos = T25CFlipped_FixedAspectRatios >= ARLimits(1) & T25CFlipped_FixedAspectRatios <= ARLimits(2);

T25CFlipped_APKeepEmbryos = T25CFlipped_FixedAPLengths >= T25CFlipped_FixedAPLimits(1) & T25CFlipped_FixedAPLengths <= T25CFlipped_FixedAPLimits(2);

T25CFlipped_CompiledEmbryos.EmbryoDimensionsGood = T25CFlipped_APKeepEmbryos & T25CFlipped_ARKeepEmbryos;


T25CFlipped_NC13Indices = T25CFlipped_CompiledEmbryos.EmbryoDimensionsGood & T25CFlipped_CompiledEmbryos.Flags == 0 &...
    T25CFlipped_CompiledEmbryos.Approved & T25CFlipped_CompiledEmbryos.nc == 13;
T25CFlipped_NC14Indices = T25CFlipped_CompiledEmbryos.EmbryoDimensionsGood & T25CFlipped_CompiledEmbryos.Flags == 0 &...
    T25CFlipped_CompiledEmbryos.Approved & T25CFlipped_CompiledEmbryos.nc == 14;


%%










