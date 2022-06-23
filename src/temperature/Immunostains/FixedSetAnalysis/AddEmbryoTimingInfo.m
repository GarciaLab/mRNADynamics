function CompiledEmbryos = AddEmbryoTimingInfo(CompiledEmbryos, exp_index)
%% 
AllSetInfo = GetFixedSetPrefixInfo;
SetPrefixes = AllSetInfo.Prefixes{exp_index};
SizeDataPath = 'S:/Gabriella/Dropbox/EmbryoSizeMeasurements/EmbryoSizeData.mat';
AllSetsProfFigPath = 'S:/Gabriella/Dropbox/ProteinProfiles/Figures/';
AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';
load(SizeDataPath,'APLengths','DVLengths','VentralDistances','DorsalDistances','NGoodEmbryos',...
    'MeanAPLength','StdAPLength','MeanDVLength','StdDVLength','AspectRatios');
MeanAspectRatio = mean(AspectRatios);
StdAspectRatio = std(AspectRatios);
[DubuisTimes, DubuisMeanProfile, ~] = getMembraneFurrowProfiles( 'dubuis');
[hisrfp25CTimes, hisrfp25CProfile, ~] = getMembraneFurrowProfiles( 'hisrfp25c_nopv');
[yw25CTimes, yw25CProfile, ~] = getMembraneFurrowProfiles( 'yw25csquished_nopv');
liveExperiments = cell(1, length(SetPrefixes));
for i = 1:length(SetPrefixes)
    liveExperiments{i} = LiveExperiment(SetPrefixes{i});
end
FixedPixelSize_um = liveExperiments{1}.pixelSize_um;

ApprovedEmbryos = CompiledEmbryos.Approved;
NEmbryos = length(ApprovedEmbryos);
%%

FixedAPLengths = NaN(1, NEmbryos);
FixedAPSlopes =  NaN(1, NEmbryos);
FixedAPIntercepts =  NaN(1, NEmbryos);
FixedDorsalDistances =  NaN(1, NEmbryos);
FixedVentralDistances =  NaN(1, NEmbryos);
FixedDVLengths =  NaN(1, NEmbryos);

for i = 1:NEmbryos
    if CompiledEmbryos.Approved(i) & ((CompiledEmbryos.Flags(i) == 0 & CompiledEmbryos.nc(i) ~= 13) | ((CompiledEmbryos.Flags(i) == 0 |CompiledEmbryos.Flags(i) == 5) & CompiledEmbryos.nc(i) == 13)) 
        FixedAPLengths(i) = FixedPixelSize_um*sqrt((CompiledEmbryos.CoordAs(i,1)-CompiledEmbryos.CoordPs(i,1))^2+...
            (CompiledEmbryos.CoordAs(i,2)-CompiledEmbryos.CoordPs(i,2))^2);
        FixedAPSlopes(i) = (CompiledEmbryos.CoordAs(i,2)-CompiledEmbryos.CoordPs(i,2))/(CompiledEmbryos.CoordAs(i,1)-CompiledEmbryos.CoordPs(i,1));
        FixedAPIntercepts(i) = CompiledEmbryos.CoordAs(i,2)-FixedAPSlopes(i)*CompiledEmbryos.CoordAs(i,1);
        FixedDorsalDistances(i) = FixedPixelSize_um*abs(FixedAPIntercepts(i)+FixedAPSlopes(i)*CompiledEmbryos.CoordDs(i,1)-CompiledEmbryos.CoordDs(i,2))/sqrt(FixedAPSlopes(i)^2+1);
        FixedVentralDistances(i) = FixedPixelSize_um*abs(FixedAPIntercepts(i)+FixedAPSlopes(i)*CompiledEmbryos.CoordVs(i,1)-CompiledEmbryos.CoordVs(i,2))/sqrt(FixedAPSlopes(i)^2+1);
        FixedDVLengths(i) = FixedDorsalDistances(i)+FixedVentralDistances(i);
    end
end

CompiledEmbryos.APLengths = FixedAPLengths;
CompiledEmbryos.APSlopes = FixedAPSlopes;
CompiledEmbryos.APIntercepts = FixedAPIntercepts;
CompiledEmbryos.DorsalDistances = FixedDorsalDistances;
CompiledEmbryos.VentralDistances = FixedVentralDistances;
CompiledEmbryos.DVLengths = FixedDVLengths;

FixedAspectRatios = FixedDVLengths./FixedAPLengths;
CompiledEmbryos.AspectRatios = FixedAspectRatios;
MeanFixedAspectRatio = mean(FixedAspectRatios, 'omitnan');
StdFixedAspectRatio = std(FixedAspectRatios, 'omitnan');

MeanFixedAPLength = mean(FixedAPLengths, 'omitnan');
StdFixedAPLength = std(FixedAPLengths, 'omitnan');
MeanFixedDVLength = mean(FixedDVLengths, 'omitnan');
StdFixedDVLength = std(FixedDVLengths, 'omitnan');


FixedAPLengthFit = fitdist(max(FixedAPLengths)-FixedAPLengths.', 'logistic');
FixedAPDistCenter = -1*FixedAPLengthFit.ParameterValues(1)+max(FixedAPLengths);
FixedAPLimits = [FixedAPDistCenter+(min(APLengths)-MeanAPLength)*FixedAPDistCenter/MeanAPLength,...
    FixedAPDistCenter+(max(APLengths)-MeanAPLength)*FixedAPDistCenter/MeanAPLength];
%
ARLimits = [min(AspectRatios,[], 'omitnan'), max(AspectRatios,[], 'omitnan')];
ARKeepEmbryos = FixedAspectRatios >= ARLimits(1) & FixedAspectRatios <= ARLimits(2);

APKeepEmbryos = FixedAPLengths >= FixedAPLimits(1) & FixedAPLengths <= FixedAPLimits(2);

CompiledEmbryos.EmbryoDimensionsGood = APKeepEmbryos & ARKeepEmbryos;

IsNC13 = CompiledEmbryos.EmbryoDimensionsGood & (CompiledEmbryos.Flags == 0 | CompiledEmbryos.Flags == 5) &...
    CompiledEmbryos.Approved & CompiledEmbryos.nc == 13;
NC13Indices = find(IsNC13);
IsNC14 = CompiledEmbryos.EmbryoDimensionsGood & (CompiledEmbryos.Flags == 0) &...
    CompiledEmbryos.Approved & CompiledEmbryos.nc == 14;
NC14Indices = find(IsNC14);

IsNC13orNC14 = IsNC13 | IsNC14;
NC13NC14Indices = find(IsNC13orNC14);
CompiledEmbryos.IsNC14 = IsNC14;
CompiledEmbryos.IsNC13 = IsNC13;
CompiledEmbryos.IsNC13orNC14 = IsNC13orNC14;

CompiledEmbryos.NumEmbryosNC13 = length(NC13Indices);
CompiledEmbryos.NumEmbryosNC14 = length(NC14Indices);
CompiledEmbryos.NumEmbryosNC13NC14 = length(NC13NC14Indices);
%%
%CompiledEmbryos.

DimensionsGoodMeanFixedAPLength = mean(FixedAPLengths(CompiledEmbryos.EmbryoDimensionsGood));
DimensionsGoodStdFixedAPLength = std(FixedAPLengths(CompiledEmbryos.EmbryoDimensionsGood));
DimensionsGoodSEFixedAPLength = DimensionsGoodStdFixedAPLength/sqrt(sum(CompiledEmbryos.EmbryoDimensionsGood));
DimensionsGoodMeanFixedDVLength = mean(FixedDVLengths(CompiledEmbryos.EmbryoDimensionsGood));
DimensionsGoodStdFixedDVLength = std(FixedDVLengths(CompiledEmbryos.EmbryoDimensionsGood));
DimensionsGoodSEFixedDVLength = DimensionsGoodStdFixedDVLength/sqrt(sum(CompiledEmbryos.EmbryoDimensionsGood));
SEAPLength = StdAPLength/sqrt(length(APLengths));
SEDVLength = StdDVLength/sqrt(length(DVLengths));
MeanAPRescalingFactor = MeanAPLength/DimensionsGoodMeanFixedAPLength;
MeanAPRescalingFactorSE = 1/DimensionsGoodMeanFixedAPLength*sqrt(SEAPLength^2+DimensionsGoodSEFixedAPLength^2*MeanAPRescalingFactor^2);
MeanDVRescalingFactor = MeanDVLength/DimensionsGoodMeanFixedDVLength;
MeanDVRescalingFactorSE = 1/DimensionsGoodMeanFixedDVLength*sqrt(SEDVLength^2+DimensionsGoodSEFixedDVLength^2*MeanDVRescalingFactor^2);
MeanRescalingFactor = (MeanAPRescalingFactor+MeanDVRescalingFactor)/2;
MeanRescalingFactorSE = sqrt(MeanAPRescalingFactorSE^2+MeanDVRescalingFactorSE^2)/2;

CompiledEmbryos.FixCorrectedDeltaFC_um = {};
CompiledEmbryos.FixCorrectedDeltaFC_um.mean = NaN(1, NEmbryos);
CompiledEmbryos.FixCorrectedDeltaFC_um.std = NaN(1, NEmbryos);
CompiledEmbryos.FixCorrectedDeltaFC_um.count = NaN(1, NEmbryos);
CompiledEmbryos.FixCorrectedDeltaFC_um.se = NaN(1, NEmbryos);
CompiledEmbryos.FixMeanCorrectedDeltaFC_um = {};
CompiledEmbryos.FixMeanCorrectedDeltaFC_um.mean = NaN(1, NEmbryos);
CompiledEmbryos.FixMeanCorrectedDeltaFC_um.std = NaN(1, NEmbryos);
CompiledEmbryos.FixMeanCorrectedDeltaFC_um.count = NaN(1, NEmbryos);
CompiledEmbryos.FixMeanCorrectedDeltaFC_um.se = NaN(1, NEmbryos);
CompiledEmbryos.EmbryoSizeScaleFactor = NaN(NEmbryos,2);
CompiledEmbryos.EmbryoAPSizeScaleFactor = NaN(NEmbryos,2);
CompiledEmbryos.EmbryoDVSizeScaleFactor = NaN(NEmbryos,2);

for embryo_index = 1:NEmbryos
    if ~isnan(CompiledEmbryos.APLengths(embryo_index)) & ~isnan(CompiledEmbryos.DeltaFC_um.mean(embryo_index))
        CompiledEmbryos.EmbryoAPSizeScaleFactor(embryo_index, 1) = MeanAPLength/CompiledEmbryos.APLengths(embryo_index);
        CompiledEmbryos.EmbryoAPSizeScaleFactor(embryo_index, 2) = SEAPLength/CompiledEmbryos.APLengths(embryo_index);
        CompiledEmbryos.EmbryoDVSizeScaleFactor(embryo_index, 1) = MeanDVLength/CompiledEmbryos.DVLengths(embryo_index);
        CompiledEmbryos.EmbryoDVSizeScaleFactor(embryo_index, 2) = SEDVLength/CompiledEmbryos.DVLengths(embryo_index);
        CompiledEmbryos.EmbryoSizeScaleFactor(embryo_index, 1) = (CompiledEmbryos.EmbryoAPSizeScaleFactor(embryo_index, 1) +CompiledEmbryos.EmbryoDVSizeScaleFactor(embryo_index, 1))/2;
        CompiledEmbryos.EmbryoSizeScaleFactor(embryo_index, 2) = sqrt(CompiledEmbryos.EmbryoAPSizeScaleFactor(embryo_index, 2)^2 +CompiledEmbryos.EmbryoDVSizeScaleFactor(embryo_index, 2)^2)/2;
        CompiledEmbryos.FixCorrectedDeltaFC_um.mean(embryo_index) = CompiledEmbryos.EmbryoSizeScaleFactor(embryo_index, 1)*CompiledEmbryos.DeltaFC_um.mean(embryo_index);
        CompiledEmbryos.FixCorrectedDeltaFC_um.std(embryo_index) = CompiledEmbryos.EmbryoSizeScaleFactor(embryo_index, 1)*CompiledEmbryos.DeltaFC_um.std(embryo_index);
        CompiledEmbryos.FixCorrectedDeltaFC_um.count(embryo_index) = CompiledEmbryos.DeltaFC_um.count(embryo_index);
        CompiledEmbryos.FixCorrectedDeltaFC_um.se(embryo_index) = sqrt(CompiledEmbryos.EmbryoSizeScaleFactor(embryo_index, 1)^2*CompiledEmbryos.DeltaFC_um.se(embryo_index)^2 +...
            CompiledEmbryos.EmbryoSizeScaleFactor(embryo_index, 2)^2*CompiledEmbryos.DeltaFC_um.mean(embryo_index)^2); 
        CompiledEmbryos.FixMeanCorrectedDeltaFC_um.mean(embryo_index) = MeanRescalingFactor*CompiledEmbryos.DeltaFC_um.mean(embryo_index);
        CompiledEmbryos.FixMeanCorrectedDeltaFC_um.std(embryo_index) = MeanRescalingFactor*CompiledEmbryos.DeltaFC_um.std(embryo_index);
        CompiledEmbryos.FixMeanCorrectedDeltaFC_um.count(embryo_index) = CompiledEmbryos.DeltaFC_um.count(embryo_index);
        CompiledEmbryos.FixMeanCorrectedDeltaFC_um.se(embryo_index) = sqrt(MeanRescalingFactor^2*CompiledEmbryos.DeltaFC_um.se(embryo_index)^2 ...
            + MeanRescalingFactorSE^2*CompiledEmbryos.DeltaFC_um.mean(embryo_index)^2);
    end
end
%% Add Dubuis Profile Times
CompiledEmbryos.DubuisEmbryoTimeBins = NaN(1, NEmbryos);
CompiledEmbryos.DubuisEmbryoTimes = NaN(1, NEmbryos);
for embryo_index=1:NEmbryos
    if ~isnan(CompiledEmbryos.FixCorrectedDeltaFC_um.mean(embryo_index))
        EmbryoDeltaFC = CompiledEmbryos.FixCorrectedDeltaFC_um.mean(embryo_index);
        BinNumber = find(DubuisMeanProfile(1:end-1) <= EmbryoDeltaFC & DubuisMeanProfile(2:end) > EmbryoDeltaFC);
        if ~isempty(BinNumber)
            CompiledEmbryos.DubuisEmbryoTimeBins(embryo_index) = BinNumber;
            
            CompiledEmbryos.DubuisEmbryoTimes(embryo_index) = DubuisTimes(BinNumber)+(EmbryoDeltaFC-DubuisMeanProfile(BinNumber))*(DubuisTimes(BinNumber+1)-DubuisTimes(BinNumber))/(DubuisMeanProfile(BinNumber+1)-DubuisMeanProfile(BinNumber));
        end
        
    end
end


%% Add hisrfp 25C Profile Times
CompiledEmbryos.HisRFP25CEmbryoTimeBins = NaN(1, NEmbryos);
CompiledEmbryos.HisRFP25CEmbryoTimes = NaN(1, NEmbryos);
for embryo_index=1:NEmbryos
    if ~isnan(CompiledEmbryos.FixCorrectedDeltaFC_um.mean(embryo_index))
        EmbryoDeltaFC = CompiledEmbryos.FixCorrectedDeltaFC_um.mean(embryo_index);
        BinNumber = find(hisrfp25CProfile(1:end-1) <= EmbryoDeltaFC & hisrfp25CProfile(2:end) > EmbryoDeltaFC);
        if ~isempty(BinNumber)
            CompiledEmbryos.HisRFP25CEmbryoTimeBins(embryo_index) = BinNumber;
            
            CompiledEmbryos.HisRFP25CEmbryoTimes(embryo_index) = hisrfp25CTimes(BinNumber)+(EmbryoDeltaFC-hisrfp25CProfile(BinNumber))*(hisrfp25CTimes(BinNumber+1)-hisrfp25CTimes(BinNumber))/(hisrfp25CProfile(BinNumber+1)-hisrfp25CProfile(BinNumber));
        end
        
    end
end

%% Add yw 25C Profile Times
CompiledEmbryos.yw25CEmbryoTimeBins = NaN(1, NEmbryos);
CompiledEmbryos.yw25CEmbryoTimes = NaN(1, NEmbryos);
for embryo_index=1:NEmbryos
    if ~isnan(CompiledEmbryos.FixCorrectedDeltaFC_um.mean(embryo_index))
        EmbryoDeltaFC = CompiledEmbryos.FixCorrectedDeltaFC_um.mean(embryo_index);
        BinNumber = find(yw25CProfile(1:end-1) <= EmbryoDeltaFC & yw25CProfile(2:end) > EmbryoDeltaFC);
        if ~isempty(BinNumber)
            CompiledEmbryos.yw25CEmbryoTimeBins(embryo_index) = BinNumber;
            
            CompiledEmbryos.yw25CEmbryoTimes(embryo_index) = yw25CTimes(BinNumber)+(EmbryoDeltaFC-yw25CProfile(BinNumber))*(yw25CTimes(BinNumber+1)-yw25CTimes(BinNumber))/(yw25CProfile(BinNumber+1)-yw25CProfile(BinNumber));
        end  
    end
end

%% Add yw 25C Profile Times
CompiledEmbryos.TsetEmbryoTimeBins = NaN(1, NEmbryos);
CompiledEmbryos.TsetEmbryoTimes = NaN(1, NEmbryos);



