function CompiledEmbryos = RescaleSlideFluos(CompiledEmbryos, exp_index)
%%
APbins = 0:0.025:1;
NumAPbins = length(APbins);
AllSetInfo = GetFixedSetPrefixInfo;
NEmbryos = length(CompiledEmbryos.Approved);
SetIsFlipped = AllSetInfo.Flipped(exp_index);
NChannels = size(CompiledEmbryos.DorsalAvgAPProfiles, 3);
% Zero level APbins for knirps: 0.775-0.875, 0.375-0.45
% find peak between 0.525 and 0.7
% KnirpsBackground < 0.5 for T25C Flipped & slide_index = 1

AllEmbryos = 1:NEmbryos;
SlideIDs = unique(CompiledEmbryos.SlideIDs);
NSlides = length(SlideIDs);
ControlChannelBackgroundBins = zeros(NChannels, length(APbins), 'logical');
ControlChannelBackgroundBins(2,:) = APbins >= 0.1 & APbins <= 0.9;%APbins >= 0.375 & APbins <= 0.6; % Hoechst
ControlChannelBackgroundBins(3,:) = APbins >= 0.1 & APbins <= 0.9;%APbins >= 0.8; % Bicoid
ControlChannelBackgroundBins(4,:) = APbins >= 0.1 & APbins <= 0.9;%APbins >= 0.775; % Knirps
ControlChannelBackgroundBins(5,:) = APbins >= 0.1 & APbins <= 0.9;%APbins >= 0.55 & APbins <= 0.625; % Hunchback
CompiledEmbryos.ControlChannelBackgroundBins = ControlChannelBackgroundBins;
TestChannelBackgroundBins = ControlChannelBackgroundBins;
if ~SetIsFlipped
    CompiledEmbryos.ControlChannelBackgroundBins = ControlChannelBackgroundBins;
    CompiledEmbryos.TestChannelBackgroundBins = TestChannelBackgroundBins;
else
    CompiledEmbryos.ControlChannelBackgroundBins = TestChannelBackgroundBins;
    CompiledEmbryos.TestChannelBackgroundBins = ControlChannelBackgroundBins;
end
CompiledEmbryos.MeanTestBackgrounds = NaN(NSlides, NChannels);
CompiledEmbryos.CountTestBackgrounds = NaN(NSlides, NChannels);
CompiledEmbryos.SlideRescalingFactors = NaN(NSlides, NChannels);
CompiledEmbryos.SlideRescalingSEs = NaN(NSlides, NChannels);
CompiledEmbryos.SlideRescalingIntercepts = NaN(NSlides, NChannels);
CompiledEmbryos.SlideRescalingInterceptSEs = NaN(NSlides, NChannels);
FitInfo = cell(NSlides, NSlides, NChannels);
ScalingFactors = NaN(NSlides, NSlides, NChannels);
ScalingSEs = NaN(NSlides, NSlides, NChannels);
SlideRescalingIntercepts = NaN(NSlides, NSlides, NChannels);
SlideRescalingInterceptSEs = NaN(NSlides, NSlides, NChannels);
CompiledEmbryos.SlideRescalingFactors(1,:) = 1;
CompiledEmbryos.SlideRescalingSEs(1,:) = 0;
CompiledEmbryos.SlideRescalingIntercepts(1,:) = 0;
CompiledEmbryos.SlideRescalingInterceptSEs(1,:) = 0;

%%
AllDeltaValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
    CompiledEmbryos.TestSetEmbryos & ~isnan([CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].');% &...
AllDeltaValidProfilesControlTF = CompiledEmbryos.IsNC14 &...
    CompiledEmbryos.ControlSetEmbryos & ~isnan([CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].');% &...

x = 0:1:45;
NumTbins = length(x);

deltafc_binwidth = 2.5;
DorsalProfiles = CompiledEmbryos.DorsalAvgAPProfiles;
WindowedProfiles.DeltaFC = {};
WindowedProfiles.DeltaFC.BinWidth = deltafc_binwidth;
WindowedProfiles.DeltaFC.x = x;
WindowedProfiles.DeltaFC.Test = {};
WindowedProfiles.DeltaFC.Test.mean = NaN(NumTbins, NumAPbins, NChannels, NSlides);
WindowedProfiles.DeltaFC.Test.std = NaN(NumTbins, NumAPbins, NChannels, NSlides);
WindowedProfiles.DeltaFC.Test.se = NaN(NumTbins, NumAPbins, NChannels, NSlides);
WindowedProfiles.DeltaFC.Test.count = zeros(NumTbins, NumAPbins, NChannels, NSlides);
WindowedProfiles.DeltaFC.Control = {};
WindowedProfiles.DeltaFC.Control.mean = NaN(NumTbins, NumAPbins, NChannels, NSlides);
WindowedProfiles.DeltaFC.Control.std = NaN(NumTbins, NumAPbins, NChannels, NSlides);
WindowedProfiles.DeltaFC.Control.se = NaN(NumTbins, NumAPbins, NChannels, NSlides);
WindowedProfiles.DeltaFC.Control.count = zeros(NumTbins, NumAPbins, NChannels, NSlides);

for slide_index = 1:NSlides
    for i = 1:NumTbins
        TFtimeBin = CompiledEmbryos.SlideIDs == slide_index & [CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].' >= x(i)-deltafc_binwidth/2 & [CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].'< x(i)+deltafc_binwidth/2;
        TFtestBin = TFtimeBin & AllDeltaValidProfilesTestTF;
        TFcontrolBin = TFtimeBin & AllDeltaValidProfilesControlTF;
        for ch_index = 2:NChannels
            if sum(TFtestBin) >= 1
                for ap_index = 1:NumAPbins
                    TFGoodAP = ~isnan(DorsalProfiles(TFtestBin,ap_index,ch_index));
                    if sum(TFGoodAP) > 1
                        WindowedProfiles.DeltaFC.Test.mean(i,ap_index,ch_index, slide_index) = mean(DorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                        WindowedProfiles.DeltaFC.Test.std(i,ap_index,ch_index, slide_index) = std(DorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                        WindowedProfiles.DeltaFC.Test.count(i,ap_index,ch_index, slide_index) = sum(~isnan(DorsalProfiles(TFtestBin,ap_index,ch_index)));
                        WindowedProfiles.DeltaFC.Test.se(i,ap_index,ch_index, slide_index) = WindowedProfiles.DeltaFC.Test.std(i,ap_index,ch_index, slide_index)/...
                            sqrt(WindowedProfiles.DeltaFC.Test.count(i,ap_index,ch_index, slide_index));
                    elseif sum(TFGoodAP)  == 1
                        TestProfs = DorsalProfiles(TFtestBin,ap_index,ch_index);
                        CompiledEmbryos.WindowedProfiles.DeltaFC.Test.mean(i,ap_index,ch_index, slide_index) = TestProfs(TFGoodAP);
                        CompiledEmbryos.WindowedProfiles.DeltaFC.Test.std(i,ap_index,ch_index, slide_index) = NaN;
                        CompiledEmbryos.WindowedProfiles.DeltaFC.Test.count(i,ap_index,ch_index, slide_index) = 1;
                        CompiledEmbryos.WindowedProfiles.DeltaFC.Test.se(i,ap_index,ch_index, slide_index) = NaN;
                    end
                end
            end
            
            if sum(TFcontrolBin) >= 1
                for ap_index = 1:NumAPbins
                    TFGoodAP = ~isnan(DorsalProfiles(TFcontrolBin,ap_index,ch_index));
                    if sum(TFGoodAP) > 1
                        WindowedProfiles.DeltaFC.Control.mean(i,ap_index,ch_index, slide_index) = mean(DorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                        WindowedProfiles.DeltaFC.Control.std(i,ap_index,ch_index, slide_index) = std(DorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                        WindowedProfiles.DeltaFC.Control.count(i,ap_index,ch_index, slide_index) = sum(~isnan(DorsalProfiles(TFcontrolBin,ap_index,ch_index)));
                        WindowedProfiles.DeltaFC.Control.se(i,ap_index,ch_index, slide_index) = WindowedProfiles.DeltaFC.Control.std(i,ap_index,ch_index)/...
                            sqrt(WindowedProfiles.DeltaFC.Control.count(i,ap_index,ch_index, slide_index));
                    elseif sum(TFGoodAP)  == 1
                        ControlProfs = DorsalProfiles(TFcontrolBin,ap_index,ch_index);
                        WindowedProfiles.DeltaFC.Control.mean(i,ap_index,ch_index, slide_index) = ControlProfs(TFGoodAP);
                        WindowedProfiles.DeltaFC.Control.std(i,ap_index,ch_index, slide_index) = NaN;
                        WindowedProfiles.DeltaFC.Control.count(i,ap_index,ch_index, slide_index) = 1;
                        WindowedProfiles.DeltaFC.Control.se(i,ap_index,ch_index, slide_index) = NaN;
                    end
                end
            end
        end
    end
end
%%
AllSlides = 1:NSlides;
CountDims = [size(WindowedProfiles.DeltaFC.Control.count, 1),...
    size(WindowedProfiles.DeltaFC.Control.count, 2), NSlides];
counts = zeros([size(WindowedProfiles.DeltaFC.Control.count) NSlides]);

ControlledProfiles = NaN([size(WindowedProfiles.DeltaFC.Control.count) NSlides]);
for i = AllSlides
    counts(:,:,:,i) =WindowedProfiles.DeltaFC.Control.count(:,:,:, i);
    ControlledProfiles(:,:,:,i) = WindowedProfiles.DeltaFC.Control.mean(:,:,:, i);
end

FitInfo = cell(NSlides, NSlides, NChannels);
ScalingFactors = NaN(NSlides, NSlides, NChannels);
ScalingSEs = NaN(NSlides, NSlides, NChannels);
ScalingIntercepts = NaN(NSlides, NSlides, NChannels);
ScalingInterceptSEs =NaN(NSlides, NSlides, NChannels);
OverlapIndices = cell(NSlides, NSlides, NChannels);

for i = 1:(NSlides-1)
    for j = i+1:NSlides
        for ch_index = 2:NChannels
            
            OverlapIndices{j, i, ch_index} = (counts(:,ControlChannelBackgroundBins(ch_index, :),ch_index,i) >= 1) &   (counts(:,ControlChannelBackgroundBins(ch_index, :),ch_index,j) >= 1);
            
            
            if sum(sum((OverlapIndices{j, i, ch_index} )) >= 1) >= 5
                
                FitSet1 = ControlledProfiles(:,ControlChannelBackgroundBins(ch_index, :),ch_index,i);
                FitSet1 = FitSet1( OverlapIndices{j, i, ch_index}).';
                
                FitSet2 = ControlledProfiles(:,ControlChannelBackgroundBins(ch_index, :),ch_index,j);
                FitSet2 = FitSet2( OverlapIndices{j, i, ch_index}).';
                
                dlm = fitlm(FitSet2, FitSet1);%,'Weights',  1./(FitSD1.^2));
                FitInfo{j, i, ch_index} = dlm;
                ScalingFactors(j, i, ch_index) = dlm.Coefficients.Estimate(2);
                ScalingSEs(j, i, ch_index) = dlm.Coefficients.SE(2);
                ScalingIntercepts(j, i, ch_index) = dlm.Coefficients.Estimate(1);
                ScalingInterceptSEs(j, i, ch_index) = dlm.Coefficients.SE(1);
                
            end
        end
        
    end
end
%%
for ch_index = 2:NChannels
    for i = 1:(NSlides-1)
        for j = i+1:NSlides
            
            ScalingFactors(i, j, ch_index) = 1/ScalingFactors(j, i, ch_index);
            ScalingSEs(i, j, ch_index) = ScalingSEs(j, i, ch_index)/(ScalingFactors(j, i, ch_index)^2);
            ScalingIntercepts(i, j, ch_index) = -ScalingIntercepts(j, i, ch_index) /ScalingFactors(j, i, ch_index);
            ScalingInterceptSEs(i, j, ch_index) = sqrt(ScalingInterceptSEs(j, i, ch_index)^2+ScalingIntercepts(i, j, ch_index)^2*ScalingInterceptSEs(j, i, ch_index)^2)/ScalingFactors(j, i, ch_index);
            
        end
    end
end

%%
for ch_index = 2:NChannels
    for j = 2:NSlides
        CompiledEmbryos.SlideRescalingFactors(j, ch_index) = ScalingFactors(j, 1, ch_index);
        CompiledEmbryos.SlideRescalingSEs(j, ch_index) = ScalingSEs(j, 1, ch_index);
        CompiledEmbryos.SlideRescalingIntercepts(j, ch_index) = ScalingIntercepts(j, 1, ch_index);
        CompiledEmbryos.SlideRescalingInterceptSEs(j, ch_index) = ScalingInterceptSEs(j, 1, ch_index);
    end
    
end






%%
CompiledEmbryos.SlideRescaledDorsalAPProfiles = CompiledEmbryos.DorsalAPProfiles;
CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles = CompiledEmbryos.DorsalAvgAPProfiles;
for embryo_index = AllEmbryos
    for ch_index = 2:NChannels
        CompiledEmbryos.SlideRescaledDorsalAPProfiles(embryo_index,:,ch_index) =...
            CompiledEmbryos.SlideRescalingFactors(CompiledEmbryos.SlideIDs(embryo_index), ch_index)*CompiledEmbryos.SlideRescaledDorsalAPProfiles(embryo_index,:,ch_index)+...
            CompiledEmbryos.SlideRescalingIntercepts(CompiledEmbryos.SlideIDs(embryo_index), ch_index);
        CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(embryo_index,:,ch_index) =...
            CompiledEmbryos.SlideRescalingFactors(CompiledEmbryos.SlideIDs(embryo_index), ch_index)*CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(embryo_index,:,ch_index)+...
            CompiledEmbryos.SlideRescalingIntercepts(CompiledEmbryos.SlideIDs(embryo_index), ch_index);
    end
end

%%
NarrowAPbins = 0:0.0125:1;
NumNarrowAPbins = length(NarrowAPbins);
% NarrowControlChannelBackgroundBins = zeros(NChannels, length(NarrowAPbins), 'logical');
% NarrowControlChannelBackgroundBins(2,:) = NarrowAPbins >= 0.375 & NarrowAPbins <= 0.6; % Hoechst
% NarrowControlChannelBackgroundBins(3,:) = NarrowAPbins >= 0.8; % Bicoid
% NarrowControlChannelBackgroundBins(4,:) = NarrowAPbins >= 0.775; % Knirps
% NarrowControlChannelBackgroundBins(5,:) = NarrowAPbins >= 0.55 & NarrowAPbins <= 0.625; % Hunchback
NarrowControlChannelBackgroundBins = zeros(NChannels, length(NarrowAPbins), 'logical');
NarrowControlChannelBackgroundBins(2,:) = NarrowAPbins >= 0.1 & NarrowAPbins <= 0.9;%APbins >= 0.375 & APbins <= 0.6; % Hoechst
NarrowControlChannelBackgroundBins(3,:) = NarrowAPbins >= 0.1 & NarrowAPbins <= 0.9;%APbins >= 0.8; % Bicoid
NarrowControlChannelBackgroundBins(4,:) = NarrowAPbins >= 0.1 & NarrowAPbins <= 0.9;%APbins >= 0.775; % Knirps
NarrowControlChannelBackgroundBins(5,:) = NarrowAPbins >= 0.1 & NarrowAPbins <= 0.9;%APbins >= 0.55 & APbins <= 0.625; % Hunchback
CompiledEmbryos.NarrowControlChannelBackgroundBins = NarrowControlChannelBackgroundBins;
NarrowTestChannelBackgroundBins = NarrowControlChannelBackgroundBins;
if ~SetIsFlipped
    CompiledEmbryos.NarrowControlChannelBackgroundBins = NarrowControlChannelBackgroundBins;
    CompiledEmbryos.NarrowTestChannelBackgroundBins = NarrowTestChannelBackgroundBins;
else
    CompiledEmbryos.NarrowControlChannelBackgroundBins = NarrowTestChannelBackgroundBins;
    CompiledEmbryos.NarrowTestChannelBackgroundBins = NarrowControlChannelBackgroundBins;
end

CompiledEmbryos.NarrowMeanTestBackgrounds = NaN(NSlides, NChannels);
CompiledEmbryos.NarrowCountTestBackgrounds = NaN(NSlides, NChannels);
CompiledEmbryos.NarrowSlideRescalingFactors = NaN(NSlides, NChannels);
CompiledEmbryos.NarrowSlideRescalingSEs = NaN(NSlides, NChannels);
CompiledEmbryos.NarrowSlideRescalingIntercepts = NaN(NSlides, NChannels);
CompiledEmbryos.NarrowSlideRescalingInterceptSEs = NaN(NSlides, NChannels);
FitInfo = cell(NSlides, NSlides, NChannels);
ScalingFactors = NaN(NSlides, NSlides, NChannels);
ScalingSEs = NaN(NSlides, NSlides, NChannels);
CompiledEmbryos.NarrowSlideRescalingFactors(1,:) = 1;
CompiledEmbryos.NarrowSlideRescalingSEs(1,:) = 0;
CompiledEmbryos.NarrowSlideRescalingIntercepts(1,:) = 1;
CompiledEmbryos.NarrowSlideRescalingInterceptSEs(1,:) = 0;
NarrowDorsalProfiles = CompiledEmbryos.DorsalAvgNarrowAPProfiles;
WindowedProfiles.DeltaFC = {};
WindowedProfiles.DeltaFC.BinWidth = deltafc_binwidth;
WindowedProfiles.DeltaFC.x = x;
WindowedProfiles.DeltaFC.Test = {};
WindowedProfiles.DeltaFC.Test.mean = NaN(NumTbins, NumNarrowAPbins, NChannels, NSlides);
WindowedProfiles.DeltaFC.Test.std = NaN(NumTbins, NumNarrowAPbins, NChannels, NSlides);
WindowedProfiles.DeltaFC.Test.se = NaN(NumTbins, NumNarrowAPbins, NChannels, NSlides);
WindowedProfiles.DeltaFC.Test.count = zeros(NumTbins, NumNarrowAPbins, NChannels, NSlides);
WindowedProfiles.DeltaFC.Control = {};
WindowedProfiles.DeltaFC.Control.mean = NaN(NumTbins, NumNarrowAPbins, NChannels, NSlides);
WindowedProfiles.DeltaFC.Control.std = NaN(NumTbins, NumNarrowAPbins, NChannels, NSlides);
WindowedProfiles.DeltaFC.Control.se = NaN(NumTbins, NumNarrowAPbins, NChannels, NSlides);
WindowedProfiles.DeltaFC.Control.count = zeros(NumTbins, NumNarrowAPbins, NChannels, NSlides);

for slide_index = 1:NSlides
    for i = 1:NumTbins
        TFtimeBin = CompiledEmbryos.SlideIDs == slide_index & [CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].' >= x(i)-deltafc_binwidth/2 & [CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].'< x(i)+deltafc_binwidth/2;
        TFtestBin = TFtimeBin & AllDeltaValidProfilesTestTF;
        TFcontrolBin = TFtimeBin & AllDeltaValidProfilesControlTF;
        for ch_index = 2:NChannels
            if sum(TFtestBin) >= 1
                for ap_index = 1:NumNarrowAPbins
                    TFGoodAP = ~isnan(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index));
                    if sum(TFGoodAP) > 1
                        WindowedProfiles.DeltaFC.Test.mean(i,ap_index,ch_index, slide_index) = mean(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                        WindowedProfiles.DeltaFC.Test.std(i,ap_index,ch_index, slide_index) = std(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index), 'omitnan');
                        WindowedProfiles.DeltaFC.Test.count(i,ap_index,ch_index, slide_index) = sum(~isnan(NarrowDorsalProfiles(TFtestBin,ap_index,ch_index)));
                        WindowedProfiles.DeltaFC.Test.se(i,ap_index,ch_index, slide_index) = WindowedProfiles.DeltaFC.Test.std(i,ap_index,ch_index, slide_index)/...
                            sqrt(WindowedProfiles.DeltaFC.Test.count(i,ap_index,ch_index, slide_index));
                    elseif sum(TFGoodAP)  == 1
                        TestProfs = NarrowDorsalProfiles(TFtestBin,ap_index,ch_index);
                        CompiledEmbryos.WindowedProfiles.DeltaFC.Test.mean(i,ap_index,ch_index, slide_index) = TestProfs(TFGoodAP);
                        CompiledEmbryos.WindowedProfiles.DeltaFC.Test.std(i,ap_index,ch_index, slide_index) = NaN;
                        CompiledEmbryos.WindowedProfiles.DeltaFC.Test.count(i,ap_index,ch_index, slide_index) = 1;
                        CompiledEmbryos.WindowedProfiles.DeltaFC.Test.se(i,ap_index,ch_index, slide_index) = NaN;
                    end
                end
            end
            
            if sum(TFcontrolBin) >= 1
                for ap_index = 1:NumAPbins
                    TFGoodAP = ~isnan(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index));
                    if sum(TFGoodAP) > 1
                        WindowedProfiles.DeltaFC.Control.mean(i,ap_index,ch_index, slide_index) = mean(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                        WindowedProfiles.DeltaFC.Control.std(i,ap_index,ch_index, slide_index) = std(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index), 'omitnan');
                        WindowedProfiles.DeltaFC.Control.count(i,ap_index,ch_index, slide_index) = sum(~isnan(NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index)));
                        WindowedProfiles.DeltaFC.Control.se(i,ap_index,ch_index, slide_index) = WindowedProfiles.DeltaFC.Control.std(i,ap_index,ch_index)/...
                            sqrt(WindowedProfiles.DeltaFC.Control.count(i,ap_index,ch_index, slide_index));
                    elseif sum(TFGoodAP)  == 1
                        ControlProfs = NarrowDorsalProfiles(TFcontrolBin,ap_index,ch_index);
                        WindowedProfiles.DeltaFC.Control.mean(i,ap_index,ch_index, slide_index) = ControlProfs(TFGoodAP);
                        WindowedProfiles.DeltaFC.Control.std(i,ap_index,ch_index, slide_index) = NaN;
                        WindowedProfiles.DeltaFC.Control.count(i,ap_index,ch_index, slide_index) = 1;
                        WindowedProfiles.DeltaFC.Control.se(i,ap_index,ch_index, slide_index) = NaN;
                    end
                end
            end
        end
    end
end
%%
AllSlides = 1:NSlides;
CountDims = [size(WindowedProfiles.DeltaFC.Control.count, 1),...
    size(WindowedProfiles.DeltaFC.Control.count, 2), NSlides];
counts = zeros([size(WindowedProfiles.DeltaFC.Control.count) NSlides]);

ControlledProfiles = NaN([size(WindowedProfiles.DeltaFC.Control.count)]);
for i = AllSlides
    counts(:,:,:,i) =WindowedProfiles.DeltaFC.Control.count(:,:,:, i);
    ControlledProfiles(:,:,:,i) = WindowedProfiles.DeltaFC.Control.mean(:,:,:, i);
end

FitInfo = cell(NSlides, NSlides, NChannels);
ScalingFactors = NaN(NSlides, NSlides, NChannels);
ScalingSEs = NaN(NSlides, NSlides, NChannels);
ScalingIntercepts = NaN(NSlides, NSlides, NChannels);
ScalingInterceptSEs =NaN(NSlides, NSlides, NChannels);
OverlapIndices = cell(NSlides, NSlides, NChannels);

for i = 1:(NSlides-1)
    for j = i+1:NSlides
        for ch_index = 2:NChannels
            
            OverlapIndices{j, i, ch_index} = (counts(:,NarrowControlChannelBackgroundBins(ch_index, :),ch_index,i) >= 1) &   (counts(:,NarrowControlChannelBackgroundBins(ch_index, :),ch_index,j) >= 1);
            
            
            if sum(sum((OverlapIndices{j, i, ch_index} )) >= 1) >= 5
                
                FitSet1 = ControlledProfiles(:,NarrowControlChannelBackgroundBins(ch_index, :),ch_index,i);
                FitSet1 = FitSet1( OverlapIndices{j, i, ch_index}).';
                
                FitSet2 = ControlledProfiles(:,NarrowControlChannelBackgroundBins(ch_index, :),ch_index,j);
                FitSet2 = FitSet2( OverlapIndices{j, i, ch_index}).';
                
                dlm = fitlm(FitSet2, FitSet1);%,'Weights',  1./(FitSD1.^2));
                FitInfo{j, i, ch_index} = dlm;
                ScalingFactors(j, i, ch_index) = dlm.Coefficients.Estimate(2);
                ScalingSEs(j, i, ch_index) = dlm.Coefficients.SE(2);
                ScalingIntercepts(j, i, ch_index) = dlm.Coefficients.Estimate(1);
                ScalingInterceptSEs(j, i, ch_index) = dlm.Coefficients.SE(1);
                
            end
        end
        
    end
end
%%
for ch_index = 2:NChannels
    for i = 1:(NSlides-1)
        for j = i+1:NSlides
            
            ScalingFactors(i, j, ch_index) = 1/ScalingFactors(j, i, ch_index);
            ScalingSEs(i, j, ch_index) = ScalingSEs(j, i, ch_index)/(ScalingFactors(j, i, ch_index)^2);
            ScalingIntercepts(i, j, ch_index) = -ScalingIntercepts(j, i, ch_index) /ScalingFactors(j, i, ch_index);
            ScalingInterceptSEs(i, j, ch_index) = sqrt(ScalingInterceptSEs(j, i, ch_index)^2+ScalingIntercepts(i, j, ch_index)^2*ScalingInterceptSEs(j, i, ch_index)^2)/ScalingFactors(j, i, ch_index);
            
        end
    end
end

%%
for ch_index = 2:NChannels
    for j = 2:NSlides
        CompiledEmbryos.NarrowSlideRescalingFactors(j, ch_index) = ScalingFactors(j, 1, ch_index);
        CompiledEmbryos.NarrowSlideRescalingSEs(j, ch_index) = ScalingSEs(j, 1, ch_index);
        CompiledEmbryos.NarrowSlideRescalingIntercepts(j, ch_index) = ScalingIntercepts(j, 1, ch_index);
        CompiledEmbryos.NarrowSlideRescalingInterceptSEs(j, ch_index) = ScalingInterceptSEs(j, 1, ch_index);
    end
    
end






%%
CompiledEmbryos.SlideRescaledDorsalNarrowAPProfiles = CompiledEmbryos.DorsalNarrowAPProfiles;
CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles = CompiledEmbryos.DorsalAvgNarrowAPProfiles;
for embryo_index = AllEmbryos
    for ch_index = 2:NChannels
        CompiledEmbryos.SlideRescaledDorsalNarrowAPProfiles(embryo_index,:,ch_index) =...
            CompiledEmbryos.NarrowSlideRescalingFactors(CompiledEmbryos.SlideIDs(embryo_index), ch_index)*CompiledEmbryos.SlideRescaledDorsalNarrowAPProfiles(embryo_index,:,ch_index)+...
            CompiledEmbryos.NarrowSlideRescalingIntercepts(CompiledEmbryos.SlideIDs(embryo_index), ch_index);
        CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(embryo_index,:,ch_index) =...
            CompiledEmbryos.NarrowSlideRescalingFactors(CompiledEmbryos.SlideIDs(embryo_index), ch_index)*CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(embryo_index,:,ch_index)+...
            CompiledEmbryos.NarrowSlideRescalingIntercepts(CompiledEmbryos.SlideIDs(embryo_index), ch_index);
    end
end



%%
CompiledEmbryos.SlideRescaledDorsalNarrowAPProfiles = CompiledEmbryos.DorsalNarrowAPProfiles;
CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles = CompiledEmbryos.DorsalAvgNarrowAPProfiles;
for embryo_index = AllEmbryos
    for ch_index = 2:NChannels
        CompiledEmbryos.SlideRescaledDorsalNarrowAPProfiles(embryo_index,:,ch_index) =...
            CompiledEmbryos.NarrowSlideRescalingFactors(CompiledEmbryos.SlideIDs(embryo_index), ch_index)*CompiledEmbryos.SlideRescaledDorsalNarrowAPProfiles(embryo_index,:,ch_index)+...
            CompiledEmbryos.NarrowSlideRescalingIntercepts(CompiledEmbryos.SlideIDs(embryo_index), ch_index);
        CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(embryo_index,:,ch_index) =...
            CompiledEmbryos.NarrowSlideRescalingFactors(CompiledEmbryos.SlideIDs(embryo_index), ch_index)*CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(embryo_index,:,ch_index)+...
            CompiledEmbryos.NarrowSlideRescalingIntercepts(CompiledEmbryos.SlideIDs(embryo_index), ch_index);
    end
end



