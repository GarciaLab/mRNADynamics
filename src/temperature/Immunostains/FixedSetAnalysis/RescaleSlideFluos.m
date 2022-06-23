function CompiledEmbryos = RescaleSlideFluos(CompiledEmbryos, exp_index)
%%
APbins = 0:0.025:1;
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
ControlChannelBackgroundBins(2,:) = APbins >= 0.375 & APbins <= 0.6; % Hoechst
ControlChannelBackgroundBins(3,:) = APbins >= 0.8; % Bicoid
ControlChannelBackgroundBins(4,:) = APbins >= 0.775; % Knirps
ControlChannelBackgroundBins(5,:) = APbins >= 0.55 & APbins <= 0.625; % Hunchback
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
FitInfo = cell(NSlides, NSlides, NChannels);
ScalingFactors = NaN(NSlides, NSlides, NChannels);
ScalingSEs = NaN(NSlides, NSlides, NChannels);
CompiledEmbryos.SlideRescalingFactors(1,:) = 1;
CompiledEmbryos.SlideRescalingSEs(1,:) = 0;

for ch_index = 2:NChannels
    if ~SetIsFlipped
        MeanChannelDorsalProfilesAllSlides = NaN( NSlides, sum(TestChannelBackgroundBins(ch_index,:)));
        StdChannelDorsalProfilesAllSlides = NaN( NSlides, sum(TestChannelBackgroundBins(ch_index,:)));
    else
        MeanChannelDorsalProfilesAllSlides = NaN( NSlides, sum(ControlChannelBackgroundBins(ch_index,:)));
        StdChannelDorsalProfilesAllSlides = NaN( NSlides, sum(ControlChannelBackgroundBins(ch_index,:)));
    end
    
    CountsDorsalProfilesAllSlides  = NaN(1, NSlides);
    for slide_index = 1:NSlides
        if ~SetIsFlipped
            ProfilesToUse = CompiledEmbryos.SlideTestSetEmbryos(slide_index,:) & CompiledEmbryos.IsNC14;
            DorsalProfiles = CompiledEmbryos.DorsalAvgAPProfiles(ProfilesToUse,TestChannelBackgroundBins(ch_index,:),ch_index);
        else
            ProfilesToUse = CompiledEmbryos.SlideControlSetEmbryos(slide_index,:) & CompiledEmbryos.IsNC14;
            DorsalProfiles = CompiledEmbryos.DorsalAvgAPProfiles(ProfilesToUse,ControlChannelBackgroundBins(ch_index,:),ch_index);
        end
        MeanChannelDorsalProfilesAllSlides(slide_index,:) =  mean(DorsalProfiles,1, 'omitnan');
        StdChannelDorsalProfilesAllSlides(slide_index,:) =  std(DorsalProfiles,1, 'omitnan');
        CountsDorsalProfilesAllSlides(slide_index) = size(DorsalProfiles, 1);
        CompiledEmbryos.MeanTestBackgrounds(slide_index, ch_index) = mean(MeanChannelDorsalProfilesAllSlides(slide_index,:), 'omitnan');
        CompiledEmbryos.CountTestBackgrounds(slide_index, ch_index) = CountsDorsalProfilesAllSlides(slide_index);
    end
    
    
    for i = 1:(NSlides-1)
        for j = i+1:NSlides
            FitSet1 = MeanChannelDorsalProfilesAllSlides(i,:);
            FitSet1 = FitSet1(:);
            FitSet2 = MeanChannelDorsalProfilesAllSlides(j,:);
            FitSet2 = FitSet2(:);
            dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
            FitInfo{j, i, ch_index} = dlm;
            ScalingFactors(j, i, ch_index) = dlm.Coefficients.Estimate;
            ScalingSEs(j, i, ch_index) = dlm.Coefficients.SE;
            
            
        end
    end
    
    for i = 1:(NSlides-1)
        for j = i+1:NSlides
            
            ScalingFactors(i, j, ch_index) = 1/ScalingFactors(j, i, ch_index);
            ScalingSEs(i, j, ch_index) = ScalingSEs(j, i, ch_index)/ScalingFactors(j, i, ch_index);
            
        end
    end
    
    for j = 2:NSlides
        CompiledEmbryos.SlideRescalingFactors(j, ch_index) = ScalingFactors(j, 1, ch_index);
        CompiledEmbryos.SlideRescalingSEs(j, ch_index) = ScalingSEs(j, 1, ch_index);
    end
    
    
    
end




%%
CompiledEmbryos.SlideRescaledDorsalAPProfiles = CompiledEmbryos.DorsalAPProfiles;
CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles = CompiledEmbryos.DorsalAvgAPProfiles;
for embryo_index = AllEmbryos
    for ch_index = 2:NChannels
        CompiledEmbryos.SlideRescaledDorsalAPProfiles(embryo_index,:,ch_index) =...
            CompiledEmbryos.SlideRescalingFactors(CompiledEmbryos.SlideIDs(embryo_index), ch_index)*CompiledEmbryos.SlideRescaledDorsalAPProfiles(embryo_index,:,ch_index);
        CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(embryo_index,:,ch_index) =...
            CompiledEmbryos.SlideRescalingFactors(CompiledEmbryos.SlideIDs(embryo_index), ch_index)*CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(embryo_index,:,ch_index);
    end
end

%%
NarrowAPbins = 0:0.0125:1;
NarrowControlChannelBackgroundBins = zeros(NChannels, length(NarrowAPbins), 'logical');
NarrowControlChannelBackgroundBins(2,:) = NarrowAPbins >= 0.375 & NarrowAPbins <= 0.6; % Hoechst
NarrowControlChannelBackgroundBins(3,:) = NarrowAPbins >= 0.8; % Bicoid
NarrowControlChannelBackgroundBins(4,:) = NarrowAPbins >= 0.775; % Knirps
NarrowControlChannelBackgroundBins(5,:) = NarrowAPbins >= 0.55 & NarrowAPbins <= 0.625; % Hunchback
NarrowTestChannelBackgroundBins = NarrowControlChannelBackgroundBins;
NarrowTestChannelBackgroundBins(4,:) =  NarrowAPbins >= 0.1 & NarrowAPbins < 0.9;
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
FitInfo = cell(NSlides, NSlides, NChannels);
ScalingFactors = NaN(NSlides, NSlides, NChannels);
ScalingSEs = NaN(NSlides, NSlides, NChannels);
CompiledEmbryos.NarrowSlideRescalingFactors(1,:) = 1;
CompiledEmbryos.NarrowSlideRescalingSEs(1,:) = 0;
for ch_index = 2:NChannels
    if ~SetIsFlipped
        MeanChannelDorsalProfilesAllSlides = NaN( NSlides, sum(NarrowTestChannelBackgroundBins(ch_index,:)));
        StdChannelDorsalProfilesAllSlides = NaN( NSlides, sum(NarrowTestChannelBackgroundBins(ch_index,:)));
    else
        MeanChannelDorsalProfilesAllSlides = NaN( NSlides, sum(NarrowControlChannelBackgroundBins(ch_index,:)));
        StdChannelDorsalProfilesAllSlides = NaN( NSlides, sum(NarrowControlChannelBackgroundBins(ch_index,:)));
    end
    CountsDorsalProfilesAllSlides  = NaN(1, NSlides);
    for slide_index = 1:NSlides
        if ~SetIsFlipped
            ProfilesToUse = CompiledEmbryos.SlideTestSetEmbryos(slide_index,:) & CompiledEmbryos.IsNC14;
            DorsalProfiles = CompiledEmbryos.DorsalAvgNarrowAPProfiles(ProfilesToUse,NarrowTestChannelBackgroundBins(ch_index,:),ch_index);
        else
            ProfilesToUse = CompiledEmbryos.SlideControlSetEmbryos(slide_index,:) & CompiledEmbryos.IsNC14;
            DorsalProfiles = CompiledEmbryos.DorsalAvgNarrowAPProfiles(ProfilesToUse,NarrowControlChannelBackgroundBins(ch_index,:),ch_index);
        end
        MeanChannelDorsalProfilesAllSlides(slide_index,:) =  mean(DorsalProfiles,1, 'omitnan');
        StdChannelDorsalProfilesAllSlides(slide_index,:) =  std(DorsalProfiles,1, 'omitnan');
        CountsDorsalProfilesAllSlides(slide_index) = size(DorsalProfiles, 1);
        CompiledEmbryos.NarrowMeanTestBackgrounds(slide_index, ch_index) = mean(MeanChannelDorsalProfilesAllSlides(slide_index,:), 'omitnan');
        CompiledEmbryos.NarrowCountTestBackgrounds(slide_index, ch_index) = CountsDorsalProfilesAllSlides(slide_index);
    end
    
    
    for i = 1:(NSlides-1)
        for j = i+1:NSlides
            FitSet1 = MeanChannelDorsalProfilesAllSlides(i,:);
            FitSet1 = FitSet1(:);
            FitSet2 = MeanChannelDorsalProfilesAllSlides(j,:);
            FitSet2 = FitSet2(:);
            dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
            FitInfo{j, i, ch_index} = dlm;
            ScalingFactors(j, i, ch_index) = dlm.Coefficients.Estimate;
            ScalingSEs(j, i, ch_index) = dlm.Coefficients.SE;
            
            
        end
    end
    
    for i = 1:(NSlides-1)
        for j = i+1:NSlides
            
            ScalingFactors(i, j, ch_index) = 1/ScalingFactors(j, i, ch_index);
            ScalingSEs(i, j, ch_index) = ScalingSEs(j, i, ch_index)/ScalingFactors(j, i, ch_index);
            
        end
    end
    
    for j = 2:NSlides
        CompiledEmbryos.NarrowSlideRescalingFactors(j, ch_index) = ScalingFactors(j, 1, ch_index);
        CompiledEmbryos.NarrowSlideRescalingSEs(j, ch_index) = ScalingSEs(j, 1, ch_index);
    end
    
    
    
end


%%
CompiledEmbryos.SlideRescaledDorsalNarrowAPProfiles = CompiledEmbryos.DorsalNarrowAPProfiles;
CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles = CompiledEmbryos.DorsalAvgNarrowAPProfiles;
for embryo_index = AllEmbryos
    for ch_index = 2:NChannels
        CompiledEmbryos.SlideRescaledDorsalNarrowAPProfiles(embryo_index,:,ch_index) =...
            CompiledEmbryos.NarrowSlideRescalingFactors(CompiledEmbryos.SlideIDs(embryo_index), ch_index)*CompiledEmbryos.SlideRescaledDorsalNarrowAPProfiles(embryo_index,:,ch_index);
        CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(embryo_index,:,ch_index) =...
            CompiledEmbryos.NarrowSlideRescalingFactors(CompiledEmbryos.SlideIDs(embryo_index), ch_index)*CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(embryo_index,:,ch_index);
    end
end



