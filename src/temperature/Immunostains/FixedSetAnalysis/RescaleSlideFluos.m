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

f = @(b,x) b(1).*x + b(2);
AllEmbryos = 1:NEmbryos;
SlideIDs = unique(CompiledEmbryos.SlideIDs);
NSlides = length(SlideIDs);
ControlChannelBackgroundBins = zeros(NChannels, length(APbins), 'logical');
ControlChannelBackgroundBins(2,:) = APbins >= 0.1 & APbins <= 0.9;%APbins >= 0.375 & APbins <= 0.6; % Hoechst
ControlChannelBackgroundBins(3,:) = APbins >= 0.1 & APbins <= 0.9;%APbins >= 0.8; % Bicoid
ControlChannelBackgroundBins(4,:) = APbins >= 0.1 & APbins <= 0.9;%APbins >= 0.775; % Knirps
ControlChannelBackgroundBins(5,:) = APbins >= 0.1 & APbins <= 0.9;%APbins >= 0.55 & APbins <= 0.625; % Hunchback
MinBins = zeros(NChannels, length(APbins), 'logical');
MinBins(2,:) = APbins >= 0.2 & APbins <= 0.8;
MinBins(3,:) = APbins >= 0.85 & APbins <= 0.9;
MinBins(4,:) = APbins >= 0.85 & APbins <= 0.9;
MinBins(5,:) = APbins >= 0.625 & APbins <= 0.675;
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

xfits = 0:1:70;
Nxfits = length(xfits);
min_2sigma_points = 3;
counts_above_limit = 1;
counts_below_limit = 1;
sigma = 5;
window_multiplier = 4;
NumBootstrappedFits = 200;
NumPoints = 100;
NumAPbins = size(CompiledEmbryos.DorsalAvgAPProfiles, 2);
NChannels =size(CompiledEmbryos.DorsalAvgAPProfiles, 3);


NValidPoints = NaN(NChannels, NSlides);


%%
for ch_index = 3
    for slide_index = 1:NSlides
        UseTF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes) & (CompiledEmbryos.SlideIDs == slide_index);
        x_sample = CompiledEmbryos.DubuisEmbryoTimes(UseTF);
        ys = CompiledEmbryos.DorsalAvgAPProfiles(UseTF, :,ch_index);
        BadTF =  sum(isnan(ys), 2).' >NumAPbins/2;
        x_sample = x_sample(~BadTF);
        ys = ys(~BadTF,:);
        
        
        
        [x_sample, sort_order] = sort(x_sample);
        ys = ys(sort_order,:);
        
        SmoothedProfiles = NaN(Nxfits, NumAPbins, NumBootstrappedFits);
        
        DiffMat = xfits.'-x_sample;
        GaussianWeights = GetGaussianWeightMat(xfits, x_sample, sigma, window_multiplier*sigma);
        counts = zeros(1, Nxfits);
        counts_above = zeros(1,Nxfits);
        counts_below = zeros(1, Nxfits);
        for x_index = 1:Nxfits
            MatchedPoints = find(GaussianWeights(x_index,:)> 0);
            Num2SigmaPoints = sum(DiffMat(x_index,MatchedPoints) > -2*sigma & DiffMat(x_index,MatchedPoints) < 2*sigma );
            
            if Num2SigmaPoints < min_2sigma_points
                GaussianWeights(x_index,:) = 0;
                continue
            end
            
            counts(x_index) = Num2SigmaPoints;
            MatchedPoints = find(GaussianWeights(x_index,:)> 0);
            
            counts_above(x_index) = sum(DiffMat(x_index,MatchedPoints) > 0 & DiffMat(x_index,MatchedPoints) < 2*sigma );
            counts_below(x_index) = sum(DiffMat(x_index,MatchedPoints) < 0 & DiffMat(x_index,MatchedPoints) > -2*sigma );
            
            if counts_above(x_index) < counts_above_limit
                GaussianWeights(x_index,:) = 0;
            end
            
            if counts_below(x_index) < counts_below_limit
                GaussianWeights(x_index,:) = 0;
            end
            
        end
        
        
        ValidCountTFs = (counts >= min_2sigma_points & counts_above >= counts_above_limit & counts_below >= counts_below_limit & xfits >= min(x_sample)+sigma/2 & xfits <= max(x_sample) - sigma/2);
        x_indices = 1:Nxfits;
        ValidXfits = x_indices(ValidCountTFs);
        NValidXfits = length(ValidXfits);
        NValidPoints(ch_index, slide_index) = NValidXfits;
    end
end
ref_slide = find(NValidPoints(3,:) == max(NValidPoints(3,:)));

%% Bootstrapping the FitSlideRescaledDorsalAvgAPProfiles for the Test Set

xfits = 0:0.1:70;
Nxfits = length(xfits);
if ~isfield(CompiledEmbryos, 'BootstrappedProfiles')
    CompiledEmbryos.BootstrappedProfiles = {};
end
if ~isfield(CompiledEmbryos.BootstrappedProfiles, 'SlideFitProfiles')
    CompiledEmbryos.BootstrappedProfiles.SlideFitProfiles = {};
end
if ~isfield(CompiledEmbryos.BootstrappedProfiles.SlideFitProfiles, 'Test')
    CompiledEmbryos.BootstrappedProfiles.SlideFitProfiles.Test = {};
end

CompiledEmbryos.BootstrappedProfiles.SlideFitProfiles.Test.mean = NaN(Nxfits, NumAPbins, NChannels, NSlides);
CompiledEmbryos.BootstrappedProfiles.SlideFitProfiles.Test.se = NaN(Nxfits, NumAPbins, NChannels, NSlides);
CompiledEmbryos.BootstrappedProfiles.SlideFitProfiles.Test.counts = NaN(Nxfits,NChannels, NSlides);
CompiledEmbryos.BootstrappedProfiles.SlideFitProfiles.Test.counts_above = NaN(Nxfits,NChannels, NSlides);
CompiledEmbryos.BootstrappedProfiles.SlideFitProfiles.Test.counts_below = NaN(Nxfits,NChannels, NSlides);
%%
for ch_index = [2 3 4 5]
    for slide_index = ref_slide
        UseTF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes) & (CompiledEmbryos.SlideIDs == slide_index);
        x_sample = CompiledEmbryos.DubuisEmbryoTimes(UseTF);
        ys = CompiledEmbryos.DorsalAvgAPProfiles(UseTF, :,ch_index);
        BadTF =  sum(isnan(ys), 2).' > NumAPbins/2;
        x_sample = x_sample(~BadTF);
        ys = ys(~BadTF,:);
        
        
        
        [x_sample, sort_order] = sort(x_sample);
        ys = ys(sort_order,:);
        
        SmoothedProfiles = NaN(Nxfits, NumAPbins, NumBootstrappedFits);
        
        DiffMat = xfits.'-x_sample;
        GaussianWeights = GetGaussianWeightMat(xfits, x_sample, sigma, window_multiplier*sigma);
        counts = zeros(1, Nxfits);
        counts_above = zeros(1,Nxfits);
        counts_below = zeros(1, Nxfits);
        for x_index = 1:Nxfits
            MatchedPoints = find(GaussianWeights(x_index,:)> 0);
            Num2SigmaPoints = sum(DiffMat(x_index,MatchedPoints) > -2*sigma & DiffMat(x_index,MatchedPoints) < 2*sigma );
            
            if Num2SigmaPoints < min_2sigma_points
                GaussianWeights(x_index,:) = 0;
                continue
            end
            
            counts(x_index) = Num2SigmaPoints;
            MatchedPoints = find(GaussianWeights(x_index,:)> 0);
            
            counts_above(x_index) = sum(DiffMat(x_index,MatchedPoints) > 0 & DiffMat(x_index,MatchedPoints) < 2*sigma );
            counts_below(x_index) = sum(DiffMat(x_index,MatchedPoints) < 0 & DiffMat(x_index,MatchedPoints) > -2*sigma );
            
            if counts_above(x_index) < counts_above_limit
                GaussianWeights(x_index,:) = 0;
            end
            
            if counts_below(x_index) < counts_below_limit
                GaussianWeights(x_index,:) = 0;
            end
            
        end
        
        
        ValidCountTFs = (counts >= min_2sigma_points & counts_above >= counts_above_limit & counts_below >= counts_below_limit & xfits >= min(x_sample)+sigma/2 & xfits <= max(x_sample) - sigma/2);
        x_indices = 1:Nxfits;
        ValidXfits = x_indices(ValidCountTFs);
        NValidXfits = length(ValidXfits);
        NValidPoints(ch_index, slide_index) = NValidXfits;
        
        for j = 1:NValidXfits
            WindowWidthLimits = [xfits(ValidXfits(j))-( min(x_sample)), max(x_sample)- xfits(ValidXfits(j))];
            for rep = 1:NumBootstrappedFits*2
                r1  = (rand(NumPoints, 1)*(2*window_multiplier*sigma)-window_multiplier*sigma).';
                Deltas = DiffMat(ValidXfits(j),:);
                Deltas(Deltas > window_multiplier*sigma | Deltas < -window_multiplier*sigma) = NaN;
                ind = NaN(1, NumPoints);
                for k = 1:NumPoints
                    ind(k) = find(min(abs(Deltas - r1(k))) == abs(Deltas - r1(k)));
                end
                SmoothedProfiles(ValidXfits(j),:,rep) =  GaussianWeights(ValidXfits(j),ind)*ys(ind, :)/(sum(GaussianWeights(ValidXfits(j),ind)));
            end
        end
        
        
        MeanSmoothedProfiles = mean(SmoothedProfiles, 3, 'omitnan');
        SmoothedProfileSEs = std(SmoothedProfiles,0,  3, 'omitnan')/sqrt(NumBootstrappedFits);
        
        SmoothedProfileSEs(MeanSmoothedProfiles  == 0) = NaN;
        MeanSmoothedProfiles(MeanSmoothedProfiles  == 0) = NaN;
        %close all
        
        CompiledEmbryos.BootstrappedProfiles.SlideFitProfiles.x = xfits;
        CompiledEmbryos.BootstrappedProfiles.SlideFitProfiles.Test.mean(:,:,ch_index, slide_index) = MeanSmoothedProfiles;
        CompiledEmbryos.BootstrappedProfiles.SlideFitProfiles.Test.se(:,:,ch_index, slide_index)  = SmoothedProfileSEs;
        CompiledEmbryos.BootstrappedProfiles.SlideFitProfiles.Test.counts(:,ch_index,slide_index) = counts;
        CompiledEmbryos.BootstrappedProfiles.SlideFitProfiles.Test.counts_above(:,ch_index,slide_index)  = counts_above;
        CompiledEmbryos.BootstrappedProfiles.SlideFitProfiles.Test.counts_below(:,ch_index,slide_index)  = counts_above;
        
    end
end
%%
CompiledEmbryos.SlideRescalingFactors(ref_slide,:) = 1;
CompiledEmbryos.SlideRescalingSEs(ref_slide,:) = 0;
CompiledEmbryos.SlideRescalingIntercepts(ref_slide,:) = 0;
CompiledEmbryos.SlideRescalingInterceptSEs(ref_slide,:) = 0;
AllSlides = 1:NSlides;
SlidesToFit = AllSlides(AllSlides ~= ref_slide);
for slide_index = SlidesToFit
    for ch_index = 2:5
        UseTF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes) & (CompiledEmbryos.SlideIDs == slide_index);
        x_sample = round(CompiledEmbryos.DubuisEmbryoTimes(UseTF), 1);
        ys = CompiledEmbryos.DorsalAvgAPProfiles(UseTF, :,ch_index);
        BadTF =  sum(isnan(ys), 2).' > NumAPbins/2;
        x_sample = x_sample(~BadTF);
        ys = ys(~BadTF,:);
        
        
        
        [x_sample, sort_order] = sort(x_sample);
        ys = ys(sort_order,:);
        
        y_pred = GetMasterProfileForEmbryoTimes(x_sample, CompiledEmbryos.BootstrappedProfiles.SlideFitProfiles.Test.mean(:,:,ch_index, ref_slide) , xfits);
        
        ys = ys(:);
        y_pred = y_pred(:);
        TFys = ~isnan(ys) & ~isnan(y_pred);
        ys = ys(TFys);
        y_pred = y_pred(TFys);
        beta0 =[1, min(ys,[], 'all')];
        if ~isempty(y_pred)
            dlm = fitnlm(y_pred(:),ys(:),f,beta0);
            CompiledEmbryos.SlideRescalingFactors(slide_index, ch_index) = 1/dlm.Coefficients.Estimate(1);
            CompiledEmbryos.SlideRescalingIntercepts(slide_index, ch_index) = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
            CompiledEmbryos.SlideRescalingSEs(slide_index, ch_index)= dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
            CompiledEmbryos.SlideRescalingInterceptSEs(slide_index, ch_index)= sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
                dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
        end
        
    end
end
   
%%
CompiledEmbryos.SlideRescaledDorsalAPProfiles = CompiledEmbryos.DorsalAPProfiles;
CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles = CompiledEmbryos.DorsalAvgAPProfiles;
for ch_index = 2:5
    for slide_index = 1:NSlides
        CompiledEmbryos.SlideRescaledDorsalAPProfiles(CompiledEmbryos.SlideIDs == slide_index,:,ch_index) = ...
            CompiledEmbryos.SlideRescalingFactors(slide_index, ch_index)*CompiledEmbryos.SlideRescaledDorsalAPProfiles(CompiledEmbryos.SlideIDs == slide_index,:,ch_index)+...
            CompiledEmbryos.SlideRescalingIntercepts(slide_index, ch_index);
        CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.SlideIDs == slide_index,:,ch_index) = ...
            CompiledEmbryos.SlideRescalingFactors(slide_index, ch_index)*CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.SlideIDs == slide_index,:,ch_index)+...
            CompiledEmbryos.SlideRescalingIntercepts(slide_index, ch_index);
    end
end


%%
NarrowAPbins = 0:0.0125:1;
NumNarrowAPbins = length(NarrowAPbins);

NarrowControlChannelBackgroundBins = zeros(NChannels, length(NarrowAPbins), 'logical');
NarrowControlChannelBackgroundBins(2,:) = NarrowAPbins >= 0.1 & NarrowAPbins <= 0.9;%APbins >= 0.375 & APbins <= 0.6; % Hoechst
NarrowControlChannelBackgroundBins(3,:) = NarrowAPbins >= 0.1 & NarrowAPbins <= 0.9;%APbins >= 0.8; % Bicoid
NarrowControlChannelBackgroundBins(4,:) = NarrowAPbins >= 0.1 & NarrowAPbins <= 0.9;%APbins >= 0.775; % Knirps
NarrowControlChannelBackgroundBins(5,:) = NarrowAPbins >= 0.1 & NarrowAPbins <= 0.9;%APbins >= 0.55 & APbins <= 0.625; % Hunchback
MinNarrowBins = zeros(NChannels, length(NarrowAPbins), 'logical');
MinNarrowBins(2,:) = NarrowAPbins >= 0.2 & NarrowAPbins <= 0.8;
MinNarrowBins(3,:) = NarrowAPbins >= 0.85 & NarrowAPbins <= 0.9;
MinNarrowBins(4,:) = NarrowAPbins >= 0.85 & NarrowAPbins <= 0.9;
MinNarrowBins(5,:) = NarrowAPbins >= 0.625 & NarrowAPbins <= 0.675;
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
CompiledEmbryos.SlideRescalingFactors(1,:) = 1;
CompiledEmbryos.SlideRescalingSEs(1,:) = 0;
CompiledEmbryos.SlideRescalingIntercepts(1,:) = 0;
CompiledEmbryos.SlideRescalingInterceptSEs(1,:) = 0;

%%

xfits = 0:1:70;
Nxfits = length(xfits);
min_2sigma_points = 3;
counts_above_limit = 1;
counts_below_limit = 1;
sigma = 5;
window_multiplier = 4;
NumBootstrappedFits = 200;
NumPoints = 100;
NumNarrowAPbins = size(CompiledEmbryos.DorsalAvgNarrowAPProfiles, 2);
NChannels =size(CompiledEmbryos.DorsalAvgNarrowAPProfiles, 3);


NValidPoints = NaN(NChannels, NSlides);


%%
for ch_index = 3
    for slide_index = 1:NSlides
        UseTF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes) & (CompiledEmbryos.SlideIDs == slide_index);
        x_sample = CompiledEmbryos.DubuisEmbryoTimes(UseTF);
        ys = CompiledEmbryos.DorsalAvgNarrowAPProfiles(UseTF, :,ch_index);
        BadTF =  sum(isnan(ys), 2).' > NumNarrowAPbins/2;
        x_sample = x_sample(~BadTF);
        ys = ys(~BadTF,:);
        
        
        
        [x_sample, sort_order] = sort(x_sample);
        ys = ys(sort_order,:);
        
        SmoothedProfiles = NaN(Nxfits, NumNarrowAPbins, NumBootstrappedFits);
        
        DiffMat = xfits.'-x_sample;
        GaussianWeights = GetGaussianWeightMat(xfits, x_sample, sigma, window_multiplier*sigma);
        counts = zeros(1, Nxfits);
        counts_above = zeros(1,Nxfits);
        counts_below = zeros(1, Nxfits);
        for x_index = 1:Nxfits
            MatchedPoints = find(GaussianWeights(x_index,:)> 0);
            Num2SigmaPoints = sum(DiffMat(x_index,MatchedPoints) > -2*sigma & DiffMat(x_index,MatchedPoints) < 2*sigma );
            
            if Num2SigmaPoints < min_2sigma_points
                GaussianWeights(x_index,:) = 0;
                continue
            end
            
            counts(x_index) = Num2SigmaPoints;
            MatchedPoints = find(GaussianWeights(x_index,:)> 0);
            
            counts_above(x_index) = sum(DiffMat(x_index,MatchedPoints) > 0 & DiffMat(x_index,MatchedPoints) < 2*sigma );
            counts_below(x_index) = sum(DiffMat(x_index,MatchedPoints) < 0 & DiffMat(x_index,MatchedPoints) > -2*sigma );
            
            if counts_above(x_index) < counts_above_limit
                GaussianWeights(x_index,:) = 0;
            end
            
            if counts_below(x_index) < counts_below_limit
                GaussianWeights(x_index,:) = 0;
            end
            
        end
        
        
        ValidCountTFs = (counts >= min_2sigma_points & counts_above >= counts_above_limit & counts_below >= counts_below_limit & xfits >= min(x_sample)+sigma/2 & xfits <= max(x_sample) - sigma/2);
        x_indices = 1:Nxfits;
        ValidXfits = x_indices(ValidCountTFs);
        NValidXfits = length(ValidXfits);
        NValidPoints(ch_index, slide_index) = NValidXfits;
    end
end
ref_slide = find(NValidPoints(3,:) == max(NValidPoints(3,:)));

%% Bootstrapping the FitSlideRescaledDorsalAvgAPProfiles for the Test Set

xfits = 0:0.1:70;
Nxfits = length(xfits);
if ~isfield(CompiledEmbryos, 'BootstrappedProfiles')
    CompiledEmbryos.BootstrappedProfiles = {};
end
if ~isfield(CompiledEmbryos.BootstrappedProfiles, 'SlideFitProfiles')
    CompiledEmbryos.BootstrappedProfiles.NarrowSlideFitProfiles = {};
end
if ~isfield(CompiledEmbryos.BootstrappedProfiles.SlideFitProfiles, 'Test')
    CompiledEmbryos.BootstrappedProfiles.NarrowSlideFitProfiles.Test = {};
end

CompiledEmbryos.BootstrappedProfiles.NarrowSlideFitProfiles.Test.mean = NaN(Nxfits, NumNarrowAPbins, NChannels, NSlides);
CompiledEmbryos.BootstrappedProfiles.NarrowSlideFitProfiles.Test.se = NaN(Nxfits, NumNarrowAPbins, NChannels, NSlides);
CompiledEmbryos.BootstrappedProfiles.NarrowSlideFitProfiles.Test.counts = NaN(Nxfits,NChannels, NSlides);
CompiledEmbryos.BootstrappedProfiles.NarrowSlideFitProfiles.Test.counts_above = NaN(Nxfits,NChannels, NSlides);
CompiledEmbryos.BootstrappedProfiles.NarrowSlideFitProfiles.Test.counts_below = NaN(Nxfits,NChannels, NSlides);
%%
for ch_index = [2 3 4 5]
    for slide_index = ref_slide
        UseTF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes) & (CompiledEmbryos.SlideIDs == slide_index);
        x_sample = CompiledEmbryos.DubuisEmbryoTimes(UseTF);
        ys = CompiledEmbryos.DorsalAvgNarrowAPProfiles(UseTF, :,ch_index);
        BadTF =  sum(isnan(ys), 2).' > NumNarrowAPbins/2;
        x_sample = x_sample(~BadTF);
        ys = ys(~BadTF,:);
        
        
        
        [x_sample, sort_order] = sort(x_sample);
        ys = ys(sort_order,:);
        
        SmoothedProfiles = NaN(Nxfits, NumNarrowAPbins, NumBootstrappedFits);
        
        DiffMat = xfits.'-x_sample;
        GaussianWeights = GetGaussianWeightMat(xfits, x_sample, sigma, window_multiplier*sigma);
        counts = zeros(1, Nxfits);
        counts_above = zeros(1,Nxfits);
        counts_below = zeros(1, Nxfits);
        for x_index = 1:Nxfits
            MatchedPoints = find(GaussianWeights(x_index,:)> 0);
            Num2SigmaPoints = sum(DiffMat(x_index,MatchedPoints) > -2*sigma & DiffMat(x_index,MatchedPoints) < 2*sigma );
            
            if Num2SigmaPoints < min_2sigma_points
                GaussianWeights(x_index,:) = 0;
                continue
            end
            
            counts(x_index) = Num2SigmaPoints;
            MatchedPoints = find(GaussianWeights(x_index,:)> 0);
            
            counts_above(x_index) = sum(DiffMat(x_index,MatchedPoints) > 0 & DiffMat(x_index,MatchedPoints) < 2*sigma );
            counts_below(x_index) = sum(DiffMat(x_index,MatchedPoints) < 0 & DiffMat(x_index,MatchedPoints) > -2*sigma );
            
            if counts_above(x_index) < counts_above_limit
                GaussianWeights(x_index,:) = 0;
            end
            
            if counts_below(x_index) < counts_below_limit
                GaussianWeights(x_index,:) = 0;
            end
            
        end
        
        
        ValidCountTFs = (counts >= min_2sigma_points & counts_above >= counts_above_limit & counts_below >= counts_below_limit & xfits >= min(x_sample)+sigma/2 & xfits <= max(x_sample) - sigma/2);
        x_indices = 1:Nxfits;
        ValidXfits = x_indices(ValidCountTFs);
        NValidXfits = length(ValidXfits);
        NValidPoints(ch_index, slide_index) = NValidXfits;
        
        for j = 1:NValidXfits
            WindowWidthLimits = [xfits(ValidXfits(j))-( min(x_sample)), max(x_sample)- xfits(ValidXfits(j))];
            for rep = 1:NumBootstrappedFits*2
                r1  = (rand(NumPoints, 1)*(2*window_multiplier*sigma)-window_multiplier*sigma).';
                Deltas = DiffMat(ValidXfits(j),:);
                Deltas(Deltas > window_multiplier*sigma | Deltas < -window_multiplier*sigma) = NaN;
                ind = NaN(1, NumPoints);
                for k = 1:NumPoints
                    ind(k) = find(min(abs(Deltas - r1(k))) == abs(Deltas - r1(k)));
                end
                SmoothedProfiles(ValidXfits(j),:,rep) =  GaussianWeights(ValidXfits(j),ind)*ys(ind, :)/(sum(GaussianWeights(ValidXfits(j),ind)));
            end
        end
        
        
        MeanSmoothedProfiles = mean(SmoothedProfiles, 3, 'omitnan');
        SmoothedProfileSEs = std(SmoothedProfiles,0,  3, 'omitnan')/sqrt(NumBootstrappedFits);
        
        SmoothedProfileSEs(MeanSmoothedProfiles  == 0) = NaN;
        MeanSmoothedProfiles(MeanSmoothedProfiles  == 0) = NaN;
        %close all
        
        CompiledEmbryos.BootstrappedProfiles.NarrowSlideFitProfiles.x = xfits;
        CompiledEmbryos.BootstrappedProfiles.NarrowSlideFitProfiles.Test.mean(:,:,ch_index, slide_index) = MeanSmoothedProfiles;
        CompiledEmbryos.BootstrappedProfiles.NarrowSlideFitProfiles.Test.se(:,:,ch_index, slide_index)  = SmoothedProfileSEs;
        CompiledEmbryos.BootstrappedProfiles.NarrowSlideFitProfiles.Test.counts(:,ch_index,slide_index) = counts;
        CompiledEmbryos.BootstrappedProfiles.NarrowSlideFitProfiles.Test.counts_above(:,ch_index,slide_index)  = counts_above;
        CompiledEmbryos.BootstrappedProfiles.NarrowSlideFitProfiles.Test.counts_below(:,ch_index,slide_index)  = counts_above;
        
    end
end
%%
CompiledEmbryos.NarrowSlideRescalingFactors(ref_slide,:) = 1;
CompiledEmbryos.NarrowSlideRescalingSEs(ref_slide,:) = 0;
CompiledEmbryos.NarrowSlideRescalingIntercepts(ref_slide,:) = 0;
CompiledEmbryos.NarrowSlideRescalingInterceptSEs(ref_slide,:) = 0;
AllSlides = 1:NSlides;
SlidesToFit = AllSlides(AllSlides ~= ref_slide);
for slide_index = SlidesToFit
    for ch_index = 2:5
        UseTF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes) & (CompiledEmbryos.SlideIDs == slide_index);
        x_sample = round(CompiledEmbryos.DubuisEmbryoTimes(UseTF), 1);
        ys = CompiledEmbryos.DorsalAvgNarrowAPProfiles(UseTF, :,ch_index);
        BadTF =  sum(isnan(ys), 2).' >NumNarrowAPbins/2;
        x_sample = x_sample(~BadTF);
        ys = ys(~BadTF,:);
        
        
        
        [x_sample, sort_order] = sort(x_sample);
        ys = ys(sort_order,:);
        
        y_pred = GetMasterProfileForEmbryoTimes(x_sample, CompiledEmbryos.BootstrappedProfiles.NarrowSlideFitProfiles.Test.mean(:,:,ch_index, ref_slide) , xfits);
        
        ys = ys(:);
        y_pred = y_pred(:);
        TFys = ~isnan(ys) & ~isnan(y_pred);
        ys = ys(TFys);
        y_pred = y_pred(TFys);
        beta0 =[1, min(ys,[], 'all')];
        if ~isempty(y_pred)
            dlm = fitnlm(y_pred(:),ys(:),f,beta0);
            CompiledEmbryos.NarrowSlideRescalingFactors(slide_index, ch_index) = 1/dlm.Coefficients.Estimate(1);
            CompiledEmbryos.NarrowSlideRescalingIntercepts(slide_index, ch_index) = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
            CompiledEmbryos.NarrowSlideRescalingSEs(slide_index, ch_index)= dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
            CompiledEmbryos.NarrowSlideRescalingInterceptSEs(slide_index, ch_index)= sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
                dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
        end
        
    end
end
   
%%
CompiledEmbryos.SlideRescaledDorsalNarrowAPProfiles = CompiledEmbryos.DorsalNarrowAPProfiles;
CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles = CompiledEmbryos.DorsalAvgNarrowAPProfiles;
for ch_index = 2:5
    for slide_index = 1:NSlides
        CompiledEmbryos.SlideRescaledDorsalNarrowAPProfiles(CompiledEmbryos.SlideIDs == slide_index,:,ch_index) = ...
            CompiledEmbryos.NarrowSlideRescalingFactors(slide_index, ch_index)*CompiledEmbryos.SlideRescaledDorsalNarrowAPProfiles(CompiledEmbryos.SlideIDs == slide_index,:,ch_index)+...
            CompiledEmbryos.NarrowSlideRescalingIntercepts(slide_index, ch_index);
        CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(CompiledEmbryos.SlideIDs == slide_index,:,ch_index) = ...
            CompiledEmbryos.NarrowSlideRescalingFactors(slide_index, ch_index)*CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles(CompiledEmbryos.SlideIDs == slide_index,:,ch_index)+...
            CompiledEmbryos.NarrowSlideRescalingIntercepts(slide_index, ch_index);
    end
end
