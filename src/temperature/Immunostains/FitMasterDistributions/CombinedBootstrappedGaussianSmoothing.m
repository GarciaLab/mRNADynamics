function [CombinedSlideMeanProfiles, CombinedSlideSEs] =  CombinedBootstrappedGaussianSmoothing(AllCompiledEmbryos)
xfits = 0:.1:70;
APbins = 0:0.025:1;
NumAPbins = length(APbins);
NChannels = 5;
ChannelNames = {'', 'Hoechst', 'Bicoid', 'Knirps', 'Hunchback'};
CombinedSlideMeanProfiles = {};
CombinedSlideSEs = {};


%% Make overall master 25C profile with knirps stained embryos
exp_indices = [1 2 3 6 9 12 15];
SetTypeStrings = {'Test', 'Test', 'Control','Control','Control','Control','Control'};

for i = 1:length(exp_indices)
    exp_index = exp_indices(i);
    if strcmpi(lower(SetTypeStrings{i}), 'test')
        TFset = AllCompiledEmbryos{exp_index}.IsNC14 & AllCompiledEmbryos{exp_index}.TestSetEmbryos  & ~isnan(AllCompiledEmbryos{exp_index}.DubuisEmbryoTimes) & AllCompiledEmbryos{exp_index}.Approved;
    else
        TFset = AllCompiledEmbryos{exp_index}.IsNC14 & AllCompiledEmbryos{exp_index}.ControlSetEmbryos & ~isnan(AllCompiledEmbryos{exp_index}.DubuisEmbryoTimes) & AllCompiledEmbryos{exp_index}.Approved;
    end
    if i == 1
        
        
        AllProfiles = AllCompiledEmbryos{exp_index}.FitSlideRescaledDorsalAvgAPProfiles(TFset,:,:);
        AllProfiles2 = AllCompiledEmbryos{exp_index}.SlideRescaledDorsalAvgAPProfiles(TFset,:,:);
        AllEmbryoTimes = AllCompiledEmbryos{exp_index}.DubuisEmbryoTimes(TFset);
        AllSetIDs = ones(1, length(AllEmbryoTimes));
    else
        SetProfile =  AllCompiledEmbryos{exp_index}.FitSlideRescaledDorsalAvgAPProfiles(TFset,:,:);
        SetProfile2 =  AllCompiledEmbryos{exp_index}.SlideRescaledDorsalAvgAPProfiles(TFset,:,:);
        SetEmbryoTimes = AllCompiledEmbryos{exp_index}.DubuisEmbryoTimes(TFset);
        SubSetIDs = i*ones(1, length(SetEmbryoTimes));
        AllProfiles = [AllProfiles ; SetProfile];
        AllProfiles2 = [AllProfiles2 ; SetProfile2];
        AllEmbryoTimes = [AllEmbryoTimes SetEmbryoTimes];
        AllSetIDs = [AllSetIDs SubSetIDs];
    end
end




f = @(b,x) b(1).*x + b(2);
NSets = length(exp_indices);
SlideRescalingFactors = NaN(NSets, NChannels);
SlideRescalingSEs = NaN(NSets, NChannels);
SlideRescalingIntercepts = NaN(NSets, NChannels);
SlideRescalingInterceptSEs = NaN(NSets, NChannels);
FitInfo = cell(NSets, NChannels);


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


NValidPoints = NaN(NChannels, NSets);


%%
for ch_index = 3
    for slide_index = 1:NSets
        UseTF = AllSetIDs == slide_index;
        x_sample = AllEmbryoTimes(UseTF);
        ys = AllProfiles(UseTF, :,ch_index);
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
SlideRescalingFactors(ref_slide,:) = 1;
SlideRescalingSEs(ref_slide,:) = 0;
SlideRescalingIntercepts(ref_slide,:) = 0;
SlideRescalingInterceptSEs(ref_slide,:) = 0;
xfits = 0:0.1:70;
Nxfits = length(xfits);
SlideFitProfiles = {};
SlideFitProfiles.mean = NaN(Nxfits, NumAPbins, NChannels, NSets);
SlideFitProfiles.se = NaN(Nxfits, NumAPbins, NChannels, NSets);
SlideFitProfiles.counts = NaN(Nxfits,NChannels, NSets);
SlideFitProfiles.counts_above = NaN(Nxfits,NChannels, NSets);
SlideFitProfiles.counts_below = NaN(Nxfits,NChannels, NSets);
%%
for ch_index = [3 5]
    for slide_index = ref_slide
        UseTF = AllSetIDs == slide_index;
        x_sample = AllEmbryoTimes(UseTF);
        if ch_index == 3
            ys = AllProfiles(UseTF, :,ch_index);
        else
            ys = AllProfiles2(UseTF, :,ch_index);
        end
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
        
        SlideFitProfiles.x = xfits;
        SlideFitProfiles.mean(:,:,ch_index, slide_index) = MeanSmoothedProfiles;
        SlideFitProfiles.se(:,:,ch_index, slide_index)  = SmoothedProfileSEs;
        SlideFitProfiles.counts(:,ch_index,slide_index) = counts;
        SlideFitProfiles.counts_above(:,ch_index,slide_index)  = counts_above;
        SlideFitProfiles.counts_below(:,ch_index,slide_index)  = counts_above;
        
    end
end
%%

AllSlides = 1:NSets;
SlidesToFit = AllSlides(AllSlides ~= ref_slide);
for idx = 1:length(SlidesToFit)
    disp(['i = ', num2str(idx)])
    slide_index = SlidesToFit(idx);
    for ch_index = [3 5]
        
        UseTF = AllSetIDs == slide_index;
        x_sample = AllEmbryoTimes(UseTF);
        if ch_index == 3
            ys = AllProfiles(UseTF, :,ch_index);
        else
            ys = AllProfiles2(UseTF, :,ch_index);
        end
        BadTF =  sum(isnan(ys), 2).' >NumAPbins/2;
        x_sample = x_sample(~BadTF);
        ys = ys(~BadTF,:);
        
        
        
        
        [x_sample, sort_order] = sort(x_sample);
        ys = ys(sort_order,:);
        
        if idx == 1
            y_pred = GetMasterProfileForEmbryoTimes(x_sample, SlideFitProfiles.mean(:,:,ch_index, ref_slide) , xfits);
        else
            y_pred = GetMasterProfileForEmbryoTimes(x_sample, SlideFitProfiles.mean(:,:,ch_index, SlidesToFit(idx-1)) , xfits);
        end    
        
        ys = ys(:);
        y_pred = y_pred(:);
        TFys = ~isnan(ys) & ~isnan(y_pred);
        ys = ys(TFys);
        y_pred = y_pred(TFys);
        beta0 =[1, min(ys,[], 'all')];
        if ~isempty(y_pred)
            dlm = fitnlm(y_pred(:),ys(:),f,beta0);
            SlideRescalingFactors(slide_index, ch_index) = 1/dlm.Coefficients.Estimate(1);
            SlideRescalingIntercepts(slide_index, ch_index) = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
            SlideRescalingSEs(slide_index, ch_index)= dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
            SlideRescalingInterceptSEs(slide_index, ch_index)= sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
                dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
        end
        
        AllProfiles(UseTF, :,ch_index) =  AllProfiles(UseTF, :,ch_index)*SlideRescalingFactors(slide_index, ch_index) + SlideRescalingInterceptSEs(slide_index, ch_index);
        AllProfiles2(UseTF, :,ch_index) =  AllProfiles2(UseTF, :,ch_index)*SlideRescalingFactors(slide_index, ch_index) + SlideRescalingInterceptSEs(slide_index, ch_index);
        
        
        UseTF2 = (AllSetIDs <= slide_index) | (AllSetIDs == ref_slide);
        x_sample = AllEmbryoTimes(UseTF2);
        if ch_index == 3
            ys = AllProfiles(UseTF2, :,ch_index);
        else
            ys = AllProfiles2(UseTF2, :,ch_index);
        end
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
        
        SlideFitProfiles.x = xfits;
        SlideFitProfiles.mean(:,:,ch_index, slide_index) = MeanSmoothedProfiles;
        SlideFitProfiles.se(:,:,ch_index, slide_index)  = SmoothedProfileSEs;
        SlideFitProfiles.counts(:,ch_index,slide_index) = counts;
        SlideFitProfiles.counts_above(:,ch_index,slide_index)  = counts_above;
        SlideFitProfiles.counts_below(:,ch_index,slide_index)  = counts_above;
        
    end
end

CombinedSlideMeanProfiles{end+1} = squeeze(SlideFitProfiles.mean(:,:,:, SlidesToFit(end)));
CombinedSlideSEs{end+1} = squeeze(SlideFitProfiles.se(:,:,:, SlidesToFit(end)));


%% Make overall master 25C profile with knirps unstained embryos
exp_indices = [3 14];
SetTypeStrings = {'Test', 'Control'};

for i = 1:length(exp_indices)
    exp_index = exp_indices(i);
    if strcmpi(lower(SetTypeStrings{i}), 'test')
        TFset = AllCompiledEmbryos{exp_index}.IsNC14 & AllCompiledEmbryos{exp_index}.TestSetEmbryos  & ~isnan(AllCompiledEmbryos{exp_index}.DubuisEmbryoTimes) & AllCompiledEmbryos{exp_index}.Approved;
    else
        TFset = AllCompiledEmbryos{exp_index}.IsNC14 & AllCompiledEmbryos{exp_index}.ControlSetEmbryos & ~isnan(AllCompiledEmbryos{exp_index}.DubuisEmbryoTimes) & AllCompiledEmbryos{exp_index}.Approved;
    end
    if i == 1
        
        
        AllProfiles = AllCompiledEmbryos{exp_index}.FitSlideRescaledDorsalAvgAPProfiles(TFset,:,:);
        AllProfiles2 = AllCompiledEmbryos{exp_index}.SlideRescaledDorsalAvgAPProfiles(TFset,:,:);
        AllEmbryoTimes = AllCompiledEmbryos{exp_index}.DubuisEmbryoTimes(TFset);
        AllSetIDs = ones(1, length(AllEmbryoTimes));
    else
        SetProfile =  AllCompiledEmbryos{exp_index}.FitSlideRescaledDorsalAvgAPProfiles(TFset,:,:);
        SetProfile2 =  AllCompiledEmbryos{exp_index}.SlideRescaledDorsalAvgAPProfiles(TFset,:,:);
        SetEmbryoTimes = AllCompiledEmbryos{exp_index}.DubuisEmbryoTimes(TFset);
        SubSetIDs = i*ones(1, length(SetEmbryoTimes));
        AllProfiles = [AllProfiles ; SetProfile];
        AllProfiles2 = [AllProfiles2 ; SetProfile2];
        AllEmbryoTimes = [AllEmbryoTimes SetEmbryoTimes];
        AllSetIDs = [AllSetIDs SubSetIDs];
    end
end




f = @(b,x) b(1).*x + b(2);
NSets = length(exp_indices);
SlideRescalingFactors = NaN(NSets, NChannels);
SlideRescalingSEs = NaN(NSets, NChannels);
SlideRescalingIntercepts = NaN(NSets, NChannels);
SlideRescalingInterceptSEs = NaN(NSets, NChannels);
FitInfo = cell(NSets, NChannels);


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


NValidPoints = NaN(NChannels, NSets);


%%
for ch_index = 3
    for slide_index = 1:NSets
        UseTF = AllSetIDs == slide_index;
        x_sample = AllEmbryoTimes(UseTF);
        ys = AllProfiles(UseTF, :,ch_index);
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
SlideRescalingFactors(ref_slide,:) = 1;
SlideRescalingSEs(ref_slide,:) = 0;
SlideRescalingIntercepts(ref_slide,:) = 0;
SlideRescalingInterceptSEs(ref_slide,:) = 0;
xfits = 0:0.1:70;
Nxfits = length(xfits);
SlideFitProfiles = {};
SlideFitProfiles.mean = NaN(Nxfits, NumAPbins, NChannels, NSets);
SlideFitProfiles.se = NaN(Nxfits, NumAPbins, NChannels, NSets);
SlideFitProfiles.counts = NaN(Nxfits,NChannels, NSets);
SlideFitProfiles.counts_above = NaN(Nxfits,NChannels, NSets);
SlideFitProfiles.counts_below = NaN(Nxfits,NChannels, NSets);
%%
for ch_index = [3 5]
    for slide_index = ref_slide
        UseTF = AllSetIDs == slide_index;
        x_sample = AllEmbryoTimes(UseTF);
        if ch_index == 3
            ys = AllProfiles(UseTF, :,ch_index);
        else
            ys = AllProfiles2(UseTF, :,ch_index);
        end
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
        
        SlideFitProfiles.x = xfits;
        SlideFitProfiles.mean(:,:,ch_index, slide_index) = MeanSmoothedProfiles;
        SlideFitProfiles.se(:,:,ch_index, slide_index)  = SmoothedProfileSEs;
        SlideFitProfiles.counts(:,ch_index,slide_index) = counts;
        SlideFitProfiles.counts_above(:,ch_index,slide_index)  = counts_above;
        SlideFitProfiles.counts_below(:,ch_index,slide_index)  = counts_above;
        
    end
end
%%

AllSlides = 1:NSets;
SlidesToFit = AllSlides(AllSlides ~= ref_slide);
for idx = 1:length(SlidesToFit)
    disp(['i = ', num2str(idx)])
    slide_index = SlidesToFit(idx);
    for ch_index = [3 5]
        
        UseTF = AllSetIDs == slide_index;
        x_sample = AllEmbryoTimes(UseTF);
        if ch_index == 3
            ys = AllProfiles(UseTF, :,ch_index);
        else
            ys = AllProfiles2(UseTF, :,ch_index);
        end
        BadTF =  sum(isnan(ys), 2).' >NumAPbins/2;
        x_sample = x_sample(~BadTF);
        ys = ys(~BadTF,:);
        
        
        
        
        [x_sample, sort_order] = sort(x_sample);
        ys = ys(sort_order,:);
        
        if idx == 1
            y_pred = GetMasterProfileForEmbryoTimes(x_sample, SlideFitProfiles.mean(:,:,ch_index, ref_slide) , xfits);
        else
            y_pred = GetMasterProfileForEmbryoTimes(x_sample, SlideFitProfiles.mean(:,:,ch_index, SlidesToFit(idx-1)) , xfits);
        end    
        
        ys = ys(:);
        y_pred = y_pred(:);
        TFys = ~isnan(ys) & ~isnan(y_pred);
        ys = ys(TFys);
        y_pred = y_pred(TFys);
        beta0 =[1, min(ys,[], 'all')];
        if ~isempty(y_pred)
            dlm = fitnlm(y_pred(:),ys(:),f,beta0);
            SlideRescalingFactors(slide_index, ch_index) = 1/dlm.Coefficients.Estimate(1);
            SlideRescalingIntercepts(slide_index, ch_index) = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
            SlideRescalingSEs(slide_index, ch_index)= dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
            SlideRescalingInterceptSEs(slide_index, ch_index)= sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
                dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
        end
        
        AllProfiles(UseTF, :,ch_index) =  AllProfiles(UseTF, :,ch_index)*SlideRescalingFactors(slide_index, ch_index) + SlideRescalingInterceptSEs(slide_index, ch_index);
        AllProfiles2(UseTF, :,ch_index) =  AllProfiles2(UseTF, :,ch_index)*SlideRescalingFactors(slide_index, ch_index) + SlideRescalingInterceptSEs(slide_index, ch_index);
        
        
        UseTF2 = (AllSetIDs <= slide_index) | (AllSetIDs == ref_slide);
        x_sample = AllEmbryoTimes(UseTF2);
        if ch_index == 3
            ys = AllProfiles(UseTF2, :,ch_index);
        else
            ys = AllProfiles2(UseTF2, :,ch_index);
        end
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
        
        SlideFitProfiles.x = xfits;
        SlideFitProfiles.mean(:,:,ch_index, slide_index) = MeanSmoothedProfiles;
        SlideFitProfiles.se(:,:,ch_index, slide_index)  = SmoothedProfileSEs;
        SlideFitProfiles.counts(:,ch_index,slide_index) = counts;
        SlideFitProfiles.counts_above(:,ch_index,slide_index)  = counts_above;
        SlideFitProfiles.counts_below(:,ch_index,slide_index)  = counts_above;
        
    end
end

CombinedSlideMeanProfiles{end+1} = squeeze(SlideFitProfiles.mean(:,:,:, SlidesToFit(end)));
CombinedSlideSEs{end+1} = squeeze(SlideFitProfiles.se(:,:,:, SlidesToFit(end)));
