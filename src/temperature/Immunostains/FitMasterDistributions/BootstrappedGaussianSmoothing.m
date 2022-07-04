% version = 6;
% AllCompiledPath = ['S:/Gabriella/Dropbox/ProteinProfiles/V',num2str(version),'Profiles/'];
% load(AllCompiledPath, 'AllCompiledEmbryos');

%%
MeanSmoothedProfiles = {};
SmoothedProfileSEs = {};
NormedSmoothProfiles = {};
NormedSmoothProfilesSE = {};

exp_indices = [1 2 3 3 3 6 9 12 14 15];
xfits = 0:.1:70;
APbins = 0:0.025:1;
NumAPbins = length(APbins);
NChannels = 5;


for i = 1:length(exp_indices)
    disp(['i = ', num2str(i)])
exp_index = exp_indices(i);
MeanSmoothedProfiles{i} = NaN(length(xfits), NumAPbins, NChannels);
SmoothedProfileSEs{i} = NaN(length(xfits), NumAPbins, NChannels);
NormedSmoothProfiles{i} = NaN(length(xfits), NumAPbins, NChannels);
NormedSmoothProfilesSE{i} = NaN(length(xfits), NumAPbins, NChannels);
for ch_index = [3 5]
CompiledEmbryos = AllCompiledEmbryos{exp_index};

if i <= 3
    UseTF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);
elseif i == 4
    UseTF = (CompiledEmbryos.TestSetEmbryos | CompiledEmbryos.ControlSetEmbryos) &...
        CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);
else
    UseTF = CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);
end
x_sample = CompiledEmbryos.DubuisEmbryoTimes(UseTF);
if ch_index == 3
    ys = CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles(UseTF, :, ch_index);
else
    ys = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(UseTF, :, ch_index);
end
%%
BadTF =  sum(isnan(ys), 2).' >min( sum(isnan(ys), 2).');
x_sample = x_sample(~BadTF);
ys = ys(~BadTF,:);
min_2sigma_points = 10;
counts_above_limit = 3;
counts_below_limit = 3;
sigma = 5;
window_multiplier = 3;
NumBootstrappedFits = 200;
NumPoints = 100;


[x_sample, sort_order] = sort(x_sample);
ys = ys(sort_order,:);

Nxfits = length(xfits);
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



MeanProf = mean(SmoothedProfiles, 3, 'omitnan');
SEProf = std(SmoothedProfiles,0,  3, 'omitnan')/sqrt(NumBootstrappedFits);
SEProf(MeanProf == 0) = NaN;
MeanProf(MeanProf == 0) = NaN;

SmoothedProfileSEs{i}(:,:,ch_index) = SEProf;
MeanSmoothedProfiles{i}(:,:,ch_index) = MeanProf;
%close all
% figure(exp_index)
% scatter(x_sample, ys(:,13))
% hold on
% % for rep = 1:NumBootstrappedFits
% %     plot(xfits, SmoothedProfiles(:,13,rep) )
% % end
% plot(xfits, MeanSmoothedProfiles{exp_index}(:,13, ch_index), 'k', 'LineWidth', 2)
% hold off

maxProf = max(MeanSmoothedProfiles{i}(:,:,ch_index) , [], 'all', 'omitnan');
minProf = min(MeanSmoothedProfiles{i}(:,:,ch_index) , [], 'all', 'omitnan');

NormedMeanSmoothedProfiles{i}(:,:,ch_index) = (MeanSmoothedProfiles{i}(:,:,ch_index) -minProf)/(maxProf-minProf);
NormedSmoothedProfileSEs{i}(:,:,ch_index) = SmoothedProfileSEs{i}(:,:,ch_index)/(maxProf-minProf);


%%

for j = 21:37
    close all
    figure(i)
scatter(x_sample, ys(:,j))
hold on
% for rep = 1:NumBootstrappedFits
%     plot(xfits, SmoothedProfiles(:,13,rep) )
% end
plot(xfits, MeanSmoothedProfiles{i}(:,j,ch_index), 'k', 'LineWidth', 2)
xlabel('Dubuis Time (m)')
title([num2str(i), ', ', ChannelNames{ch_index}, ', AP: ', num2str(APbins(j))])
hold off
end
end
end
%%

outpath = 'S:/Gabriella/Dropbox/BootstrappedTestData/25CBootstrappedGaussianSmoothedProfiles.mat';
save(outpath, 'MeanSmoothedProfiles', 'SmoothedProfileSEs', 'xfits', 'NormedSmoothedProfileSEs', 'NormedMeanSmoothedProfiles');


