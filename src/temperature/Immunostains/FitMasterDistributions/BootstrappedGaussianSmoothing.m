version = 6;
AllCompiledPath = ['S:/Gabriella/Dropbox/ProteinProfiles/V',num2str(version),'Profiles/'];
save(AllCompiledPath, 'AllCompiledEmbryos');
%%
MeanSmoothedProfiles = {};
SmoothedProfileSEs = {};

for exp_index = 1:3
APbin = 13;
UseTF = TestSetTFs{exp_index} & NC14TFs{exp_index} & ~isnan(x_samples{exp_index});
x_sample = x_samples{exp_index}(UseTF);
ys = FullFitBcdProfs2{exp_index}(UseTF, :);
NumAPbins = length(APbins);
BadTF =  sum(isnan(ys), 2).' >min( sum(isnan(ys), 2).');
x_sample = x_sample(~BadTF);
ys = ys(~BadTF,:);
min_2sigma_points = 10;
sigma = 5;
window_multiplier = 3;
NumBootstrappedFits = 1000;
NumPoints = 100;


[x_sample, sort_order] = sort(x_sample);
ys = ys(sort_order,:);
xfits = 0:.1:70;
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
    
    if counts_above(x_index) < 5
        GaussianWeights(x_index,:) = 0;
    end
    
    if counts_below(x_index) < 5
        GaussianWeights(x_index,:) = 0;
    end

end


ValidCountTFs = (counts >= min_2sigma_points & counts_above >= 5 & counts_below >= 5 & xfits >= min(x_sample)+2*sigma & xfits <= max(x_sample) - 2*sigma);
x_indices = 1:Nxfits;
ValidXfits = x_indices(ValidCountTFs);
NValidXfits = length(ValidXfits);

for j = 1:NValidXfits
    WindowWidthLimits = [xfits(ValidXfits(j))-( min(x_sample)), max(x_sample)- xfits(ValidXfits(j))];
for rep = 1:NumBootstrappedFits*2
r1  = (rand(NumPoints, 1)*(2*min([window_multiplier*sigma WindowWidthLimits ]))-min([window_multiplier*sigma WindowWidthLimits ])).';
Deltas = DiffMat(ValidXfits(j),:);
Deltas(Deltas > window_multiplier*sigma | Deltas < -window_multiplier*sigma) = NaN;
ind = NaN(1, NumPoints);
for k = 1:NumPoints
ind(k) = find(min(abs(Deltas - r1(k))) == abs(Deltas - r1(k)));
end


SmoothedProfiles(ValidXfits(j),:,rep) =  GaussianWeights(ValidXfits(j),ind)*ys(ind, :)/(sum(GaussianWeights(ValidXfits(j),ind)));
end
end



MeanSmoothedProfiles{exp_index} = mean(SmoothedProfiles, 3, 'omitnan');
SmoothedProfileSEs{exp_index} = std(SmoothedProfiles,0,  3, 'omitnan')/sqrt(NumBootstrappedFits);

SmoothedProfileSEs{exp_index}(MeanSmoothedProfiles{exp_index}  == 0) = NaN;
MeanSmoothedProfiles{exp_index}(MeanSmoothedProfiles{exp_index}  == 0) = NaN;
%close all
figure(exp_index)
scatter(x_sample, ys(:,13))
hold on
% for rep = 1:NumBootstrappedFits
%     plot(xfits, SmoothedProfiles(:,13,rep) )
% end
plot(xfits, MeanSmoothedProfiles{exp_index}(:,13), 'k', 'LineWidth', 2)
hold off

maxProf = max(MeanSmoothedProfiles{exp_index}, [], 'all', 'omitnan');
MeanSmoothedProfiles{exp_index} = MeanSmoothedProfiles{exp_index}/maxProf;
SmoothedProfileSEs{exp_index} = SmoothedProfileSEs{exp_index}/maxProf;
end

%%
exp_index = 3;
for i = 5:37
    figure(exp_index)
scatter(x_sample, ys(:,i))
hold on
% for rep = 1:NumBootstrappedFits
%     plot(xfits, SmoothedProfiles(:,13,rep) )
% end
plot(xfits, MeanSmoothedProfiles{exp_index}(:,i), 'k', 'LineWidth', 2)
xlabel('Dubuis Time (m)')
title(['AP: ', num2str(APbins(i))])
hold off
end

%%

outpath = 'S:/Gabriella/Dropbox/BootstrappedTestData/Stricter25CBootstrappedGaussianSmoothedProfiles.mat';
save(outpath, 'MeanSmoothedProfiles', 'SmoothedProfileSEs', 'xfits');


%%
MeanSmoothedProfiles = {};
SmoothedProfileSEs = {};
NormedSmoothProfiles = {};
NormedSmoothProfilesSE = {};

for exp_index = 1:3

UseTF = TestSetTFs{exp_index} & NC14TFs{exp_index} & ~isnan(x_samples{exp_index});
x_sample = x_samples{exp_index}(UseTF);
ys = FullFitBcdProfs2{exp_index}(UseTF, :);
NumAPbins = length(APbins);
BadTF =  sum(isnan(ys), 2).' >min( sum(isnan(ys), 2).');
x_sample = x_sample(~BadTF);
ys = ys(~BadTF,:);
min_2sigma_points = 10;
sigma = 5;
window_multiplier = 3;
NumBootstrappedFits = 1000;
NumPoints = 100;


[x_sample, sort_order] = sort(x_sample);
ys = ys(sort_order,:);
xfits = 0:.1:70;
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
    
    if counts_above(x_index) < 5
        GaussianWeights(x_index,:) = 0;
    end
    
    if counts_below(x_index) < 5
        GaussianWeights(x_index,:) = 0;
    end

end


ValidCountTFs = (counts >= min_2sigma_points & counts_above >= 3 & counts_below >= 3 & xfits >= min(x_sample)+sigma/2 & xfits <= max(x_sample) - sigma/2);
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



MeanSmoothedProfiles{exp_index} = mean(SmoothedProfiles, 3, 'omitnan');
SmoothedProfileSEs{exp_index} = std(SmoothedProfiles,0,  3, 'omitnan')/sqrt(NumBootstrappedFits);

SmoothedProfileSEs{exp_index}(MeanSmoothedProfiles{exp_index}  == 0) = NaN;
MeanSmoothedProfiles{exp_index}(MeanSmoothedProfiles{exp_index}  == 0) = NaN;
%close all
figure(exp_index)
scatter(x_sample, ys(:,13))
hold on
% for rep = 1:NumBootstrappedFits
%     plot(xfits, SmoothedProfiles(:,13,rep) )
% end
plot(xfits, MeanSmoothedProfiles{exp_index}(:,13), 'k', 'LineWidth', 2)
hold off

maxProf = max(MeanSmoothedProfiles{exp_index}, [], 'all', 'omitnan');
minProf = min(MeanSmoothedProfiles{exp_index}, [], 'all', 'omitnan');
MeanSmoothedProfiles{exp_index} = MeanSmoothedProfiles{exp_index};
SmoothedProfileSEs{exp_index} = SmoothedProfileSEs{exp_index};
NormedMeanSmoothedProfiles{exp_index} = (MeanSmoothedProfiles{exp_index}-minProf)/(maxProf-minProf);
NormedSmoothedProfileSEs{exp_index} = SmoothedProfileSEs{exp_index}/(maxProf-minProf);
end

%%
exp_index = 3;
for i = 5:37
    figure(exp_index)
scatter(x_sample, ys(:,i))
hold on
% for rep = 1:NumBootstrappedFits
%     plot(xfits, SmoothedProfiles(:,13,rep) )
% end
plot(xfits, MeanSmoothedProfiles{exp_index}(:,i), 'k', 'LineWidth', 2)
xlabel('Dubuis Time (m)')
title(['AP: ', num2str(APbins(i))])
hold off
end

%%

outpath = 'S:/Gabriella/Dropbox/BootstrappedTestData/25CBootstrappedGaussianSmoothedProfiles.mat';
save(outpath, 'MeanSmoothedProfiles', 'SmoothedProfileSEs', 'xfits', 'NormedSmoothedProfileSEs', 'NormedMeanSmoothedProfiles');


