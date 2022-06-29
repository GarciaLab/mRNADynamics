% data$weight <- (abs(data$Time-(data$group_num*0.05+0.025)))^(-1)
% 
% data$Norm<-ave(data$weight,data$group_num,FUN=function(x) x/sum(x))
% 
% data$Time2<- data$Time*data$Norm
% 
% data$Right2<- data$Right*data$Norm
% 
% data$Left2<- data$Left*data$Norm
% 
% data2$Time<- tapply(data$Time2, data$group_num, sum)
% 
% data2$Right<- tapply(data$Right2, data$group_num, sum)
% 
% data2$Left<- tapply(data$Left2, data$group_num, sum)

load('S:/Gabriella/Dropbox/Petkova2019Data/data_raw_dorsal_151220_Cic_Bcd_Hb.mat');


f = @(b,x) b(1).*exp(b(2).*x) + b(3);
NumObs = length(data);
fitx = (.5:999.5)/1000;
mdl = {};
FitBcdProfs = NaN(NumObs,1000);
RawBcdProfs = NaN(NumObs, 1000);
for i = 1:NumObs
    beta0 =[ max(data(i).Bcd), -3, min(data(i).Bcd)];
    
    mdl{i} = fitnlm(fitx(201:900),data(i).Bcd(201:900),f,beta0);
    FitBcdProfs(i,201:900) = mdl{i}.Fitted.';
    RawBcdProfs(i,:) = data(i).Bcd;
end



%%
AllSetsCombinedEmbryosPath = ['S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/'];
AllSetInfo = GetFixedSetPrefixInfo;
 
AllCompiledEmbryos = cell(1, 3);
for exp_index = 1:15
SetLabel = AllSetInfo.SetLabels{exp_index};
PlotLabel = AllSetInfo.PlotLabels{exp_index};
SetPrefixes = AllSetInfo.Prefixes{exp_index};
SetIsFlipped = AllSetInfo.Flipped(exp_index);
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];

if ~isdir(OutEmbryoPath)
    mkdir(OutEmbryoPath)
end
liveExperiments = cell(1, length(SetPrefixes));
for i = 1:length(SetPrefixes)
    liveExperiments{i} = LiveExperiment(SetPrefixes{i});
end
FixedPixelSize_um = liveExperiments{1}.pixelSize_um;
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
%load(CEoutpath, 'CompiledEmbryos');
CompiledEmbryos = CombineCompiledEmbryos(SetPrefixes);

CompiledEmbryos = AddEmbryoTimingInfo(CompiledEmbryos, exp_index);

NEmbryos = length(CompiledEmbryos.Approved);
AllEmbryos = 1:NEmbryos;
NC14Indices = find(CompiledEmbryos.IsNC14);
NC13Indices = find(CompiledEmbryos.IsNC13);
NC13NC14Indices = find(CompiledEmbryos.IsNC13orNC14);
NumEmbryosNC14 = CompiledEmbryos.NumEmbryosNC14;

KnirpsIndex = 4;
CompiledEmbryos = PartitionEmbryosTestControl(CompiledEmbryos, exp_index, KnirpsIndex);


CompiledEmbryos = RescaleSlideFluos(CompiledEmbryos, exp_index);
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
end

%%
fitx2 = 0:0.025:1;
mdl2 = cell(1,4);
FitBcdProfs2 = cell(1,4);
NumObs2 = cell(1,4);
RawBcdProfs2 = cell(1,4);
TFPlot = cell(1,4);
TFPlotTest = cell(1,4);
TFPlotControl = cell(1,4);
for exp_index = 1:3
NumObs2{exp_index}  = size(AllCompiledEmbryos{exp_index}.SlideRescaledDorsalAvgAPProfiles, 1);
FitBcdProfs2{exp_index} = NaN(NumObs2{exp_index},41);
RawBcdProfs2{exp_index} = NaN(NumObs2{exp_index}, 41);
mdl2{exp_index} = {};
for i = 1:NumObs2{exp_index}
    beta0 =[ max(AllCompiledEmbryos{exp_index}.SlideRescaledDorsalAvgAPProfiles(i,:,3)), -3, min(AllCompiledEmbryos{exp_index}.SlideRescaledDorsalAvgAPProfiles(i,:,3))];
    if ~isnan(beta0(1))
    mdl2{exp_index}{i} = fitnlm(fitx2(9:37),AllCompiledEmbryos{exp_index}.SlideRescaledDorsalAvgAPProfiles(i,9:37,3),f,beta0);
    FitBcdProfs2{exp_index}(i,9:37) = mdl2{exp_index}{i}.Fitted.';
    RawBcdProfs2{exp_index}(i,:) =AllCompiledEmbryos{exp_index}.SlideRescaledDorsalAvgAPProfiles(i,:,3);
    else 
        mdl2{exp_index}{i}  = {};
    end
end
TFPlot{exp_index} = AllCompiledEmbryos{exp_index}.IsNC14 & ~isnan(AllCompiledEmbryos{exp_index}.DubuisEmbryoTimes);
TFPlotTest{exp_index} = TFPlot{exp_index} & AllCompiledEmbryos{exp_index}.TestSetEmbryos;
TFPlotControl{exp_index} = TFPlot{exp_index} & AllCompiledEmbryos{exp_index}.ControlSetEmbryos;

end

%%
x_samples = {};
TestSetTFs = {};
ControlSetTFs = {};
NC14TFs = {};
APbins = (0.0125:0.025:(1+0.0125));
FullFitBcdProfs = cell(1, 3);
NormedBcdProfs = cell(1, 3);
for exp_index = 1:3
    FullFitBcdProfs2{exp_index} = NaN(NumObs2{exp_index},size(FitBcdProfs2{exp_index}, 2));
    NormedBcdProfs{exp_index} = RawBcdProfs2{exp_index};
    for i = 1:NumObs2{exp_index}
        if ~isempty(mdl2{exp_index}{i})
            FullFitBcdProfs2{exp_index}(i,:) =FitBcdProfs2{exp_index}(i,:)- mdl2{exp_index}{i}.Coefficients.Estimate(3);
            NormedBcdProfs{exp_index}(i,:) = RawBcdProfs2{exp_index}(i,:) -mdl2{exp_index}{i}.Coefficients.Estimate(3);
        else
            disp(['i = ', num2str(i)])
            disp('pause')
        end
    end
    x_samples{exp_index} = AllCompiledEmbryos{exp_index}.DubuisEmbryoTimes;
    TestSetTFs{exp_index}= AllCompiledEmbryos{exp_index}.TestSetEmbryos;
    ControlSetTFs{exp_index}= AllCompiledEmbryos{exp_index}.ControlSetEmbryos;
    NC14TFs{exp_index} = AllCompiledEmbryos{exp_index}.IsNC14;
    
end

APbins = fitx2;
save('S:/Gabriella/Dropbox/BootstrappedTestData/TestSet.mat', 'APbins', 'FullFitBcdProfs2','NormedBcdProfs', 'RawBcdProfs2', 'mdl2', 'x_samples',...
    'TestSetTFs', 'ControlSetTFs', 'NC14TFs')
%%
SmoothedProfiles = {};
BinData = {};
RefBin = 13;

APbins = 0:0.025:1;
NumAPbins = length(APbins);
NChannels = 5;
sigma = 5;%floor(MaxT*sigma_time_multiplier);

for exp_index = 1:3

 AllCompiledEmbryos{exp_index}.AllDubuisValidProfilesTestTF =  AllCompiledEmbryos{exp_index}.IsNC14 &...
         AllCompiledEmbryos{exp_index}.TestSetEmbryos & ~isnan(AllCompiledEmbryos{exp_index}.DubuisEmbryoTimes);% &...
y =  AllCompiledEmbryos{exp_index}.DubuisEmbryoTimes( AllCompiledEmbryos{exp_index}.AllDubuisValidProfilesTestTF );
% ThreeMinDeltas = mink(y, 3);
% minx = floor(ThreeMinDeltas(end)/.1)*.1;
%MaxT = ceil(max([CompiledEmbryos.DubuisEmbryoTimes(:)]));
MaxT = 70;
x = 0:0.1:MaxT;
DorsalProfiles =  AllCompiledEmbryos{exp_index}.SlideRescaledDorsalAvgAPProfiles;


min_2sigma_points = 10;
DiffMat = x.'-y;
%[nn_idx, nn_D] = knnsearch(y.', x.', 'K', 50);
GaussianWeights = GetGaussianWeightMat(x, y, sigma, 2*sigma);
keepx = ones(1, length(x), 'logical');
% counts_within_halfsigma = zeros(1, length(x));
% counts_within_1sigma = zeros(1, length(x));
% counts_within_2sigma = zeros(1, length(x));
% counts_within_3sigma = zeros(1, length(x));
% counts_within_4sigma = zeros(1, length(x));
% counts_within_2sigma_abovecenter = zeros(1, length(x));
% counts_within_2sigma_belowcenter = zeros(1, length(x));
counts = zeros(1, length(x));
counts_above = zeros(1, length(x));
counts_below = zeros(1, length(x));
for x_index = 1:length(x)
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
    
    if counts_above(x_index) < min_2sigma_points/2
        GaussianWeights(x_index,:) = 0;
    end
    
    if counts_below(x_index) < min_2sigma_points/2
        GaussianWeights(x_index,:) = 0;
    end

end

 [val, idx] = max(sum(GaussianWeights> 0, 2).')
 TestWeights = GaussianWeights(idx,:)
 DiffRow = DiffMat(idx,:)
 
%%


SmoothedProfiles{exp_index} = NaN(length(x), NumAPbins, NChannels);
ValidIndices = find(TFPlotTest{exp_index});
BinData{exp_index} = NaN(sum(TFPlotTest{exp_index}), NumAPbins, NChannels);
for ch_index = [3]%2:NChannels % 2:5
    
for RefBin = 1:NumAPbins
        BinData{exp_index}(:,RefBin,ch_index) = FitBcdProfs2{exp_index}(TFPlotTest{exp_index},RefBin);
        for i = 1:size(BinData{exp_index}, 1)
            if ~isempty(mdl2{exp_index}{ValidIndices(i)})
            BinData{exp_index}(i,RefBin,ch_index)  = BinData{exp_index}(i,RefBin,ch_index) - ...
                mdl2{exp_index}{ValidIndices(i)}.Coefficients.Estimate(3);
            else
              BinData{exp_index}(i,RefBin,ch_index) = NaN;  
            end
        end
        BinIsNaN = isnan(BinData{exp_index}(:,RefBin,ch_index));
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        SubBinData = BinData{exp_index}(~BinIsNaN,RefBin,ch_index);
        SmoothedBinProfile = ((BinGaussianWeights*SubBinData)./WeightSums).';
        SmoothedProfiles{exp_index}(:,RefBin,ch_index) = SmoothedBinProfile.';
end

end


RefBin = 13;
figure(exp_index)
hold off
scatter(AllCompiledEmbryos{exp_index}.DubuisEmbryoTimes(TFPlotTest{exp_index}), BinData{exp_index}(:,RefBin,ch_index), 100, '.')
hold on 
plot(x, SmoothedProfiles{exp_index}(:,RefBin, ch_index))
hold off


end
%%
exp_index = 4;
y =  [data(:).age];
% ThreeMinDeltas = mink(y, 3);
% minx = floor(ThreeMinDeltas(end)/.1)*.1;
%MaxT = ceil(max([CompiledEmbryos.DubuisEmbryoTimes(:)]));
MaxT = 70;
x = 0:0.1:MaxT;


min_2sigma_points = 10;
DiffMat = x.'-y;
%[nn_idx, nn_D] = knnsearch(y.', x.', 'K', 50);
GaussianWeights = GetGaussianWeightMat(x, y, sigma, 4*sigma);
% counts_within_halfsigma = zeros(1, length(x));
% counts_within_1sigma = zeros(1, length(x));
% counts_within_2sigma = zeros(1, length(x));
% counts_within_3sigma = zeros(1, length(x));
% counts_within_4sigma = zeros(1, length(x));
% counts_within_2sigma_abovecenter = zeros(1, length(x));
% counts_within_2sigma_belowcenter = zeros(1, length(x));
counts = zeros(1, length(x));
counts_above = zeros(1, length(x));
counts_below = zeros(1, length(x));
for x_index = 1:length(x)
    Num2SigmaPoints = sum(GaussianWeights(x_index,:)> 0);
    
    if Num2SigmaPoints < min_2sigma_points
        GaussianWeights(x_index,:) = 0;
        continue
    end
    
    counts(x_index) = Num2SigmaPoints;
    MatchedPoints = find(GaussianWeights(x_index,:)> 0);
    
    counts_above(x_index) = sum(DiffMat(x_index,MatchedPoints) > 0);
    counts_below(x_index) = sum(DiffMat(x_index,MatchedPoints) < 0);
    
    if counts_above(x_index) < min_2sigma_points/2
        GaussianWeights(x_index,:) = 0;
    end
    
    if counts_below(x_index) < min_2sigma_points/2
        GaussianWeights(x_index,:) = 0;
    end

end




SmoothedProfiles{exp_index} = NaN(length(x), size(FitBcdProfs, 2), NChannels);
BinData{exp_index} = NaN(size(FitBcdProfs, 1), size(FitBcdProfs, 2), NChannels);
for ch_index = [3]%2:NChannels % 2:5
    
for RefBin = 1:size(FitBcdProfs, 2)
        BinData{exp_index}(:,RefBin,ch_index) = FitBcdProfs(:,RefBin);
        for i = 1:size(BinData{exp_index}, 1)
            if ~isempty(mdl{i})
            BinData{exp_index}(i,RefBin,ch_index)  = BinData{exp_index}(i,RefBin,ch_index) - ...
                mdl{i}.Coefficients.Estimate(3);
            else
              BinData{exp_index}(i,RefBin,ch_index) = NaN;  
            end
        end
        BinIsNaN = isnan(BinData{exp_index}(:,RefBin,ch_index));
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        SubBinData = BinData{exp_index}(~BinIsNaN,RefBin,ch_index);
        SmoothedBinProfile = ((BinGaussianWeights*SubBinData)./WeightSums).';
        SmoothedProfiles{exp_index}(:,RefBin,ch_index) = SmoothedBinProfile.';
end

end


RefBin = 301;
figure(exp_index)
hold off
scatter(y, BinData{exp_index}(:,RefBin,ch_index), 100, '.')
hold on 
plot(x, SmoothedProfiles{exp_index}(:,RefBin, ch_index))
hold off

%%



%%

for exp_index = 5:7

AllDubuisValidProfilesTF =  AllCompiledEmbryos{exp_index-4}.IsNC14  & ~isnan(AllCompiledEmbryos{exp_index-4}.DubuisEmbryoTimes);% &...
y =  AllCompiledEmbryos{exp_index-4}.DubuisEmbryoTimes( AllDubuisValidProfilesTF);
% ThreeMinDeltas = mink(y, 3);
% minx = floor(ThreeMinDeltas(end)/.1)*.1;
%MaxT = ceil(max([CompiledEmbryos.DubuisEmbryoTimes(:)]));
MaxT = 70;
x = 0:0.1:MaxT;
DorsalProfiles =  AllCompiledEmbryos{exp_index-4}.SlideRescaledDorsalAvgAPProfiles;


min_2sigma_points = 10;
DiffMat = x.'-y;
%[nn_idx, nn_D] = knnsearch(y.', x.', 'K', 50);
GaussianWeights = GetGaussianWeightMat(x, y, sigma, 4*sigma);
keepx = ones(1, length(x), 'logical');
% counts_within_halfsigma = zeros(1, length(x));
% counts_within_1sigma = zeros(1, length(x));
% counts_within_2sigma = zeros(1, length(x));
% counts_within_3sigma = zeros(1, length(x));
% counts_within_4sigma = zeros(1, length(x));
% counts_within_2sigma_abovecenter = zeros(1, length(x));
% counts_within_2sigma_belowcenter = zeros(1, length(x));
counts = zeros(1, length(x));
counts_above = zeros(1, length(x));
counts_below = zeros(1, length(x));
for x_index = 1:length(x)
    Num2SigmaPoints = sum(GaussianWeights(x_index,:)> 0);
    
    if Num2SigmaPoints < min_2sigma_points
        GaussianWeights(x_index,:) = 0;
        continue
    end
    
    counts(x_index) = Num2SigmaPoints;
    MatchedPoints = find(GaussianWeights(x_index,:)> 0);
    
    counts_above(x_index) = sum(DiffMat(x_index,MatchedPoints) > 0);
    counts_below(x_index) = sum(DiffMat(x_index,MatchedPoints) < 0);
    
    if counts_above(x_index) < min_2sigma_points/2
        GaussianWeights(x_index,:) = 0;
    end
    
    if counts_below(x_index) < min_2sigma_points/2
        GaussianWeights(x_index,:) = 0;
    end

end




SmoothedProfiles{exp_index} = NaN(length(x), NumAPbins, NChannels);
ValidIndices = find(TFPlot{exp_index-4});
BinData{exp_index} = NaN(sum(AllDubuisValidProfilesTF), NumAPbins, NChannels);
for ch_index = [3]%2:NChannels % 2:5
    
for RefBin = 1:NumAPbins
        BinData{exp_index}(:,RefBin,ch_index) = FitBcdProfs2{exp_index-4}(AllDubuisValidProfilesTF,RefBin);
        for i = 1:size(BinData{exp_index}, 1)
            if ~isempty(mdl2{exp_index-4}{ValidIndices(i)})
            BinData{exp_index}(i,RefBin,ch_index)  = BinData{exp_index}(i,RefBin,ch_index) - ...
                mdl2{exp_index-4}{ValidIndices(i)}.Coefficients.Estimate(3);
            else
              BinData{exp_index}(i,RefBin,ch_index) = NaN;  
            end
        end
        NotBinIsNaN = ~isnan(BinData{exp_index}(:,RefBin,ch_index)).';
        BinGaussianWeights = GaussianWeights(:,NotBinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        SubBinData = BinData{exp_index}(NotBinIsNaN,RefBin,ch_index);
        SmoothedBinProfile = ((BinGaussianWeights*SubBinData)./WeightSums).';
        SmoothedProfiles{exp_index}(:,RefBin,ch_index) = SmoothedBinProfile.';
end

end


RefBin = 13;
figure(exp_index)
hold off
scatter(AllCompiledEmbryos{exp_index-4}.DubuisEmbryoTimes(TFPlot{exp_index-4}), BinData{exp_index}(:,RefBin,ch_index), 100, '.')
hold on 
plot(x, SmoothedProfiles{exp_index}(:,RefBin, ch_index))
hold off


end