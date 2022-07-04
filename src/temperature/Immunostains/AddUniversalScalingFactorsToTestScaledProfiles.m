function CompiledEmbryos =  AddUniversalScalingFactorsToTestScaledProfiles(CompiledEmbryos, exp_index)
AllSetInfo = GetFixedSetPrefixInfo;

AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';

SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];


[NumEmbryos, NumAPbins, NChannels] = size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles);
NumMasterProfSets = size(CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Control, 1);




%%

xfits = 0:1:70;
Nxfits = length(xfits);
min_2sigma_points = 5;
counts_above_limit = 3;
counts_below_limit = 3;
sigma = 5;
window_multiplier = 3;
NumBootstrappedFits = 200;
NumPoints = 100;
NumEmbryos = size(CompiledEmbryos.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling, 1);
NumAPbins = size(CompiledEmbryos.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling, 2);
NChannels =size(CompiledEmbryos.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling, 3);
NumMasterProfSets = size(CompiledEmbryos.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling, 4);


%% Bootstrapping the FitSlideRescaledDorsalAvgAPProfiles for the Test Set
CompiledEmbryos.BootstrappedUnivScaledProfiles = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x = xfits;
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.se = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.counts = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.counts = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.counts_above = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.counts_above = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.counts_below = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.counts_below = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.TestScaling = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.TestScaling.x = xfits;
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.se = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.counts = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.counts = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.counts_above = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.counts_above = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.counts_below = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.counts_below = NaN(Nxfits, NChannels, NumMasterProfSets);


CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.x = xfits;
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.se = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.counts = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.counts = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.counts_above = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.counts_above = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.counts_below = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.counts_below = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling.x = xfits;
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.se = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.counts = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.counts = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.counts_above = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.counts_above = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.counts_below = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.counts_below = NaN(Nxfits, NChannels, NumMasterProfSets);


CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.x = xfits;
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.se = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.counts = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.counts = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.counts_above = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.counts_above = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.counts_below = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.counts_below = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling.x = xfits;
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.se = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.counts = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.counts = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.counts_above = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.counts_above = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.counts_below = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.counts_below = NaN(Nxfits, NChannels, NumMasterProfSets);


CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.ControlScaling = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.ControlScaling.x = xfits;
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.se = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.counts = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.counts = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.counts_above = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.counts_above = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.counts_below = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.counts_below = NaN(Nxfits, NChannels, NumMasterProfSets);

CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.TestScaling = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.TestScaling.x = xfits;
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.TestScaling.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.TestScaling.TestSet.se = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.counts = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.TestScaling.TestSet.counts = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.counts_above = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.TestScaling.TestSet.counts_above = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.counts_below = NaN(Nxfits, NChannels, NumMasterProfSets);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridRescaledDorsalAvgAPProfiles.TestScaling.TestSet.counts_below = NaN(Nxfits, NChannels, NumMasterProfSets);


ProfString1s = {'SlideRescaledDorsalAvgAPProfiles', 'ZeroCorrectedSlideRescaledDorsalAvgAPProfiles',...
    'FitSlideRescaledDorsalAvgAPProfiles','HybridRescaledDorsalAvgAPProfiles',...
    'SlideRescaledDorsalAvgAPProfiles', 'ZeroCorrectedSlideRescaledDorsalAvgAPProfiles',...
    'FitSlideRescaledDorsalAvgAPProfiles', 'HybridRescaledDorsalAvgAPProfiles'};

ProfString2s = {'ControlScaling', 'ControlScaling',...
    'ControlScaling','ControlScaling',...
    'TestScaling', 'TestScaling',...
    'TestScaling', 'TestScaling'};

if AllSetInfo.Temperatures(exp_index) == 25
    NumTypes = length(ProfString1s);
else
    NumTypes = length(ProfString1s)/2;
end

%%
UseTestTF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);
UseControlTF = CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);

for type_index = 1:NumTypes
    ProfString1 = ProfString1s{type_index};
    ProfString2 = ProfString2s{type_index};
    
    
    for set_index = 1:NumMasterProfSets
        for ch_index = [3 5]
            disp(['Type: ', num2str(type_index), ', Master Set: ', num2str(set_index), ', Channel: ', num2str(ch_index)]);
            x_sample = CompiledEmbryos.DubuisEmbryoTimes(UseTestTF);
            ys = CompiledEmbryos.UnivScaledProfiles.(ProfString1).(ProfString2);
            
            ys = ys(UseTestTF, :,ch_index, set_index);
            BadTF =  sum(isnan(ys), 2).' >min( sum(isnan(ys), 2).');
            
            x_sample = x_sample(~BadTF);
            ys = ys(~BadTF,:);
            [x_sample, sort_order] = sort(x_sample);
            ys = ys(sort_order,:);
            
            
            DiffMat = xfits.'-x_sample;
            GaussianWeights = GetGaussianWeightMat(xfits, x_sample, sigma, window_multiplier*sigma);
            
            SmoothedProfiles = NaN(Nxfits, NumAPbins, NumBootstrappedFits);
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
            
            CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).x = xfits;
            CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).TestSet.mean(:,:,ch_index, set_index) = MeanProf;
            CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).TestSet.se(:,:,ch_index, set_index) = SEProf;
            CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).TestSet.counts(:,ch_index, set_index) = counts;
            CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).TestSet.counts_above(:,ch_index, set_index) = counts_above;
            CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).TestSet.counts_below(:,ch_index, set_index) = counts_below;
            
            % Control Set Smoothing
            x_sample = CompiledEmbryos.DubuisEmbryoTimes(UseControlTF);
            ys = CompiledEmbryos.UnivScaledProfiles.(ProfString1).(ProfString2);
            
            ys = ys(UseControlTF, :,ch_index, set_index);
            BadTF =  sum(isnan(ys), 2).' >min( sum(isnan(ys), 2).');
            
            x_sample = x_sample(~BadTF);
            ys = ys(~BadTF,:);
            [x_sample, sort_order] = sort(x_sample);
            ys = ys(sort_order,:);
            
            
            DiffMat = xfits.'-x_sample;
            GaussianWeights = GetGaussianWeightMat(xfits, x_sample, sigma, window_multiplier*sigma);
            
            SmoothedProfiles = NaN(Nxfits, NumAPbins, NumBootstrappedFits);
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
            
            CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).ControlSet.mean(:,:,ch_index, set_index) = MeanProf;
            CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).ControlSet.se(:,:,ch_index, set_index) = SEProf;
            CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).ControlSet.counts(:,ch_index, set_index) = counts;
            CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).ControlSet.counts_above(:,ch_index, set_index) = counts_above;
            CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).ControlSet.counts_below(:,ch_index, set_index) = counts_below;
            
            
        end

    end
end
%%
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
save(CEoutpath, 'CompiledEmbryos');