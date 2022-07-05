function CompiledEmbryos =  AddBootstrappedTestFitProfiles(CompiledEmbryos, exp_index)
AllSetInfo = GetFixedSetPrefixInfo;

AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';

SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];


[NumEmbryos, NumAPbins, NChannels] = size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles);




%%

xfits = 0:.1:70;
Nxfits = length(xfits);
min_2sigma_points = 5;
counts_above_limit = 3;
counts_below_limit = 3;
sigma = 5;
window_multiplier = 3;
NumBootstrappedFits = 200;
NumPoints = 100;
NumEmbryos = size(CompiledEmbryos.UnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1, 1);
NumAPbins = size(CompiledEmbryos.UnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1, 2);
NChannels =size(CompiledEmbryos.UnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1, 3);


%% Bootstrapping the FitSlideRescaledDorsalAvgAPProfiles for the Test Set]
if ~isfield(CompiledEmbryos, 'BootstrappedUnivScaledProfiles')
    CompiledEmbryos.BootstrappedUnivScaledProfiles = {};
end
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1 = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1.x = xfits;
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1.TestSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1.ControlSet.counts =NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1.TestSet.counts = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1.ControlSet.counts_above =NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1.TestSet.counts_above = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1.ControlSet.counts_below = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1.TestSet.counts_below = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2 = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.counts =NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.counts = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.counts_above =NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.counts_above = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.counts_below = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.counts_below = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped.TestSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped.ControlSet.counts =NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped.TestSet.counts = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped.ControlSet.counts_above =NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped.TestSet.counts_above = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped.ControlSet.counts_below = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped.TestSet.counts_below = NaN(Nxfits, NumAPbins, NChannels);


CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1 = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1.x = xfits;
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1.TestSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1.ControlSet.counts =NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1.TestSet.counts = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1.ControlSet.counts_above =NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1.TestSet.counts_above = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1.ControlSet.counts_below = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1.TestSet.counts_below = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2 = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.counts =NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.counts = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.counts_above =NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.counts_above = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.counts_below = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.counts_below = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.TestSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.ControlSet.counts =NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.TestSet.counts = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.ControlSet.counts_above =NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.TestSet.counts_above = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.ControlSet.counts_below = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.TestSet.counts_below = NaN(Nxfits, NumAPbins, NChannels);

CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1 = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1.x = xfits;
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1.TestSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1.ControlSet.counts =NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1.TestSet.counts = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1.ControlSet.counts_above =NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1.TestSet.counts_above = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1.ControlSet.counts_below = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1.TestSet.counts_below = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2 = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.counts =NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.counts = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.counts_above =NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.counts_above = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.counts_below = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.counts_below = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped.TestSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped.ControlSet.counts =NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped.TestSet.counts = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped.ControlSet.counts_above =NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped.TestSet.counts_above = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped.ControlSet.counts_below = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped.TestSet.counts_below = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles = {};

CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep1 = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep1.x = xfits;
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep1.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep1.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep1.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep1.TestSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep1.ControlSet.counts =NaN(1, Nxfits);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep1.TestSet.counts =NaN(1, Nxfits);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep1.ControlSet.counts_above =NaN(1, Nxfits);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep1.TestSet.counts_above = NaN(1, Nxfits);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep1.ControlSet.counts_below = NaN(1, Nxfits);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep1.TestSet.counts_below = NaN(1, Nxfits);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep2 = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep2.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep2.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep2.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep2.TestSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep2.ControlSet.counts =NaN(1, Nxfits);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep2.TestSet.counts = NaN(1, Nxfits);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep2.ControlSet.counts_above =NaN(1, Nxfits);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep2.TestSet.counts_above =NaN(1, Nxfits);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep2.ControlSet.counts_below =NaN(1, Nxfits);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Rep2.TestSet.counts_below = NaN(1, Nxfits);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Flipped = {};
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Flipped.ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Flipped.TestSet.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Flipped.ControlSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Flipped.TestSet.se = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Flipped.ControlSet.counts =NaN(1, Nxfits);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Flipped.TestSet.counts = NaN(1, Nxfits);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Flipped.ControlSet.counts_above =NaN(1, Nxfits);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Flipped.TestSet.counts_above = NaN(1, Nxfits);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Flipped.ControlSet.counts_below = NaN(1, Nxfits);
CompiledEmbryos.BootstrappedUnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.Flipped.TestSet.counts_below = NaN(1, Nxfits);

ProfString1s = {'TestSetSlideRescaledDorsalAvgAPProfiles', 'TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles',...
    'TestSetFitSlideRescaledDorsalAvgAPProfiles','HybridTestSetRescaledDorsalAvgAPProfiles',...
    'TestSetSlideRescaledDorsalAvgAPProfiles', 'TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles',...
    'TestSetFitSlideRescaledDorsalAvgAPProfiles', 'HybridTestSetRescaledDorsalAvgAPProfiles',...
    'TestSetSlideRescaledDorsalAvgAPProfiles', 'TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles',...
    'TestSetFitSlideRescaledDorsalAvgAPProfiles', 'HybridTestSetRescaledDorsalAvgAPProfiles'};

ProfString2s = {'Rep1', 'Rep1',...
    'Rep1','Rep1',...
    'Rep2', 'Rep2',...
    'Rep2','Rep2',...
    'Flipped', 'Flipped',...
    'Flipped','Flipped'};


ProfString3s = {'SlideRescaledDorsalAvgAP', 'ZeroCorrectedSlideRescaledDorsalAvgAP',...
    'FitSlideRescaledDorsalAvgAP','ZeroCorrectedSlideRescaledDorsalAvgAP',...
    'SlideRescaledDorsalAvgAP', 'ZeroCorrectedSlideRescaledDorsalAvgAP',...
    'FitSlideRescaledDorsalAvgAP', 'ZeroCorrectedSlideRescaledDorsalAvgAP',...
    'SlideRescaledDorsalAvgAP', 'ZeroCorrectedSlideRescaledDorsalAvgAP',...
    'FitSlideRescaledDorsalAvgAP', 'ZeroCorrectedSlideRescaledDorsalAvgAP'};


NumTypes = length(ProfString1s);
IsFlipped = false;
IsRep1 = false;
IsRep2 = false;
if mod(exp_index, 3) == 0
    IsFlipped = true;
    Rep1 = exp_index-2;
    Rep2 = exp_index-1;
    Flipped = exp_index;
elseif mod(exp_index, 3) == 1
    IsRep1 = true;
    Rep1 = exp_index;
    Rep2 = exp_index + 1;
    Flipped = exp_index + 2;
else
    IsRep2 = true;
    Rep1 = exp_index-1;
    Rep2 = exp_index;
    Flipped = exp_index+1;
end

%%
UseTestTF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);
UseControlTF = CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);

for type_index = 1:NumTypes
    ProfString1 = ProfString1s{type_index};
    ProfString2 = ProfString2s{type_index};
    ProfString3 = ProfString3s{type_index};
    if (IsRep1 & type_index < 5) | (IsRep2 & type_index >= 5 & type_index < 9) | (IsFlipped & type_index >= 9)
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).x = ...
            CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.x;
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).ControlSet.mean(:,:,:) = ...
            CompiledEmbryos.BootstrappedProfiles.(ProfString3).Control.mean;
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).ControlSet.se(:,:,:) = ...
            CompiledEmbryos.BootstrappedProfiles.(ProfString3).Control.se;
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).ControlSet.counts = ...
            CompiledEmbryos.BootstrappedProfiles.(ProfString3).Control.counts;
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).ControlSet.counts_above = ...
            CompiledEmbryos.BootstrappedProfiles.(ProfString3).Control.counts_above;
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).ControlSet.counts_below = ...
            CompiledEmbryos.BootstrappedProfiles.(ProfString3).Control.counts_below;
        
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).TestSet.mean(:,:,:) = ...
            CompiledEmbryos.BootstrappedProfiles.(ProfString3).Test.mean;
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).TestSet.se(:,:,:) = ...
            CompiledEmbryos.BootstrappedProfiles.(ProfString3).Test.se;
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).TestSet.counts = ...
            CompiledEmbryos.BootstrappedProfiles.(ProfString3).Test.counts;
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).TestSet.counts_above = ...
            CompiledEmbryos.BootstrappedProfiles.(ProfString3).Test.counts_above;
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).TestSet.counts_below = ...
            CompiledEmbryos.BootstrappedProfiles.(ProfString3).Test.counts_below;

        continue
    end
    
    for ch_index = [3 5]
        disp(['Type: ', num2str(type_index), ', Channel: ', num2str(ch_index)]);
        x_sample = CompiledEmbryos.DubuisEmbryoTimes(UseTestTF);
        ys = CompiledEmbryos.UnivScaledProfiles.(ProfString1).(ProfString2);
        
        ys = ys(UseTestTF, :,ch_index);
        BadTF =  sum(isnan(ys), 2).' > NumAPbins/2;
        
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
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).TestSet.mean(:,:,ch_index) = MeanProf;
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).TestSet.se(:,:,ch_index) = SEProf;
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).TestSet.counts = counts;
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).TestSet.counts_above = counts_above;
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).TestSet.counts_below = counts_below;
        
        % Control Set Smoothing
        x_sample = CompiledEmbryos.DubuisEmbryoTimes(UseControlTF);
        ys = CompiledEmbryos.UnivScaledProfiles.(ProfString1).(ProfString2);
        
        ys = ys(UseControlTF, :,ch_index);
        BadTF =  sum(isnan(ys), 2).' >NumAPbins/2;
        
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
        
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).ControlSet.mean(:,:,ch_index) = MeanProf;
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).ControlSet.se(:,:,ch_index) = SEProf;
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).ControlSet.counts = counts;
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).ControlSet.counts_above = counts_above;
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(ProfString1).(ProfString2).ControlSet.counts_below = counts_below;
        
        
    end
    
end
%%
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
save(CEoutpath, 'CompiledEmbryos');