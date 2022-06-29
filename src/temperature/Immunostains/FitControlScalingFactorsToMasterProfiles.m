function CompiledEmbryos = FitControlScalingFactorsToMasterProfiles(CompiledEmbryos, exp_index)
%%
NChannels = 5;
AllSetInfo = GetFixedSetPrefixInfo;

AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';

SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];

outpath = 'S:/Gabriella/Dropbox/BootstrappedTestData/25CBootstrappedGaussianSmoothedProfiles.mat';
load(outpath, 'MeanSmoothedProfiles', 'SmoothedProfileSEs', 'xfits');

%%
ch_index = 3;
ControlSetTF = CompiledEmbryos.IsNC14 & CompiledEmbryos.ControlSetEmbryos & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);


CompiledEmbryos.mdl = cell(3, NChannels);
CompiledEmbryos.BootstrappedScaleFactors = NaN(3, NChannels);
CompiledEmbryos.BootstrappedScaleFactorSEs = NaN(3, NChannels);
CompiledEmbryos.BootstrappedScaleIntercepts = NaN(3, NChannels);
CompiledEmbryos.BootstrappedScaleInterceptSEs = NaN(3, NChannels);

f = @(b,x) b(1).*x + b(2);

for ProfIndex = 1:3
MasterProfile = MeanSmoothedProfiles{ProfIndex};
xsample = CompiledEmbryos.DubuisEmbryoTimes(ControlSetTF);

ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);

yset = CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles(ControlSetTF,:,ch_index);

beta0 =[1, min(yset,[], 'all')];
    if ~isnan(beta0(2))
    CompiledEmbryos.mdl{ProfIndex, ch_index} = fitnlm(ysample(:),yset(:),f,beta0);
    CompiledEmbryos.BootstrappedScaleFactors(ProfIndex, ch_index) = 1/CompiledEmbryos.mdl{ProfIndex, ch_index}.Coefficients.Estimate(1);
    CompiledEmbryos.BootstrappedScaleIntercepts(ProfIndex, ch_index) = -CompiledEmbryos.mdl{ProfIndex, ch_index}.Coefficients.Estimate(2)/CompiledEmbryos.mdl{ProfIndex, ch_index}.Coefficients.Estimate(1);
    CompiledEmbryos.BootstrappedScaleFactorSEs(ProfIndex, ch_index) = CompiledEmbryos.mdl{ProfIndex, ch_index}.Coefficients.SE(1)/(CompiledEmbryos.mdl{ProfIndex, ch_index}.Coefficients.Estimate(1)^2);
    CompiledEmbryos.BootstrappedScaleInterceptSEs(ProfIndex, ch_index) = sqrt(CompiledEmbryos.mdl{ProfIndex, ch_index}.Coefficients.SE(1)^2*CompiledEmbryos.mdl{ProfIndex, ch_index}.Coefficients.Estimate(2)^2/(CompiledEmbryos.mdl{ProfIndex, ch_index}.Coefficients.Estimate(1)^4)+...
        CompiledEmbryos.mdl{ProfIndex, ch_index}.Coefficients.SE(2)^2/(CompiledEmbryos.mdl{ProfIndex, ch_index}.Coefficients.Estimate(1)^2));
    else 
        CompiledEmbryos.mdl{ProfIndex, ch_index} = {};
    end
end
    
    
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
save(CEoutpath, 'CompiledEmbryos');