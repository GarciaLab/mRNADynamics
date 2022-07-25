function AllCompiledEmbryos = FitScalingFactorsToTestProfiles(AllCompiledEmbryos)
%%
NChannels = 5;
AllSetInfo = GetFixedSetPrefixInfo;

AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';




f = @(b,x) b(1).*x + b(2);

NumSets = length(AllSetInfo.Temperatures);

%% Fit to Slide Rescal Dorsal Avg AP Profiles
for exp_index = 1:NumSets
    CompiledEmbryos = AllCompiledEmbryos{exp_index};
    SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];

if ~isfield(CompiledEmbryos, 'mdl')
    CompiledEmbryos.mdl = {};
end
if ~isfield(CompiledEmbryos, 'ScaleFits')
    CompiledEmbryos.ScaleFits  = {};
end

ControlSetTF = CompiledEmbryos.IsNC14 & CompiledEmbryos.ControlSetEmbryos & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);
TestSetTF = CompiledEmbryos.IsNC14 & CompiledEmbryos.TestSetEmbryos & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);

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
    
SetIndices = [Rep1, Rep2, Flipped];
SetStrings = {'Rep1', 'Rep2', 'Flipped'};
%% Fit to Exponential-Fit Bicoid Profiles

CompiledEmbryos.mdl.TestSetFitSlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles = {};

CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1 = {};
CompiledEmbryos.mdl.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1 = cell(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1 = {};
CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1.ScaleEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1.ScaleSE = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1.InterceptEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1.InterceptSE = NaN(1, NChannels);

CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2 = {};
CompiledEmbryos.mdl.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2 = cell(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2 = {};
CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2.ScaleEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2.ScaleSE = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2.InterceptEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2.InterceptSE = NaN(1, NChannels);


CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped = {};
CompiledEmbryos.mdl.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped = cell(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped = {};
CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped.ScaleEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped.ScaleSE = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped.InterceptEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.Flipped.InterceptSE = NaN(1, NChannels);




for ch_index = [3]



for master_index = 1:3
    SetString = SetStrings{master_index};
    if master_index == 1 & IsRep1
        continue
    elseif master_index == 2 & IsRep2
        continue
    elseif master_index == 3 & IsFlipped
        continue
    end
    MasterProfile = AllCompiledEmbryos{SetIndices(master_index)}.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Test.mean(:,:,ch_index);
    xfits =  AllCompiledEmbryos{SetIndices(master_index)}.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.x;
    xsample = CompiledEmbryos.DubuisEmbryoTimes(TestSetTF);
    
    ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
    
    yset = CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles(TestSetTF,:,ch_index);
    
    beta0 =[1, min(yset,[], 'all')];
    if ~isnan(beta0(2))
        dlm = fitnlm(ysample(:),yset(:),f,beta0);
        CompiledEmbryos.mdl.TestSetFitSlideRescaledDorsalAvgAPProfiles.(SetString){ch_index} = dlm;
        CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.(SetString).ScaleEstimate(ch_index) = 1/dlm.Coefficients.Estimate(1);
        CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.(SetString).InterceptEstimate(ch_index) = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
        CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.(SetString).ScaleSE(ch_index) = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
        CompiledEmbryos.ScaleFits.TestSetFitSlideRescaledDorsalAvgAPProfiles.(SetString).InterceptSE(ch_index) = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
            dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
    end
    

end

end

%% Fit to Exponential-Fit Zero-Corrected Bicoid Profiles

CompiledEmbryos.mdl.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles = {};

CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep1 = {};
CompiledEmbryos.mdl.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep1 = cell(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep1 = {};
CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep1.ScaleEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep1.ScaleSE = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep1.InterceptEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep1.InterceptSE = NaN(1, NChannels);

CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep2 = {};
CompiledEmbryos.mdl.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep2 = cell(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep2 = {};
CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep2.ScaleEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep2.ScaleSE = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep2.InterceptEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Rep2.InterceptSE = NaN(1, NChannels);


CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Flipped = {};
CompiledEmbryos.mdl.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Flipped = cell(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Flipped = {};
CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Flipped.ScaleEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Flipped.ScaleSE = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Flipped.InterceptEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.Flipped.InterceptSE = NaN(1, NChannels);




for ch_index = [3]



for master_index = 1:3
    SetString = SetStrings{master_index};
    if master_index == 1 & IsRep1
        continue
    elseif master_index == 2 & IsRep2
        continue
    elseif master_index == 3 & IsFlipped
        continue
    end
    MasterProfile = AllCompiledEmbryos{SetIndices(master_index)}.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Test.mean(:,:,ch_index);
    xfits =  AllCompiledEmbryos{SetIndices(master_index)}.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.x;
    xsample = CompiledEmbryos.DubuisEmbryoTimes(TestSetTF);
    
    ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
    
    yset = CompiledEmbryos.FitZeroedSlideRescaledDorsalAvgAPProfiles(TestSetTF,:,ch_index);
    
    beta0 =[1, min(yset,[], 'all')];
    if ~isnan(beta0(2))
        dlm = fitnlm(ysample(:),yset(:),f,beta0);
        CompiledEmbryos.mdl.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.(SetString){ch_index} = dlm;
        CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.(SetString).ScaleEstimate(ch_index) = 1/dlm.Coefficients.Estimate(1);
        CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.(SetString).InterceptEstimate(ch_index) = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
        CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.(SetString).ScaleSE(ch_index) = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
        CompiledEmbryos.ScaleFits.TestSetFitZeroedSlideRescaledDorsalAvgAPProfiles.(SetString).InterceptSE(ch_index) = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
            dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
    end
    

end

end

%% Fit to Zero Corrected Slide Rescaled Dorsal Avg AP Profiles

CompiledEmbryos.mdl.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles = {};

CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1 = {};
CompiledEmbryos.mdl.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1 = cell(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1 = {};
CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1.ScaleEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1.ScaleSE = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1.InterceptEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep1.InterceptSE = NaN(1, NChannels);

CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2 = {};
CompiledEmbryos.mdl.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2 = cell(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2 = {};
CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.ScaleEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.ScaleSE = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.InterceptEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.InterceptSE = NaN(1, NChannels);


CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped = {};
CompiledEmbryos.mdl.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped = cell(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped = {};
CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.ScaleEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.ScaleSE = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.InterceptEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.InterceptSE = NaN(1, NChannels);




for ch_index = [3 5]
for master_index = 1:3
    SetString = SetStrings{master_index};
    if master_index == 1 & IsRep1
        continue
    elseif master_index == 2 & IsRep2
        continue
    elseif master_index == 3 & IsFlipped
        continue
    end
    MasterProfile = AllCompiledEmbryos{SetIndices(master_index)}.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.Test.mean(:,:,ch_index);
    xfits =  AllCompiledEmbryos{SetIndices(master_index)}.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.x;
    xsample = CompiledEmbryos.DubuisEmbryoTimes(TestSetTF);
    
    ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
    
    yset = CompiledEmbryos.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(TestSetTF,:,ch_index);
    
    beta0 =[1, min(yset,[], 'all')];
    if ~isnan(beta0(2))
        dlm = fitnlm(ysample(:),yset(:),f,beta0);
        CompiledEmbryos.mdl.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.(SetString){ch_index} = dlm;
        CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.(SetString).ScaleEstimate(ch_index) = 1/dlm.Coefficients.Estimate(1);
        CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.(SetString).InterceptEstimate(ch_index) = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
        CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.(SetString).ScaleSE(ch_index) = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
        CompiledEmbryos.ScaleFits.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.(SetString).InterceptSE(ch_index) = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
            dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
    end
    

end

end


%% Fit to Slide Rescaled Dorsal Average AP Corrected Slide Rescaled Dorsal Avg AP Profiles

CompiledEmbryos.mdl.TestSetSlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles = {};

CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1 = {};
CompiledEmbryos.mdl.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1 = cell(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1 = {};
CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1.ScaleEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1.ScaleSE = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1.InterceptEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Rep1.InterceptSE = NaN(1, NChannels);

CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2 = {};
CompiledEmbryos.mdl.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2 = cell(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2 = {};
CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2.ScaleEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2.ScaleSE = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2.InterceptEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Rep2.InterceptSE = NaN(1, NChannels);


CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped = {};
CompiledEmbryos.mdl.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped = cell(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped = {};
CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped.ScaleEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped.ScaleSE = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped.InterceptEstimate = NaN(1, NChannels);
CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.Flipped.InterceptSE = NaN(1, NChannels);




for ch_index = [3 5]
for master_index = 1:3
    SetString = SetStrings{master_index};
    if master_index == 1 & IsRep1
        continue
    elseif master_index == 2 & IsRep2
        continue
    elseif master_index == 3 & IsFlipped
        continue
    end
    MasterProfile = AllCompiledEmbryos{SetIndices(master_index)}.BootstrappedProfiles.SlideRescaledDorsalAvgAP.Test.mean(:,:,ch_index);
    xfits =  AllCompiledEmbryos{SetIndices(master_index)}.BootstrappedProfiles.SlideRescaledDorsalAvgAP.x;
    xsample = CompiledEmbryos.DubuisEmbryoTimes(TestSetTF);
    
    ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
    
    yset = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(TestSetTF,:,ch_index);
    
    beta0 =[1, min(yset,[], 'all')];
    if ~isnan(beta0(2))
        dlm = fitnlm(ysample(:),yset(:),f,beta0);
        CompiledEmbryos.mdl.TestSetSlideRescaledDorsalAvgAPProfiles.(SetString){ch_index} = dlm;
        CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.(SetString).ScaleEstimate(ch_index) = 1/dlm.Coefficients.Estimate(1);
        CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.(SetString).InterceptEstimate(ch_index) = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
        CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.(SetString).ScaleSE(ch_index) = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
        CompiledEmbryos.ScaleFits.TestSetSlideRescaledDorsalAvgAPProfiles.(SetString).InterceptSE(ch_index) = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
            dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
    end
    

end

end
%%
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
save(CEoutpath, 'CompiledEmbryos');
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
end
%%
