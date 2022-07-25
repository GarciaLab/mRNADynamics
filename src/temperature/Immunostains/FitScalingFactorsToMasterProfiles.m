function CompiledEmbryos = FitScalingFactorsToMasterProfiles(CompiledEmbryos, exp_index)
%%
NChannels = 5;
AllSetInfo = GetFixedSetPrefixInfo;

AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';

SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];

BootstrappedProfilePath = 'S:/Gabriella/Dropbox/BootstrappedTestData/25CBootstrappedGaussianSmoothedProfiles.mat';
load(BootstrappedProfilePath, 'MeanSmoothedProfiles', 'SmoothedProfileSEs', 'xfits');

NumMasterProfs = length(MeanSmoothedProfiles);


f = @(b,x) b(1).*x + b(2);

%% Fit to Slide Rescal Dorsal Avg AP Profiles
CompiledEmbryos.mdl = {};
CompiledEmbryos.ScaleFits  = {};

CompiledEmbryos.mdl.SlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles = {};

CompiledEmbryos.mdl.SlideRescaledDorsalAvgAPProfiles.Control = cell(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Control = {};
CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Control.ScaleEstimate = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Control.ScaleSE = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Control.InterceptEstimate = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Control.InterceptSE = NaN(NumMasterProfs, NChannels);

CompiledEmbryos.mdl.SlideRescaledDorsalAvgAPProfiles.Test = cell(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Test = {};
CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Test.ScaleEstimate = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Test.ScaleSE = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Test.InterceptEstimate = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Test.InterceptSE = NaN(NumMasterProfs, NChannels);



for ch_index = [3 5]
ControlSetTF = CompiledEmbryos.IsNC14 & CompiledEmbryos.ControlSetEmbryos & ~isnan(CompiledEmbryos.DubuisEmbryoTimes) & CompiledEmbryos.DubuisEmbryoTimes >= 30;
TestSetTF = CompiledEmbryos.IsNC14 & CompiledEmbryos.TestSetEmbryos & ~isnan(CompiledEmbryos.DubuisEmbryoTimes) & CompiledEmbryos.DubuisEmbryoTimes >= 30;


for ProfIndex = 1:NumMasterProfs
    MasterProfile = MeanSmoothedProfiles{ProfIndex}(:,:,ch_index);
    xsample = CompiledEmbryos.DubuisEmbryoTimes(ControlSetTF);
    
    ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
    
    yset = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(ControlSetTF,:,ch_index);
    
    ysample = ysample(:);
    yset = yset(:);
    
    TFys = ~isnan(ysample) & ~isnan(yset);
    ysample = ysample(TFys);
    yset = yset(TFys);
    beta0 =[1, min(yset,[], 'all')];
    if ~isempty(yset)
        dlm = fitnlm(ysample,yset,f,beta0);
        CompiledEmbryos.mdl.SlideRescaledDorsalAvgAPProfiles.Control{ProfIndex, ch_index} = dlm;
        CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Control.ScaleEstimate(ProfIndex, ch_index) = 1/dlm.Coefficients.Estimate(1);
        CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Control.InterceptEstimate(ProfIndex, ch_index) = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
        CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Control.ScaleSE(ProfIndex, ch_index) = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
        CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Control.InterceptSE(ProfIndex, ch_index) = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
            dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
    end
    
    if (AllSetInfo.Temperatures(exp_index) == 25)
        xsample = CompiledEmbryos.DubuisEmbryoTimes(TestSetTF);
        
        ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
        
        yset = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(TestSetTF,:,ch_index);
        
        beta0 =[1, min(yset,[], 'all')];
        if ~isnan(beta0(2))
            dlm = fitnlm(ysample(:),yset(:),f,beta0);
            CompiledEmbryos.mdl.SlideRescaledDorsalAvgAPProfiles.Test{ProfIndex, ch_index} = dlm;
            CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Test.ScaleEstimate(ProfIndex, ch_index) = 1/dlm.Coefficients.Estimate(1);
            CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Test.InterceptEstimate(ProfIndex, ch_index) = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
            CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Test.ScaleSE(ProfIndex, ch_index) = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
            CompiledEmbryos.ScaleFits.SlideRescaledDorsalAvgAPProfiles.Test.InterceptSE(ProfIndex, ch_index) = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
                dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
        end
    end
end

end

%% Fit to Zero Corrected Slide Rescaled Dorsal Avg AP Profiles

CompiledEmbryos.mdl.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles = {};

CompiledEmbryos.mdl.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Control = cell(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Control = {};
CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Control.ScaleEstimate = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Control.ScaleSE = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Control.InterceptEstimate = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Control.InterceptSE = NaN(NumMasterProfs, NChannels);

CompiledEmbryos.mdl.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Test = cell(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Test = {};
CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Test.ScaleEstimate = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Test.ScaleSE = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Test.InterceptEstimate = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Test.InterceptSE = NaN(NumMasterProfs, NChannels);



for ch_index = [3 5]
ControlSetTF = CompiledEmbryos.IsNC14 & CompiledEmbryos.ControlSetEmbryos & ~isnan(CompiledEmbryos.DubuisEmbryoTimes) & CompiledEmbryos.DubuisEmbryoTimes >= 30;
TestSetTF = CompiledEmbryos.IsNC14 & CompiledEmbryos.TestSetEmbryos & ~isnan(CompiledEmbryos.DubuisEmbryoTimes)& CompiledEmbryos.DubuisEmbryoTimes >= 30;


for ProfIndex = 1:NumMasterProfs
    MasterProfile = MeanSmoothedProfiles{ProfIndex}(:,:,ch_index);
    xsample = CompiledEmbryos.DubuisEmbryoTimes(ControlSetTF);
    
    ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
    
    yset = CompiledEmbryos.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(ControlSetTF,:,ch_index);
    ysample = ysample(:);
    yset = yset(:);
    
    TFys = ~isnan(ysample) & ~isnan(yset);
    ysample = ysample(TFys);
    yset = yset(TFys);
    beta0 =[1, min(yset,[], 'all')];
    if ~isempty(yset)
        dlm = fitnlm(ysample,yset,f,beta0);
        CompiledEmbryos.mdl.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Control{ProfIndex, ch_index} = dlm;
        CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Control.ScaleEstimate(ProfIndex, ch_index) = 1/dlm.Coefficients.Estimate(1);
        CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Control.InterceptEstimate(ProfIndex, ch_index) = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
        CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Control.ScaleSE(ProfIndex, ch_index) = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
        CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Control.InterceptSE(ProfIndex, ch_index) = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
            dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
    end
    
    if (AllSetInfo.Temperatures(exp_index) == 25)
        xsample = CompiledEmbryos.DubuisEmbryoTimes(TestSetTF);
        
        ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
        
        yset = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(TestSetTF,:,ch_index);
        
        beta0 =[1, min(yset,[], 'all')];
        if ~isnan(beta0(2))
            dlm = fitnlm(ysample(:),yset(:),f,beta0);
            CompiledEmbryos.mdl.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Test{ProfIndex, ch_index} = dlm;
            CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Test.ScaleEstimate(ProfIndex, ch_index) = 1/dlm.Coefficients.Estimate(1);
            CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Test.InterceptEstimate(ProfIndex, ch_index) = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
            CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Test.ScaleSE(ProfIndex, ch_index) = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
            CompiledEmbryos.ScaleFits.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Test.InterceptSE(ProfIndex, ch_index) = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
                dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
        end
    end
end

end

%% Fit to Exponential-Fit Bicoid Profiles

CompiledEmbryos.mdl.FitSlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles = {};

CompiledEmbryos.mdl.FitSlideRescaledDorsalAvgAPProfiles.Control = cell(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Control = {};
CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Control.ScaleEstimate = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Control.ScaleSE = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Control.InterceptEstimate = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Control.InterceptSE = NaN(NumMasterProfs, NChannels);

CompiledEmbryos.mdl.FitSlideRescaledDorsalAvgAPProfiles.Test = cell(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Test = {};
CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Test.ScaleEstimate = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Test.ScaleSE = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Test.InterceptEstimate = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Test.InterceptSE = NaN(NumMasterProfs, NChannels);



for ch_index = [3]
ControlSetTF = CompiledEmbryos.IsNC14 & CompiledEmbryos.ControlSetEmbryos & ~isnan(CompiledEmbryos.DubuisEmbryoTimes)& CompiledEmbryos.DubuisEmbryoTimes >= 30;
TestSetTF = CompiledEmbryos.IsNC14 & CompiledEmbryos.TestSetEmbryos & ~isnan(CompiledEmbryos.DubuisEmbryoTimes)& CompiledEmbryos.DubuisEmbryoTimes >= 30;


for ProfIndex = 1:NumMasterProfs
    MasterProfile = MeanSmoothedProfiles{ProfIndex}(:,:,ch_index);
    xsample = CompiledEmbryos.DubuisEmbryoTimes(ControlSetTF);
    
    ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
    
    yset = CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles(ControlSetTF,:,ch_index);
    ysample = ysample(:);
    yset = yset(:);
    
    TFys = ~isnan(ysample) & ~isnan(yset);
    ysample = ysample(TFys);
    yset = yset(TFys);
    beta0 =[1, min(yset,[], 'all')];
    if ~isempty(yset)
        dlm = fitnlm(ysample,yset,f,beta0);
        CompiledEmbryos.mdl.FitSlideRescaledDorsalAvgAPProfiles.Control{ProfIndex, ch_index} = dlm;
        CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Control.ScaleEstimate(ProfIndex, ch_index) = 1/dlm.Coefficients.Estimate(1);
        CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Control.InterceptEstimate(ProfIndex, ch_index) = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
        CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Control.ScaleSE(ProfIndex, ch_index) = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
        CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Control.InterceptSE(ProfIndex, ch_index) = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
            dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
    end
    
    if (AllSetInfo.Temperatures(exp_index) == 25)
        xsample = CompiledEmbryos.DubuisEmbryoTimes(TestSetTF);
        
        ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
        
        yset = CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles(TestSetTF,:,ch_index);
        
        beta0 =[1, min(yset,[], 'all')];
        if ~isnan(beta0(2))
            dlm = fitnlm(ysample(:),yset(:),f,beta0);
            CompiledEmbryos.mdl.FitSlideRescaledDorsalAvgAPProfiles.Test{ProfIndex, ch_index} = dlm;
            CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Test.ScaleEstimate(ProfIndex, ch_index) = 1/dlm.Coefficients.Estimate(1);
            CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Test.InterceptEstimate(ProfIndex, ch_index) = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
            CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Test.ScaleSE(ProfIndex, ch_index) = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
            CompiledEmbryos.ScaleFits.FitSlideRescaledDorsalAvgAPProfiles.Test.InterceptSE(ProfIndex, ch_index) = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
                dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
        end
    end
end

end

%% Fit to Exponential-Fit Zero-Corrected Bicoid Profiles

CompiledEmbryos.mdl.FitZeroedSlideRescaledDorsalAvgAPProfiles = {};
CompiledEmbryos.ScaleFits.FitZeroedSlideRescaledDorsalAvgAPProfiles = {};

CompiledEmbryos.mdl.FitZeroedSlideRescaledDorsalAvgAPProfiles.Control = cell(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.FitZeroedSlideRescaledDorsalAvgAPProfiles.Control = {};
CompiledEmbryos.ScaleFits.FitZeroedSlideRescaledDorsalAvgAPProfiles.Control.ScaleEstimate = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.FitZeroedSlideRescaledDorsalAvgAPProfiles.Control.ScaleSE = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.FitZeroedSlideRescaledDorsalAvgAPProfiles.Control.InterceptEstimate = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.FitZeroedSlideRescaledDorsalAvgAPProfiles.Control.InterceptSE = NaN(NumMasterProfs, NChannels);

CompiledEmbryos.mdl.FitZeroedSlideRescaledDorsalAvgAPProfiles.Test = cell(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.FitZeroedSlideRescaledDorsalAvgAPProfiles.Test = {};
CompiledEmbryos.ScaleFits.FitZeroedSlideRescaledDorsalAvgAPProfiles.Test.ScaleEstimate = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.FitZeroedSlideRescaledDorsalAvgAPProfiles.Test.ScaleSE = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.FitZeroedSlideRescaledDorsalAvgAPProfiles.Test.InterceptEstimate = NaN(NumMasterProfs, NChannels);
CompiledEmbryos.ScaleFits.FitZeroedSlideRescaledDorsalAvgAPProfiles.Test.InterceptSE = NaN(NumMasterProfs, NChannels);



for ch_index = [3]
ControlSetTF = CompiledEmbryos.IsNC14 & CompiledEmbryos.ControlSetEmbryos & ~isnan(CompiledEmbryos.DubuisEmbryoTimes)& CompiledEmbryos.DubuisEmbryoTimes >= 30;
TestSetTF = CompiledEmbryos.IsNC14 & CompiledEmbryos.TestSetEmbryos & ~isnan(CompiledEmbryos.DubuisEmbryoTimes)& CompiledEmbryos.DubuisEmbryoTimes >= 30;


for ProfIndex = 1:NumMasterProfs
    MasterProfile = MeanSmoothedProfiles{ProfIndex}(:,:,ch_index);
    xsample = CompiledEmbryos.DubuisEmbryoTimes(ControlSetTF);
    
    ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
    
    yset = CompiledEmbryos.FitZeroedSlideRescaledDorsalAvgAPProfiles(ControlSetTF,:,ch_index);
    ysample = ysample(:);
    yset = yset(:);
    
    TFys = ~isnan(ysample) & ~isnan(yset);
    ysample = ysample(TFys);
    yset = yset(TFys);
    beta0 =[1, min(yset,[], 'all')];
    if ~isempty(yset)
        dlm = fitnlm(ysample,yset,f,beta0);
        CompiledEmbryos.mdl.FitZeroedSlideRescaledDorsalAvgAPProfiles.Control{ProfIndex, ch_index} = dlm;
        CompiledEmbryos.ScaleFits.FitZeroedSlideRescaledDorsalAvgAPProfiles.Control.ScaleEstimate(ProfIndex, ch_index) = 1/dlm.Coefficients.Estimate(1);
        CompiledEmbryos.ScaleFits.FitZeroedSlideRescaledDorsalAvgAPProfiles.Control.InterceptEstimate(ProfIndex, ch_index) = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
        CompiledEmbryos.ScaleFits.FitZeroedSlideRescaledDorsalAvgAPProfiles.Control.ScaleSE(ProfIndex, ch_index) = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
        CompiledEmbryos.ScaleFits.FitZeroedSlideRescaledDorsalAvgAPProfiles.Control.InterceptSE(ProfIndex, ch_index) = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
            dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
    end
    
    if (AllSetInfo.Temperatures(exp_index) == 25)
        xsample = CompiledEmbryos.DubuisEmbryoTimes(TestSetTF);
        
        ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
        
        yset = CompiledEmbryos.FitZeroedSlideRescaledDorsalAvgAPProfiles(TestSetTF,:,ch_index);
        
        beta0 =[1, min(yset,[], 'all')];
        if ~isnan(beta0(2))
            dlm = fitnlm(ysample(:),yset(:),f,beta0);
            CompiledEmbryos.mdl.FitZeroedSlideRescaledDorsalAvgAPProfiles.Test{ProfIndex, ch_index} = dlm;
            CompiledEmbryos.ScaleFits.FitZeroedSlideRescaledDorsalAvgAPProfiles.Test.ScaleEstimate(ProfIndex, ch_index) = 1/dlm.Coefficients.Estimate(1);
            CompiledEmbryos.ScaleFits.FitZeroedSlideRescaledDorsalAvgAPProfiles.Test.InterceptEstimate(ProfIndex, ch_index) = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
            CompiledEmbryos.ScaleFits.FitZeroedSlideRescaledDorsalAvgAPProfiles.Test.ScaleSE(ProfIndex, ch_index) = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
            CompiledEmbryos.ScaleFits.FitZeroedSlideRescaledDorsalAvgAPProfiles.Test.InterceptSE(ProfIndex, ch_index) = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
                dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
        end
    end
end

end
%%
% CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
% save(CEoutpath, 'CompiledEmbryos');