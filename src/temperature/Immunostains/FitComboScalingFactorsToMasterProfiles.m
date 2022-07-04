function AllCompiledEmbryos = FitComboScalingFactorsToMasterProfiles(AllCompiledEmbryos)
%%
NChannels = 5;
AllSetInfo = GetFixedSetPrefixInfo;

AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';

BootstrappedProfilePath = 'S:/Gabriella/Dropbox/BootstrappedTestData/25CBootstrappedGaussianSmoothedProfiles.mat';
load(BootstrappedProfilePath, 'MeanSmoothedProfiles', 'SmoothedProfileSEs', 'xfits');
NumMasterProfs = length(MeanSmoothedProfiles);


f = @(b,x) b(1).*x + b(2);

NumSets = length(AllSetInfo.Temperatures);


SetStrings = {'Rep1', 'Rep2', 'Flipped'};
%% Fit to Slide Rescal Dorsal Avg AP Profiles
for temp = [25 27.5 22.5 20 17.5]
    Rep1Index = find(AllSetInfo.Temperatures == temp & AllSetInfo.Replicates == 1);
    Rep2Index = find(AllSetInfo.Temperatures == temp & AllSetInfo.Replicates == 2);
    FlippedIndex = find(AllSetInfo.Temperatures == temp & AllSetInfo.Replicates == 0);
    SetIndices = [Rep1Index, Rep2Index, FlippedIndex];
    CCEs = {};
    CCEs{1} = AllCompiledEmbryos{Rep1Index};
    CCEs{2} = AllCompiledEmbryos{Rep2Index};
    CCEs{3} = AllCompiledEmbryos{FlippedIndex};
    SetLabels = {};
    SetLabels{1} = AllSetInfo.SetLabels{Rep1Index};
    SetLabels{2} = AllSetInfo.SetLabels{Rep2Index};
    SetLabels{3} = AllSetInfo.SetLabels{FlippedIndex};
    OutEmbryoPaths = {};
    OutEmbryoPaths{1} = [AllSetsCombinedEmbryosPath, SetLabels{1}];
    OutEmbryoPaths{2} = [AllSetsCombinedEmbryosPath, SetLabels{2}];
    OutEmbryoPaths{3} = [AllSetsCombinedEmbryosPath, SetLabels{3}];

    ControlSetTFs = {};
    TestSetTFs = {};
    FitTypeStrings = {'Rep1Rep2ComboFitSlideRescaledDorsalAvgAPProfiles', 'Rep1FlippedComboFitSlideRescaledDorsalAvgAPProfiles',...
        'Rep2FlippedComboFitSlideRescaledDorsalAvgAPProfiles', 'AllCombinedFitSlideRescaledDorsalAvgAPProfiles'};
    FitTypeIncludedSets = zeros(4, 3, 'logical');
    FitTypeIncludedSets(1, 1:2) = 1;
    FitTypeIncludedSets(2, [1 3]) = 1;
    FitTypeIncludedSets(3, 2:3) = 1;
    FitTypeIncludedSets(4, 1:3) = 1;
    for i = 1:3
        if ~isfield('CompiledEmbryos', 'mdl')
            CCEs{i}.mdl = {};
        end
        if ~isfield('CompiledEmbryos', 'ScaleFits')
            CCEs{i}.ScaleFits  = {};
        end
        ControlSetTFs{i} = CCEs{i}.IsNC14 & CCEs{i}.ControlSetEmbryos & ~isnan(CCEs{i}.DubuisEmbryoTimes);
        TestSetTFs{i} = CCEs{i}.IsNC14 & CCEs{i}.TestSetEmbryos & ~isnan(CCEs{i}.DubuisEmbryoTimes);
        
        for j = 1:length(FitTypeStrings)
            CCEs{i}.mdl.(FitTypeStrings{j}) = {};
            CCEs{i}.ScaleFits.(FitTypeStrings{j}) = {};
            for k = 1:length(SetStrings)
                CCEs{i}.mdl.(FitTypeStrings{j}).(SetStrings{k}) = cell(NumMasterProfs, NChannels);
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}) = {};
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}).ScaleEstimate = NaN(NumMasterProfs, NChannels);
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}).ScaleSE = NaN(NumMasterProfs, NChannels);
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}).InterceptEstimate = NaN(NumMasterProfs, NChannels);
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}).InterceptSE = NaN(NumMasterProfs, NChannels);
            end
            
        end

    end

%% Fit to Exponential-Fit Bicoid Profiles
    for ch_index = 3
        for fittype_index = 1:size(FitTypeIncludedSets, 1)
            for set_index = 1:length(SetStrings)
                ControlSetTF = [];
                AllDubuisTimes = [];
                Allyset = [];
                for i = 1:3
                    if FitTypeIncludedSets(fittype_index, i)
                        ControlSetTF = [ControlSetTF ...
                            CCEs{i}.IsNC14 & CCEs{i}.ControlSetEmbryos & ~isnan(CCEs{i}.DubuisEmbryoTimes)];
                    else
                        ControlSetTF = [ControlSetTF  zeros(1, length(CCEs{i}.IsNC14), 'logical')];
                    end
                    if set_index ~= i
                        Allyset = [Allyset ; CCEs{i}.UnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.(SetStrings{set_index})(:,:,ch_index)];
                    else
                        Allyset = [Allyset ; CCEs{i}.FitSlideRescaledDorsalAvgAPProfiles(:,:,ch_index)];
                        CCEs{i}.UnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.(SetStrings{set_index})(:,:,ch_index) = ...
                            CCEs{i}.FitSlideRescaledDorsalAvgAPProfiles(:,:,ch_index);
                    end
                    AllDubuisTimes = [AllDubuisTimes CCEs{i}.DubuisEmbryoTimes];
                end
                
                ControlSetTF = logical(ControlSetTF);
                xsample = AllDubuisTimes(ControlSetTF);
                yset = Allyset(ControlSetTF,:);
                yset_flat = yset(:);
                beta0 =[1, min(yset,[], 'all')];
                if ~isnan(beta0(2))
                    for ProfIndex = 1:NumMasterProfs
                        MasterProfile = MeanSmoothedProfiles{ProfIndex}(:,:,ch_index);
                        ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
                        ysample_flat = ysample(:);
                       
                        TFys = ~isnan(ysample_flat) & ~isnan(yset_flat);
                        dlm =  fitnlm(ysample_flat(TFys),yset_flat(TFys),f,beta0);
                        for i = 1:3
                            if FitTypeIncludedSets(fittype_index, i)
                                
                                CCEs{i}.mdl.(FitTypeStrings{fittype_index}).(SetStrings{set_index}){ProfIndex, ch_index} = dlm;
                                CCEs{i}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).ScaleEstimate(ProfIndex, ch_index) = 1/dlm.Coefficients.Estimate(1);
                                CCEs{i}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).InterceptEstimate(ProfIndex, ch_index)  = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
                                CCEs{i}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).ScaleSE(ProfIndex, ch_index)  = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
                                CCEs{i}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).InterceptSE(ProfIndex, ch_index)  = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
                                    dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
                            end
                        end
                    end
                end
                
                if sum( FitTypeIncludedSets(fittype_index, :)) < 3
                    NotSetIndex = find(~FitTypeIncludedSets(fittype_index, :));
                    ControlSetTF = CCEs{NotSetIndex}.IsNC14 & ...
                        CCEs{NotSetIndex}.ControlSetEmbryos & ~isnan(CCEs{NotSetIndex}.DubuisEmbryoTimes) ;
                    xsample = CCEs{NotSetIndex}.DubuisEmbryoTimes(ControlSetTF);
                    if set_index ~= NotSetIndex
                        yset =  CCEs{NotSetIndex}.UnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.(SetStrings{set_index})(ControlSetTF,:,ch_index);
                    else
                        yset = CCEs{NotSetIndex}.FitSlideRescaledDorsalAvgAPProfiles(ControlSetTF,:,ch_index);
                    end
                    
                    beta0 =[1, min(yset,[], 'all')];
                    if ~isnan(beta0(2))
                        for ProfIndex = 1:NumMasterProfs
                            MasterProfile = MeanSmoothedProfiles{ProfIndex}(:,:,ch_index);
                            ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
                            
                            dlm =  fitnlm(ysample(:),yset(:),f,beta0);
                      
                            
                            CCEs{NotSetIndex}.mdl.(FitTypeStrings{fittype_index}).(SetStrings{set_index}){ProfIndex, ch_index} = dlm;
                            CCEs{NotSetIndex}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).ScaleEstimate(ProfIndex, ch_index) = 1/dlm.Coefficients.Estimate(1);
                            CCEs{NotSetIndex}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).InterceptEstimate(ProfIndex, ch_index)  = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
                            CCEs{NotSetIndex}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).ScaleSE(ProfIndex, ch_index)  = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
                            CCEs{NotSetIndex}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).InterceptSE(ProfIndex, ch_index)  = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
                                dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
                            
                        end
                    end
                    
                end
                %
                
            end
        end
    end
    
%% Add Zero Corrected Fits
    FitTypeStrings = {'Rep1Rep2ComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles', 'Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles',...
        'Rep2FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles', 'AllCombinedZeroCorrectedSlideRescaledDorsalAvgAPProfiles'};
    for i = 1:3
        for j = 1:length(FitTypeStrings)
            CCEs{i}.mdl.(FitTypeStrings{j}) = {};
            CCEs{i}.ScaleFits.(FitTypeStrings{j}) = {};
            for k = 1:length(SetStrings)
                CCEs{i}.mdl.(FitTypeStrings{j}).(SetStrings{k}) = cell(NumMasterProfs, NChannels);
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}) = {};
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}).ScaleEstimate = NaN(NumMasterProfs, NChannels);
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}).ScaleSE = NaN(NumMasterProfs, NChannels);
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}).InterceptEstimate = NaN(NumMasterProfs, NChannels);
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}).InterceptSE = NaN(NumMasterProfs, NChannels);
            end
            
        end
    end

    for ch_index = [3 5]
        for fittype_index = 1:size(FitTypeIncludedSets, 1)
            for set_index = 1:length(SetStrings)
                ControlSetTF = [];
                AllDubuisTimes = [];
                Allyset = [];
                for i = 1:3
                    if FitTypeIncludedSets(fittype_index, i)
                        ControlSetTF = [ControlSetTF ...
                            CCEs{i}.IsNC14 & CCEs{i}.ControlSetEmbryos & ~isnan(CCEs{i}.DubuisEmbryoTimes)];
                    else
                        ControlSetTF = [ControlSetTF  zeros(1, length(CCEs{i}.IsNC14), 'logical')];
                    end
                    if set_index ~= i
                        Allyset = [Allyset ; CCEs{i}.UnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.(SetStrings{set_index})(:,:,ch_index)];
                    else
                        Allyset = [Allyset ; CCEs{i}.FitSlideRescaledDorsalAvgAPProfiles(:,:,ch_index)];
                        CCEs{i}.UnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.(SetStrings{set_index})(:,:,ch_index) = ...
                            CCEs{i}.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(:,:,ch_index);
                    end
                    AllDubuisTimes = [AllDubuisTimes CCEs{i}.DubuisEmbryoTimes];
                end
                
                ControlSetTF = logical(ControlSetTF);
                xsample = AllDubuisTimes(ControlSetTF);
                yset = Allyset(ControlSetTF,:);
                yset_flat = yset(:);
                beta0 =[1, min(yset,[], 'all')];
                if ~isnan(beta0(2))
                    for ProfIndex = 1:NumMasterProfs
                        MasterProfile = MeanSmoothedProfiles{ProfIndex}(:,:,ch_index);
                        ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
                        ysample_flat = ysample(:);
                       
                        TFys = ~isnan(ysample_flat) & ~isnan(yset_flat);
                        dlm =  fitnlm(ysample_flat(TFys),yset_flat(TFys),f,beta0);
                        for i = 1:3
                            if FitTypeIncludedSets(fittype_index, i)
                                
                                CCEs{i}.mdl.(FitTypeStrings{fittype_index}).(SetStrings{set_index}){ProfIndex, ch_index} = dlm;
                                CCEs{i}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).ScaleEstimate(ProfIndex, ch_index) = 1/dlm.Coefficients.Estimate(1);
                                CCEs{i}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).InterceptEstimate(ProfIndex, ch_index)  = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
                                CCEs{i}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).ScaleSE(ProfIndex, ch_index)  = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
                                CCEs{i}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).InterceptSE(ProfIndex, ch_index)  = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
                                    dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
                            end
                        end
                    end
                end
                
                if sum( FitTypeIncludedSets(fittype_index, :)) < 3
                    NotSetIndex = find(~FitTypeIncludedSets(fittype_index, :));
                    ControlSetTF = CCEs{NotSetIndex}.IsNC14 & ...
                        CCEs{NotSetIndex}.ControlSetEmbryos & ~isnan(CCEs{NotSetIndex}.DubuisEmbryoTimes) ;
                    xsample = CCEs{NotSetIndex}.DubuisEmbryoTimes(ControlSetTF);
                    if set_index ~= NotSetIndex
                        yset =  CCEs{NotSetIndex}.UnivScaledProfiles.TestSetZeroCorrectedSlideRescaledDorsalAvgAPProfiles.(SetStrings{set_index})(ControlSetTF,:,ch_index);
                    else
                        yset = CCEs{NotSetIndex}.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(ControlSetTF,:,ch_index);
                    end
                    
                    beta0 =[1, min(yset,[], 'all')];
                    if ~isnan(beta0(2))
                        for ProfIndex = 1:NumMasterProfs
                            MasterProfile = MeanSmoothedProfiles{ProfIndex}(:,:,ch_index);
                            ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
                            
                            dlm =  fitnlm(ysample(:),yset(:),f,beta0);
                      
                            
                            CCEs{NotSetIndex}.mdl.(FitTypeStrings{fittype_index}).(SetStrings{set_index}){ProfIndex, ch_index} = dlm;
                            CCEs{NotSetIndex}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).ScaleEstimate(ProfIndex, ch_index) = 1/dlm.Coefficients.Estimate(1);
                            CCEs{NotSetIndex}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).InterceptEstimate(ProfIndex, ch_index)  = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
                            CCEs{NotSetIndex}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).ScaleSE(ProfIndex, ch_index)  = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
                            CCEs{NotSetIndex}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).InterceptSE(ProfIndex, ch_index)  = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
                                dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
                            
                        end
                    end
                    
                end
                %
                
            end
        end
    end
    
    %% Add Slide Rescaled Avg Profile 
    FitTypeStrings = {'Rep1Rep2ComboSlideRescaledDorsalAvgAPProfiles', 'Rep1FlippedComboSlideRescaledDorsalAvgAPProfiles',...
        'Rep2FlippedComboSlideRescaledDorsalAvgAPProfiles', 'AllCombinedSlideRescaledDorsalAvgAPProfiles'};
    for i = 1:3
        for j = 1:length(FitTypeStrings)
            CCEs{i}.mdl.(FitTypeStrings{j}) = {};
            CCEs{i}.ScaleFits.(FitTypeStrings{j}) = {};
            for k = 1:length(SetStrings)
                CCEs{i}.mdl.(FitTypeStrings{j}).(SetStrings{k}) = cell(NumMasterProfs, NChannels);
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}) = {};
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}).ScaleEstimate = NaN(NumMasterProfs, NChannels);
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}).ScaleSE = NaN(NumMasterProfs, NChannels);
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}).InterceptEstimate = NaN(NumMasterProfs, NChannels);
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}).InterceptSE = NaN(NumMasterProfs, NChannels);
            end 
        end
    end

    for ch_index = [3 5]
        for fittype_index = 1:size(FitTypeIncludedSets, 1)
            for set_index = 1:length(SetStrings)
                ControlSetTF = [];
                AllDubuisTimes = [];
                Allyset = [];
                for i = 1:3
                    if FitTypeIncludedSets(fittype_index, i)
                        ControlSetTF = [ControlSetTF ...
                            CCEs{i}.IsNC14 & CCEs{i}.ControlSetEmbryos & ~isnan(CCEs{i}.DubuisEmbryoTimes)];
                    else
                        ControlSetTF = [ControlSetTF  zeros(1, length(CCEs{i}.IsNC14), 'logical')];
                    end
                    if set_index ~= i
                        Allyset = [Allyset ; CCEs{i}.UnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.(SetStrings{set_index})(:,:,ch_index)];
                    else
                        Allyset = [Allyset ; CCEs{i}.SlideRescaledDorsalAvgAPProfiles(:,:,ch_index)];
                        CCEs{i}.UnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.(SetStrings{set_index})(:,:,ch_index) = ...
                            CCEs{i}.SlideRescaledDorsalAvgAPProfiles(:,:,ch_index);
                    end
                    AllDubuisTimes = [AllDubuisTimes CCEs{i}.DubuisEmbryoTimes];
                end
                
                ControlSetTF = logical(ControlSetTF);
                xsample = AllDubuisTimes(ControlSetTF);
                yset = Allyset(ControlSetTF,:);
                yset_flat = yset(:);
                beta0 =[1, min(yset,[], 'all')];
                if ~isnan(beta0(2))
                    for ProfIndex = 1:NumMasterProfs
                        MasterProfile = MeanSmoothedProfiles{ProfIndex}(:,:,ch_index);
                        ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
                        ysample_flat = ysample(:);
                       
                        TFys = ~isnan(ysample_flat) & ~isnan(yset_flat);
                        dlm =  fitnlm(ysample_flat(TFys),yset_flat(TFys),f,beta0);
                        for i = 1:3
                            if FitTypeIncludedSets(fittype_index, i)
                                
                                CCEs{i}.mdl.(FitTypeStrings{fittype_index}).(SetStrings{set_index}){ProfIndex, ch_index} = dlm;
                                CCEs{i}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).ScaleEstimate(ProfIndex, ch_index) = 1/dlm.Coefficients.Estimate(1);
                                CCEs{i}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).InterceptEstimate(ProfIndex, ch_index)  = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
                                CCEs{i}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).ScaleSE(ProfIndex, ch_index)  = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
                                CCEs{i}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).InterceptSE(ProfIndex, ch_index)  = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
                                    dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
                            end
                        end
                    end
                end
                
                if sum( FitTypeIncludedSets(fittype_index, :)) < 3
                    NotSetIndex = find(~FitTypeIncludedSets(fittype_index, :));
                    ControlSetTF = CCEs{NotSetIndex}.IsNC14 & ...
                        CCEs{NotSetIndex}.ControlSetEmbryos & ~isnan(CCEs{NotSetIndex}.DubuisEmbryoTimes) ;
                    xsample = CCEs{NotSetIndex}.DubuisEmbryoTimes(ControlSetTF);
                    if set_index ~= NotSetIndex
                        yset =  CCEs{NotSetIndex}.UnivScaledProfiles.TestSetSlideRescaledDorsalAvgAPProfiles.(SetStrings{set_index})(ControlSetTF,:,ch_index);
                    else
                        yset = CCEs{NotSetIndex}.SlideRescaledDorsalAvgAPProfiles(ControlSetTF,:,ch_index);
                    end
                    
                    beta0 =[1, min(yset,[], 'all')];
                    if ~isnan(beta0(2))
                        for ProfIndex = 1:NumMasterProfs
                            MasterProfile = MeanSmoothedProfiles{ProfIndex}(:,:,ch_index);
                            ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
                            
                            dlm =  fitnlm(ysample(:),yset(:),f,beta0);
                      
                            
                            CCEs{NotSetIndex}.mdl.(FitTypeStrings{fittype_index}).(SetStrings{set_index}){ProfIndex, ch_index} = dlm;
                            CCEs{NotSetIndex}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).ScaleEstimate(ProfIndex, ch_index) = 1/dlm.Coefficients.Estimate(1);
                            CCEs{NotSetIndex}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).InterceptEstimate(ProfIndex, ch_index)  = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
                            CCEs{NotSetIndex}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).ScaleSE(ProfIndex, ch_index)  = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
                            CCEs{NotSetIndex}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).InterceptSE(ProfIndex, ch_index)  = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
                                dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
                            
                        end
                    end
                    
                end
                %
                
            end
        end
    end
    
        %% Add Hybrid 
    FitTypeStrings = {'Rep1Rep2ComboHybridRescaledDorsalAvgAPProfiles', 'Rep1FlippedComboHybridRescaledDorsalAvgAPProfiles',...
        'Rep2FlippedComboHybridRescaledDorsalAvgAPProfiles', 'AllCombinedHybridRescaledDorsalAvgAPProfiles'};
    for i = 1:3
        for j = 1:length(FitTypeStrings)
            CCEs{i}.mdl.(FitTypeStrings{j}) = {};
            CCEs{i}.ScaleFits.(FitTypeStrings{j}) = {};
            for k = 1:length(SetStrings)
                CCEs{i}.mdl.(FitTypeStrings{j}).(SetStrings{k}) = cell(NumMasterProfs, NChannels);
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}) = {};
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}).ScaleEstimate = NaN(NumMasterProfs, NChannels);
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}).ScaleSE = NaN(NumMasterProfs, NChannels);
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}).InterceptEstimate = NaN(NumMasterProfs, NChannels);
                CCEs{i}.ScaleFits.(FitTypeStrings{j}).(SetStrings{k}).InterceptSE = NaN(NumMasterProfs, NChannels);
            end 
        end
    end

    for ch_index = 3
        for fittype_index = 1:size(FitTypeIncludedSets, 1)
            for set_index = 1:length(SetStrings)
                ControlSetTF = [];
                AllDubuisTimes = [];
                Allyset = [];
                for i = 1:3
                    if FitTypeIncludedSets(fittype_index, i)
                        ControlSetTF = [ControlSetTF ...
                            CCEs{i}.IsNC14 & CCEs{i}.ControlSetEmbryos & ~isnan(CCEs{i}.DubuisEmbryoTimes)];
                    else
                        ControlSetTF = [ControlSetTF  zeros(1, length(CCEs{i}.IsNC14), 'logical')];
                    end
                    if set_index ~= i
                        Allyset = [Allyset ; CCEs{i}.UnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.(SetStrings{set_index})(:,:,ch_index)];
                    else
                        Allyset = [Allyset ; CCEs{i}.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(:,:,ch_index)];
                        CCEs{i}.UnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.(SetStrings{set_index})(:,:,ch_index) = ...
                            CCEs{i}.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(:,:,ch_index);
                    end
                    AllDubuisTimes = [AllDubuisTimes CCEs{i}.DubuisEmbryoTimes];
                end
                
                ControlSetTF = logical(ControlSetTF);
                xsample = AllDubuisTimes(ControlSetTF);
                yset = Allyset(ControlSetTF,:);
                yset_flat = yset(:);
                beta0 =[1, min(yset,[], 'all')];
                if ~isnan(beta0(2))
                    for ProfIndex = 1:NumMasterProfs
                        MasterProfile = MeanSmoothedProfiles{ProfIndex}(:,:,ch_index);
                        ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
                        ysample_flat = ysample(:);
                       
                        TFys = ~isnan(ysample_flat) & ~isnan(yset_flat);
                        dlm =  fitnlm(ysample_flat(TFys),yset_flat(TFys),f,beta0);
                        for i = 1:3
                            if FitTypeIncludedSets(fittype_index, i)
                                
                                CCEs{i}.mdl.(FitTypeStrings{fittype_index}).(SetStrings{set_index}){ProfIndex, ch_index} = dlm;
                                CCEs{i}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).ScaleEstimate(ProfIndex, ch_index) = 1/dlm.Coefficients.Estimate(1);
                                CCEs{i}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).InterceptEstimate(ProfIndex, ch_index)  = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
                                CCEs{i}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).ScaleSE(ProfIndex, ch_index)  = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
                                CCEs{i}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).InterceptSE(ProfIndex, ch_index)  = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
                                    dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
                            end
                        end
                    end
                end
                
                if sum( FitTypeIncludedSets(fittype_index, :)) < 3
                    NotSetIndex = find(~FitTypeIncludedSets(fittype_index, :));
                    ControlSetTF = CCEs{NotSetIndex}.IsNC14 & ...
                        CCEs{NotSetIndex}.ControlSetEmbryos & ~isnan(CCEs{NotSetIndex}.DubuisEmbryoTimes) ;
                    xsample = CCEs{NotSetIndex}.DubuisEmbryoTimes(ControlSetTF);
                    if set_index ~= NotSetIndex
                        yset =  CCEs{NotSetIndex}.UnivScaledProfiles.HybridTestSetRescaledDorsalAvgAPProfiles.(SetStrings{set_index})(ControlSetTF,:,ch_index);
                    else
                        yset = CCEs{NotSetIndex}.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(ControlSetTF,:,ch_index);
                    end
                    
                    beta0 =[1, min(yset,[], 'all')];
                    if ~isnan(beta0(2))
                        for ProfIndex = 1:NumMasterProfs
                            MasterProfile = MeanSmoothedProfiles{ProfIndex}(:,:,ch_index);
                            ysample = GetMasterProfileForEmbryoTimes(xsample, MasterProfile, xfits);
                            
                            dlm =  fitnlm(ysample(:),yset(:),f,beta0);
                      
                            
                            CCEs{NotSetIndex}.mdl.(FitTypeStrings{fittype_index}).(SetStrings{set_index}){ProfIndex, ch_index} = dlm;
                            CCEs{NotSetIndex}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).ScaleEstimate(ProfIndex, ch_index) = 1/dlm.Coefficients.Estimate(1);
                            CCEs{NotSetIndex}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).InterceptEstimate(ProfIndex, ch_index)  = -dlm.Coefficients.Estimate(2)/dlm.Coefficients.Estimate(1);
                            CCEs{NotSetIndex}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).ScaleSE(ProfIndex, ch_index)  = dlm.Coefficients.SE(1)/(dlm.Coefficients.Estimate(1)^2);
                            CCEs{NotSetIndex}.ScaleFits.(FitTypeStrings{fittype_index}).(SetStrings{set_index}).InterceptSE(ProfIndex, ch_index)  = sqrt(dlm.Coefficients.SE(1)^2*dlm.Coefficients.Estimate(2)^2/(dlm.Coefficients.Estimate(1)^4)+...
                                dlm.Coefficients.SE(2)^2/(dlm.Coefficients.Estimate(1)^2));
                            
                        end
                    end
                    
                end
                %
                
            end
        end
    end
    
    for i = 1:3
        CompiledEmbryos = CCEs{i};
        exp_index = SetIndices(i);
        AllCompiledEmbryos{exp_index} = CCEs{i};
        SetLabel = AllSetInfo.SetLabels{exp_index};
        OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];

        CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
        save(CEoutpath, 'CompiledEmbryos');
        clear CompiledEmbryos  exp_index
    end

    disp(['Temperature: ', num2str(temp)])
    
end


%%

end
%%
