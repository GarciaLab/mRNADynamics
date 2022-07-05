function CompiledEmbryos = AddFitBicoidProfiles(CompiledEmbryos, exp_index)
%%
AllSetInfo = GetFixedSetPrefixInfo;

AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';

SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];

%%
f = @(b,x) b(1).*exp(b(2).*x) + b(3);
fitx2 = (0.0125:0.025:(1+0.0125));
UseFitX2 = fitx2 >= 0.2 & fitx2 < 0.925;
UsePredX2 = fitx2 >= 0.1 & fitx2 < 0.925;
CompiledEmbryos.RescaledAvgFitModels = cell(1, size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles, 1));
CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles = NaN(size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles));

CompiledEmbryos.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles = NaN(size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles));

NumAPbins = size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles, 2);
NChannels = size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles, 3);
ch_index = 3;

for i = 1:size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles, 1)
    beta0 =[ max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(i,:,ch_index)), -3, min(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(i,:,ch_index))];
    if ~isnan(beta0(1))
    CompiledEmbryos.RescaledAvgFitModels{i} = fitnlm(fitx2(UseFitX2),CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(i,UseFitX2,ch_index),f,beta0);

    CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles(i,UsePredX2,ch_index) = f(CompiledEmbryos.RescaledAvgFitModels{i}.Coefficients.Estimate, fitx2(UsePredX2))-...
        CompiledEmbryos.RescaledAvgFitModels{i}.Coefficients.Estimate(3);
    CompiledEmbryos.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(i,:,ch_index) = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(i,:,ch_index)-...
        CompiledEmbryos.RescaledAvgFitModels{i}.Coefficients.Estimate(3);
    
        if CompiledEmbryos.RescaledAvgFitModels{i}.Rsquared.Ordinary < 0.9
            CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles(i,UsePredX2,ch_index) = NaN(1, sum(UsePredX2), 1);
            CompiledEmbryos.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(i,UsePredX2,ch_index) = NaN(1, sum(UsePredX2), 1);
            CompiledEmbryos.Approved(i) = false;
            CompiledEmbryos.IsNC14(i) = false;
            CompiledEmbryos.IsNC13(i) = false;
            CompiledEmbryos.IsNC13orNC14(i) = false;
            CompiledEmbryos.NumEmbryosNC13 = sum(CompiledEmbryos.IsNC13 > 0 );
            CompiledEmbryos.NumEmbryosNC14 = sum(CompiledEmbryos.IsNC14 > 0 );
            CompiledEmbryos.NumEmbryosNC13NC14 = sum(CompiledEmbryos.IsNC13orNC14 > 0 );
            CompiledEmbryos.ControlSetEmbryos(i) = false;
            CompiledEmbryos.TestSetEmbryos(i) = false;
            SlideIDs = unique(CompiledEmbryos.SlideIDs);
            for slide_index = SlideIDs
                CompiledEmbryos.NControlEmbryos(slide_index) = sum(CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.SlideIDs == slide_index);
                CompiledEmbryos.NTestEmbryos(slide_index) = sum(CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.SlideIDs == slide_index);
                CompiledEmbryos.SlideControlSetEmbryos(slide_index,:) = CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.SlideIDs == slide_index;
                CompiledEmbryos.SlideTestSetEmbryos(slide_index,:) = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.SlideIDs == slide_index;
            end
        end
    end
end



%%


CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
save(CEoutpath, 'CompiledEmbryos');
