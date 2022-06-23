function CompiledEmbryos = PartitionEmbryosTestControl(CompiledEmbryos, exp_index, KnirpsIndex)
%%
APbins = 0:0.025:1;
AllSetInfo = GetFixedSetPrefixInfo;
NEmbryos = length(CompiledEmbryos.Approved);
SetIsFlipped = AllSetInfo.Flipped(exp_index);
NC13NC14Indices = find(CompiledEmbryos.IsNC13orNC14);
NumEmbryosNC13NC14 = length(NC13NC14Indices);
% Zero level APbins for knirps: 0.775-0.875, 0.375-0.45
% find peak between 0.525 and 0.7
% KnirpsBackground < 0.5 for T25C Flipped & slide_index = 1

AllEmbryos = 1:NEmbryos;
MinimalKnirpsAPbins = APbins >= 0.775;% | (APbins >= 0.375 & APbins <= 0.45);
MaxKnirpsAPbins = (APbins >= 0.525 & APbins <= 0.7);
SlideIDs = unique(CompiledEmbryos.SlideIDs);

KnirpsBackground = NaN(1, NEmbryos);
KnirpsPeakValue = NaN(1, NEmbryos);
KnirpsMeanValue = NaN(1, NEmbryos);
for i =1:NumEmbryosNC13NC14
    idx = NC13NC14Indices(i);
    KnirpsBackground(idx) = mean(CompiledEmbryos.DorsalAvgAPProfiles(idx,MinimalKnirpsAPbins,KnirpsIndex).', 'omitnan');
    KnirpsPeakValue(idx) = max(CompiledEmbryos.DorsalAvgAPProfiles(idx,MaxKnirpsAPbins,KnirpsIndex),[],2, 'omitnan');
    KnirpsMeanValue(idx) = mean(CompiledEmbryos.DorsalAvgAPProfiles(idx,:,KnirpsIndex).', 'omitnan');
end
KnirpsDorsalRatios = KnirpsPeakValue./KnirpsBackground;
CompiledEmbryos.KnirpsBackground = KnirpsBackground;
CompiledEmbryos.KnirpsDorsalRatios = KnirpsDorsalRatios;
CompiledEmbryos.KnirpsPeakValue = KnirpsPeakValue;
CompiledEmbryos.KnirpsMeanValue = KnirpsMeanValue;

NSlides = length(unique(CompiledEmbryos.SlideIDs));
CompiledEmbryos.ControlSetEmbryos = zeros(1, NEmbryos, 'logical');
CompiledEmbryos.TestSetEmbryos = zeros(1, NEmbryos, 'logical');
CompiledEmbryos.NControlEmbryos = NaN(1,NSlides );
CompiledEmbryos.NTestEmbryos = NaN(1, NSlides);
CompiledEmbryos.SlideControlSetEmbryos = zeros(NSlides, NEmbryos, 'logical');

for slide_index=SlideIDs
   IsControlSlideEmbryos = ismember(1:NEmbryos, NC13NC14Indices) & CompiledEmbryos.SlideIDs == slide_index;
   IsTestSlideEmbryos = ismember(1:NEmbryos, NC13NC14Indices) & CompiledEmbryos.SlideIDs == slide_index;
   if ~SetIsFlipped
       
       if AllSetInfo.UsePeakRatioToIdentifyControls{exp_index}(slide_index)
           IsControlSlideEmbryos = IsControlSlideEmbryos & KnirpsBackground > AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff{exp_index}(slide_index);
           IsTestSlideEmbryos = IsTestSlideEmbryos & KnirpsBackground < AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff{exp_index}(slide_index);
       end
       if AllSetInfo.UsePeakValueToIdentifyControls{exp_index}(slide_index)
           IsControlSlideEmbryos = IsControlSlideEmbryos & KnirpsPeakValue > AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff{exp_index}(slide_index);
           IsTestSlideEmbryos = IsTestSlideEmbryos & KnirpsPeakValue < AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff{exp_index}(slide_index);
       end
       
       if AllSetInfo.UseMeanValueToIdentifyControls{exp_index}(slide_index)
           IsControlSlideEmbryos = IsControlSlideEmbryos & KnirpsMeanValue > AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff{exp_index}(slide_index);
           IsTestSlideEmbryos = IsTestSlideEmbryos & KnirpsMeanValue < AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff{exp_index}(slide_index);
       end
   else
        if AllSetInfo.UsePeakRatioToIdentifyControls{exp_index}(slide_index)
           IsControlSlideEmbryos = IsControlSlideEmbryos & KnirpsBackground < AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff{exp_index}(slide_index);
           IsTestSlideEmbryos = IsTestSlideEmbryos & KnirpsBackground > AllSetInfo.DorsalAvgAPKnirpsBackgroundCutoff{exp_index}(slide_index);
       end
       if AllSetInfo.UsePeakValueToIdentifyControls{exp_index}(slide_index)
           IsControlSlideEmbryos = IsControlSlideEmbryos & KnirpsPeakValue < AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff{exp_index}(slide_index);
           IsTestSlideEmbryos = IsTestSlideEmbryos & KnirpsPeakValue > AllSetInfo.DorsalAvgAPKnirpsPeakValueCutoff{exp_index}(slide_index);
       end
       
       if AllSetInfo.UseMeanValueToIdentifyControls{exp_index}(slide_index)
           IsControlSlideEmbryos = IsControlSlideEmbryos & KnirpsMeanValue < AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff{exp_index}(slide_index);
           IsTestSlideEmbryos = IsTestSlideEmbryos & KnirpsMeanValue > AllSetInfo.DorsalAvgAPKnirpsMeanValueCutoff{exp_index}(slide_index);
       end
       
   end
   
   if AllSetInfo.MinimumFixCorrectedDeltaFC{exp_index}(slide_index) > 0
   IsControlSlideEmbryos = IsControlSlideEmbryos & ([CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].' >= AllSetInfo.MinimumFixCorrectedDeltaFC{exp_index}(slide_index));
   IsTestSlideEmbryos = IsTestSlideEmbryos & ([CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].' >= AllSetInfo.MinimumFixCorrectedDeltaFC{exp_index}(slide_index));
   end
   CompiledEmbryos.SlideControlSetEmbryos(slide_index,:) = IsControlSlideEmbryos;
   CompiledEmbryos.SlideTestSetEmbryos(slide_index,:) = IsTestSlideEmbryos;
   CompiledEmbryos.NControlEmbryos(slide_index) = sum(IsControlSlideEmbryos);
   CompiledEmbryos.NTestEmbryos(slide_index) = sum(IsTestSlideEmbryos);
   CompiledEmbryos.ControlSetEmbryos  = CompiledEmbryos.ControlSetEmbryos | IsControlSlideEmbryos;
   CompiledEmbryos.TestSetEmbryos  = CompiledEmbryos.TestSetEmbryos | IsTestSlideEmbryos;
end
