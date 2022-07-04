function CompiledEmbryos = AddNonBicoidZeroCorrection(CompiledEmbryos, exp_index)
%%
AllSetInfo = GetFixedSetPrefixInfo;

AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';

SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
APbins = 0:0.025:1;
NumAPbins = length(APbins);
NChannels = 5;
ChannelZeroBins = ones(NChannels, NumAPbins, 'logical');
ChannelZeroBins(3, :) = APbins >= 0.85 & APbins <= 0.9;
ChannelZeroBins(4, :) = APbins >= 0.85 & APbins <= 0.9;
ChannelZeroBins(5, :) = APbins > 0.6 & APbins < 0.7;

%%
if ~isfield(CompiledEmbryos, 'ZeroCorrectedSlideRescaledDorsalAvgAPProfiles')
CompiledEmbryos.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles = NaN(size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles));
end

Channels  = 4:5;
for ch_index = Channels
    for i = 1:size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles, 1)
        

        CompiledEmbryos.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(i,:,ch_index) = ...
            CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(i,:,ch_index)-...
            mean(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(i,ChannelZeroBins(ch_index, :),ch_index), 'omitnan');
        
        
        
    end
end



%%


CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
save(CEoutpath, 'CompiledEmbryos');
