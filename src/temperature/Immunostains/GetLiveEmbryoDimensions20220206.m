LivePrefix = '2022-02-05-yw_Unsquished_FullEmbryoImages';

FixedPrefixes = {'2021-08-31-yw_25C_MixedAgingPractice_ThreeStackStaging',...
    '2021-09-01-ywGM3_30m25Ccoll_80-100m25Cinc_Slide1',...
    '2021-09-01-ywGM5_30m25Ccoll_95-110m25Cinc_Slide1',...
    '2021-09-02-ywGM2_30m25Ccoll_110-130m25Cinc_Slide1',...
    '2021-09-03-ywGM4_30m25Ccoll_135-150m25Cinc_Slide1',...
    };


FixedCompiledEmbryos = {};
CompiledEmbryos = DefineAPAxesForAllEmbryos(FixedPrefixes{1});
for j = 1:length(FixedPrefixes)
   FixedCompiledEmbryos{j} =  DefineAPAxesForAllEmbryos(FixedPrefixes{j});
end
    

CompiledLiveEmbryos = DefineAPAxesForAllEmbryos(LivePrefix);

%%
LiveAPLengths = CompiledLiveEmbryos.APLengths(CompiledLiveEmbryos.Approved);
LiveDVLengths = CompiledLiveEmbryos.DVLengths(CompiledLiveEmbryos.Approved);

FixedAPLengths = [];
FixedDVLengths = [];
for j = 1:length(FixedPrefixes)
    FixedAPLengths = [FixedAPLengths FixedCompiledEmbryos{j}.APLengths(FixedCompiledEmbryos{j}.Approved)];
    FixedDVLengths = [FixedDVLengths FixedCompiledEmbryos{j}.DVLengths(FixedCompiledEmbryos{j}.Approved)];
end


