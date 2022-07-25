exp_index = 13;
CompiledEmbryos = AllCompiledEmbryos{13};

UniqueSlideIDs = unique(CompiledEmbryos.SlideIDs);
figure(1)
plot(APbins, CompiledEmbryos.DorsalAvgAPProfiles(CompiledEmbryos.IsNC14 & CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.SlideIDs == 3,:,3).')
