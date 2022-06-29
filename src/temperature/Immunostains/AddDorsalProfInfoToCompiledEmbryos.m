function CompiledEmbryos = AddDorsalProfInfoToCompiledEmbryos(Prefix)


liveExperiment = LiveExperiment(Prefix);

CompiledEmbryoPath = [liveExperiment.resultsFolder, filesep, 'CompiledEmbryos.Mat'];
load(CompiledEmbryoPath);

DorsalNuclearPath = [liveExperiment.resultsFolder, filesep, 'DorsalNuclearProfiles.mat'];
    load(DorsalNuclearPath, 'DorsalAvgAPProfiles', 'DorsalAvgNarrowAPProfiles', 'NarrowProfileNucleiFluoInfo',...
    'DorsalAPProfiles', 'DorsalNarrowAPProfiles','ProfileNucleiFluoInfo',...
    'ProfileNarrowNucleiFluoInfo','AllDorsalNucleiFluoInfo','DorsalStdAPProfiles',...
    'DorsalStdNarrowAPProfiles','DorsalCountAPProfiles','DorsalCountNarrowAPProfiles');


CompiledEmbryos.DorsalAPProfiles = DorsalAPProfiles;


CompiledEmbryos.DorsalAvgAPProfiles = DorsalAvgAPProfiles;


CompiledEmbryos.NarrowProfileNucleiFluoInfo = NarrowProfileNucleiFluoInfo;

CompiledEmbryos.DorsalNarrowAPProfiles = DorsalNarrowAPProfiles;


CompiledEmbryos.DorsalAvgNarrowAPProfiles = DorsalAvgNarrowAPProfiles;
CompiledEmbryos.ProfileNucleiFluoInfo = ProfileNucleiFluoInfo;
CompiledEmbryos.ProfileNarrowNucleiFluoInfo = ProfileNarrowNucleiFluoInfo;
CompiledEmbryos.AllDorsalNucleiFluoInfo = AllDorsalNucleiFluoInfo;
CompiledEmbryos.DorsalStdAPProfiles = DorsalStdAPProfiles;
CompiledEmbryos.DorsalStdNarrowAPProfiles = DorsalStdNarrowAPProfiles;
CompiledEmbryos.DorsalCountAPProfiles = DorsalCountAPProfiles;
CompiledEmbryos.DorsalCountNarrowAPProfiles = DorsalCountNarrowAPProfiles;

save(CompiledEmbryoPath, 'CompiledEmbryos')
