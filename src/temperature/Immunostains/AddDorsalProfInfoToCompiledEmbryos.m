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


% 
% DorsalNuclearPath = [liveExperiment.resultsFolder, filesep, 'EdgeCorrectedDorsalNuclearProfiles.mat'];
% load(DorsalNuclearPath, 'DorsalAvgAPProfiles', 'DorsalAvgNarrowAPProfiles', 'NarrowProfileNucleiFluoInfo',...
%     'DorsalAPProfiles', 'DorsalNarrowAPProfiles','ProfileNucleiFluoInfo',...
%     'ProfileNarrowNucleiFluoInfo','AllDorsalNucleiFluoInfo','DorsalStdAPProfiles',...
%     'DorsalStdNarrowAPProfiles','DorsalCountAPProfiles','DorsalCountNarrowAPProfiles',...
%     'MiddleMeanDorsalAvgAPProfiles', 'MiddleMeanDorsalAvgNarrowAPProfiles', 'MiddleMeanDorsalAPProfiles',...
%     'MiddleMeanDorsalNarrowAPProfiles','MiddleMeanDorsalStdAPProfiles','MiddleMeanDorsalStdNarrowAPProfiles',...
%     'MiddleMeanDorsalCountAPProfiles','MiddleMeanDorsalCountNarrowAPProfiles');
% 
% 
% CompiledEmbryos.EdgeCorrectedDorsalAPProfiles = DorsalAPProfiles;
% CompiledEmbryos.EdgeCorrectedDorsalAvgAPProfiles = DorsalAvgAPProfiles;
% CompiledEmbryos.EdgeCorrectedDorsalNarrowAPProfiles = DorsalNarrowAPProfiles;
% CompiledEmbryos.EdgeCorrectedDorsalAvgNarrowAPProfiles = DorsalAvgNarrowAPProfiles;
% CompiledEmbryos.EdgeCorrectedDorsalStdAPProfiles = DorsalStdAPProfiles;
% CompiledEmbryos.EdgeCorrectedDorsalStdNarrowAPProfiles = DorsalStdNarrowAPProfiles;
% CompiledEmbryos.EdgeCorrectedDorsalCountAPProfiles = DorsalCountAPProfiles;
% CompiledEmbryos.EdgeCorrectedDorsalCountNarrowAPProfiles = DorsalCountNarrowAPProfiles;
% 
% CompiledEmbryos.MiddleMeanDorsalAPProfiles = MiddleMeanDorsalAPProfiles;
% CompiledEmbryos.MiddleMeanDorsalAvgAPProfiles = MiddleMeanDorsalAvgAPProfiles;
% CompiledEmbryos.MiddleMeanDorsalNarrowAPProfiles = MiddleMeanDorsalNarrowAPProfiles;
% CompiledEmbryos.MiddleMeanDorsalAvgNarrowAPProfiles = MiddleMeanDorsalAvgNarrowAPProfiles;
% CompiledEmbryos.MiddleMeanDorsalStdAPProfiles = MiddleMeanDorsalStdAPProfiles;
% CompiledEmbryos.MiddleMeanDorsalStdNarrowAPProfiles = MiddleMeanDorsalStdNarrowAPProfiles;
% CompiledEmbryos.MiddleMeanDorsalCountAPProfiles = MiddleMeanDorsalCountAPProfiles;
% CompiledEmbryos.MiddleMeanDorsalCountNarrowAPProfiles = MiddleMeanDorsalCountNarrowAPProfiles;


save(CompiledEmbryoPath, 'CompiledEmbryos')
