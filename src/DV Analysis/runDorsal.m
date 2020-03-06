
dataType = '1DG';
[allData, Prefixes, resultsFolder] = LoadMS2Sets(dataType, 'justPrefixes');
% parpool(8)
for i= 1:length(Prefixes)
    Prefixes{i}
    i
    fit3DGaussiansToAllSpots(Prefixes{i})
%     CompileNuclearProtein(Prefixes{i});
%     plotByDorsalConc(Prefixes{i});
    AddParticlePosition(Prefixes{i},  'yToManualAlignmentPrompt');

    CompileParticles(Prefixes{i}, 'SkipAll', 'ApproveAll', 'minBinSize', .3, 'MinParticles', 0, 'yToManualAlignmentPrompt');

end
% [allData, Prefixes, resultsFolder] = LoadMS2Sets(dataType);
% averageAcrossEmbryosDV(allData);
% AverageDatasets_syntheticsDV(dataType, 'savePath','[E:\Armando\LivemRNA\Data\Dropbox\DorsalSyntheticsDropbox\analysis\figure\',dataType]);
% plotFractionOnDV(['E:\Armando\LivemRNA\Data\Dropbox\DorsalSyntheticsDropbox\analysis\figure\',dataType,'\',dataType,'.mat'])
addDVStuffToSchnitzCells(dataType)
plotFracByDlFluo(dataType)