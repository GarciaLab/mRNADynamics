function dorsalResults = compileAllProjects(DataType)

[~, ~, prefixes] = getDorsalPrefixes(DataType);

compiledProjects = cell(1, length(prefixes));
for k = 1:length(prefixes)
    CompileParticles(prefixes{k}, 'SkipAll', 'ApproveAll', 'minBinSize', 0, 'MinParticles', 0, 'yToManualAlignmentPrompt');
    addDVStuffToSchnitzCells(DataType)
    alignCompiledParticlesByAnaphase(prefixes{k});
    compiledProjects{k} = makeCompiledProject(prefixes{k});
end

binDorsal(DataType, true)

dorsalResults = plotFracByDlFluo2(DataType); 

activity = '';
plotDorsalResultsLoop(DataType, activity)

end