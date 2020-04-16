function dorsalResults = compileAllProjects(DataType)


thisProject = liveProject(DataType);


%% Validate experiments included in the analysis
hasAllPushed = [thisProject.hasSpots, thisProject.hasParticles,...
    thisProject.hasSchnitzcells, thisProject.hasCompiledParticles,...
    thisProject.hasAnaphaseFramesAnnotated];

assert( all(hasAllPushed) ); 

%%

prefixes = thisProject.includedExperimentNames;

compiledProjects = cell(1, length(prefixes));

for k = 1:length(prefixes)
    TrackmRNADynamics(prefixes{k}, 'noretrack');
    CompileParticles(prefixes{k},  'minBinSize', 0, 'MinParticles', 0,...
        'yToManualAlignmentPrompt');
    alignCompiledParticlesByAnaphase(prefixes{k});
end

addDVStuffToSchnitzCells(DataType)

binDorsal(DataType, false)

for k = 1:length(prefixes)
    compiledProjects{k} = makeCompiledProject(prefixes{k});
end

dorsalResults = plotFracByDlFluo2(DataType); 

activity = '';
plotDorsalResultsLoop(DataType, activity)

end