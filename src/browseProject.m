function browseProject(Project)

[~, experiment, ~] = LoadMS2Sets(Project, 'justPrefixes', 'noCompiledNuclei');

for i = 1:length(experiment)
    
     [rawExportedDirectory,...
    processedDirectory,...
    resultsDirectory] = browseExperiment(experiment{i});
    
end


end