function folder = findExperimentFolder(Prefix)

liveExperiment = LiveExperiment(Prefix);

if isfolder([liveExperiment.userExperimentsFolder, filesep, Prefix])
    
    preFolder = '';
    procFolder = '';
    rawFolder = '';
    
end