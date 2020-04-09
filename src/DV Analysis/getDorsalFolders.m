function [dataFolder, resultsFolder, PreProcessedFolder, ProcessedFolder] = getDorsalFolders()

    configValues = csv2cell('ComputerFolders.csv', 'fromfile');
    dataFolder = getConfigValue(configValues, 'DorsalSynthetics');
    resultsFolder = getConfigValue(configValues, 'DorsalSyntheticsDropbox');
    PreProcessedFolder = [dataFolder, filesep, 'PreProcessedData'];
    ProcessedFolder = [dataFolder, filesep, 'ProcessedData'];
    
end
