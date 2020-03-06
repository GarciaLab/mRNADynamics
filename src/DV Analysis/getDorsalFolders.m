function [dataFolder, resultsFolder] = getDorsalFolders()

    configValues = csv2cell('ComputerFolders.csv', 'fromfile');
    dataFolder = getConfigValue(configValues, 'DorsalSynthetics');
    resultsFolder = getConfigValue(configValues, 'DorsalSyntheticsDropbox');

end
