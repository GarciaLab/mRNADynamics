function copyFilterMovieTifsData(Prefix, testPath, dynamicsResultsPath, preprocessedDataPath)
  dynamicsResultsExpectedPath = [testPath, '/filterMovieTifs/DynamicsResults/', Prefix];
  validateExpectedDataFolderExists(dynamicsResultsExpectedPath, Prefix)

  preProcessedDataExpectedPath = [testPath, '/filterMovieTifs/PreProcessedData/', Prefix];
  validateExpectedDataFolderExists(preProcessedDataExpectedPath, Prefix)

  deleteDirectory(dynamicsResultsPath, Prefix);
  mkdir(dynamicsResultsPath);

  deleteDirectory(preprocessedDataPath, Prefix);
  mkdir(preprocessedDataPath);

  copyExpectedDataToData(Prefix, dynamicsResultsPath, dynamicsResultsExpectedPath, preprocessedDataPath,... 
    preProcessedDataExpectedPath, [], []);
end
