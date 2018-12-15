function copyFilterMovieWekaData(Prefix, testPath, dynamicsResultsPath, preprocessedDataPath, processedDataPath)
  dynamicsResultsExpectedPath = [testPath, '/filterMovieTifs/DynamicsResults/', Prefix];
  validateExpectedDataFolderExists(dynamicsResultsExpectedPath, Prefix)

  preProcessedDataExpectedPath = [testPath, '/filterMovieTifs/PreProcessedData/', Prefix];
  validateExpectedDataFolderExists(preProcessedDataExpectedPath, Prefix)

  processedDataExpectedPath = [testPath, '/filterMovieWeka/ProcessedData/', Prefix, '_'];
  validateExpectedDataFolderExists(processedDataExpectedPath, Prefix)

  deleteDirectory(dynamicsResultsPath, Prefix);
  mkdir(dynamicsResultsPath);

  deleteDirectory(preprocessedDataPath, Prefix);
  mkdir(preprocessedDataPath);

  deleteDirectory(processedDataPath, Prefix);
  mkdir(processedDataPath);

  copyExpectedDataToData(Prefix, dynamicsResultsPath, dynamicsResultsExpectedPath, preprocessedDataPath,... 
    preProcessedDataExpectedPath, processedDataPath, processedDataExpectedPath);
end
