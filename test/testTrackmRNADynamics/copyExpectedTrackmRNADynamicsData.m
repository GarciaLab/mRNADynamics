function copyExpectedTrackmRNADynamicsData(Prefix, testPath, dynamicsResultsPath, preprocessedDataPath, processedDataPath)
  preProcessedDataExpectedPath = [testPath, filesep, 'TrackmRNADynamics', filesep, 'PreProcessedData', filesep, Prefix];
  validateExpectedDataFolderExists(preProcessedDataExpectedPath, Prefix)

  dynamicsResultsExpectedPath = [testPath, filesep, 'TrackmRNADynamics', filesep, 'DynamicsResults', filesep, Prefix]; 
  validateExpectedDataFolderExists(dynamicsResultsExpectedPath, Prefix)

  processedDataExpectedPath = [testPath, filesep, 'TrackmRNADynamics', filesep, 'ProcessedData', filesep, Prefix, '_'];
  validateExpectedDataFolderExists(processedDataExpectedPath, Prefix)

  deleteDirectory(dynamicsResultsPath, Prefix);
  mkdir(dynamicsResultsPath);

  deleteDirectory(preprocessedDataPath, Prefix);
  mkdir(preprocessedDataPath);

  deleteDirectory(processedDataExpectedPath, Prefix);
  mkdir(processedDataExpectedPath);

  copyExpectedDataToData(Prefix, dynamicsResultsPath, dynamicsResultsExpectedPath, preprocessedDataPath,...
    preProcessedDataExpectedPath, processedDataPath, processedDataExpectedPath);
end
