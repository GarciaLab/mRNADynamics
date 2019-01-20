function copyExpectedTrackmRNADynamicsData(Prefix, testPath, dynamicsResultsPath, preprocessedDataPath, processedDataPath)
  preProcessedDataExpectedPath = [testPath, filesep, 'TrackmRNADynamics', filesep, 'PreProcessedData', filesep, Prefix];
  validateExpectedDataFolderExists(preProcessedDataExpectedPath, Prefix)

  dynamicsResultsExpectedPath = [testPath, filesep, 'TrackmRNADynamics', filesep, 'DynamicsResults', filesep, Prefix]; 
  validateExpectedDataFolderExists(dynamicsResultsExpectedPath, Prefix)

  deleteDirectory(dynamicsResultsPath, Prefix);
  mkdir(dynamicsResultsPath);

  deleteDirectory(preprocessedDataPath, Prefix);
  mkdir(preprocessedDataPath);

  copyExpectedDataToData(Prefix, dynamicsResultsPath, dynamicsResultsExpectedPath, preprocessedDataPath,...
    preProcessedDataExpectedPath, [], []);
end
