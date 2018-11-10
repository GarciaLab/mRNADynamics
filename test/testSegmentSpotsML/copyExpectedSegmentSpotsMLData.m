function copyExpectedSegmentSpotsMLData(Prefix, testPath, dynamicsResultsPath, preprocessedDataPath, processedDataPath)
  dynamicsResultsExpectedPath = [testPath, '/SegmentSpotsML/DynamicsResults/', Prefix];
  validateExpectedDataFolderExists(dynamicsResultsExpectedPath, Prefix)

  preProcessedDataExpectedPath = [testPath, '/SegmentSpotsML/PreProcessedData/', Prefix];
  validateExpectedDataFolderExists(preProcessedDataExpectedPath, Prefix)

  processedDataExpectedPath = [testPath, '/SegmentSpotsML/ProcessedData/', Prefix, '_'];
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
