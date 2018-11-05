function copyExpectedSegmentSpotsMLData(Prefix, testPath, dynamicsResultsPath, preprocessedDataPath, processedDataPath)
  preProcessedDataExpectedPath = getExpectedDataFolderNameSegmentSpotsML(testPath, 'Tifs', 'PreProcessedData', Prefix);
  validateExpectedDataFolderExists(preProcessedDataExpectedPath, Prefix)

  dynamicsResultsExpectedPath = getExpectedDataFolderNameSegmentSpotsML(testPath, 'SegmentSpotsML', 'DynamicsResults', Prefix);
  validateExpectedDataFolderExists(dynamicsResultsExpectedPath, Prefix)

  processedDataExpectedPath = getExpectedDataFolderNameSegmentSpotsML(testPath, 'SegmentSpotsML', 'ProcessedData', Prefix);
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
