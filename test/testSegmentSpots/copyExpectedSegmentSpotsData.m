function copyExpectedSegmentSpotsData(Prefix, testPath, dynamicsResultsPath, preprocessedDataPath)
  dynamicsResultsExpectedPath = [testPath, '/SegmentSpots/DynamicsResults/', Prefix];
  validateExpectedDataFolderExists(dynamicsResultsExpectedPath, Prefix)

  preProcessedDataExpectedPath = [testPath, '/SegmentSpots/PreProcessedData/', Prefix];
  validateExpectedDataFolderExists(preProcessedDataExpectedPath, Prefix)

  deleteDirectory(dynamicsResultsPath, Prefix);
  mkdir(dynamicsResultsPath);

  deleteDirectory(preprocessedDataPath, Prefix);
  mkdir(preprocessedDataPath);

  copyExpectedDataToData(Prefix, dynamicsResultsPath, dynamicsResultsExpectedPath, preprocessedDataPath,...
    preProcessedDataExpectedPath, [], []);
end
