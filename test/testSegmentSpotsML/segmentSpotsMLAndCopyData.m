function segmentSpotsMLAndCopyData(testCase, testPath, dynamicsResultsPath, preprocessedDataPath, processedDataPath)
  expectedDataFolder = [testPath, '/SegmentSpotsML'];
  if 7 ~= exist(expectedDataFolder, 'dir')
    mkdir(expectedDataFolder);
  end

  dynamicsResultsExpectedPath = [expectedDataFolder, '/DynamicsResults/', testCase.Prefix];
  preProcessedDataExpectedPath = [expectedDataFolder, '/PreProcessedData/', testCase.Prefix];
  processedDataExpectedPath = [expectedDataFolder, '/ProcessedData/', testCase.Prefix, '_'];

  % Executes segment spots ML with known DoG
  segmentSpotsML(testCase.Prefix, testCase.Threshold);

  copyDataToExpectedData(testCase.Prefix, dynamicsResultsPath, dynamicsResultsExpectedPath, preprocessedDataPath,...
    preProcessedDataExpectedPath, processedDataPath, processedDataExpectedPath);

  disp(['Test data for segment spots completed for Prefix ', testCase.Prefix]);
end
