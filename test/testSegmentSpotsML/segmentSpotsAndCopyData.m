function segmentSpotsAndCopyData(testCase, testPath, dynamicsResultsPath, preprocessedDataPath)
  expectedDataFolder = [testPath, '/SegmentSpots'];
  if 7 ~= exist(expectedDataFolder, 'dir')
    mkdir(expectedDataFolder);
  end

  dynamicsResultsExpectedPath = [expectedDataFolder, '/DynamicsResults/', testCase.Prefix];
  preProcessedDataExpectedPath = [expectedDataFolder, '/PreProcessedData/', testCase.Prefix];

  % Executes segment spots with known DoG
  segmentSpots(testCase.Prefix, testCase.Threshold);

  copyDataToExpectedData(testCase.Prefix, dynamicsResultsPath, dynamicsResultsExpectedPath, preprocessedDataPath,...
    preProcessedDataExpectedPath, [], []);

  disp(['Test data for segment spots completed for Prefix ', testCase.Prefix]);
end
