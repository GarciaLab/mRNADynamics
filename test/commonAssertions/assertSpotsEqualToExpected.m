% Given the path of the dynamic results folder and the expected data path (1st pass or 2nd pass),
% asserts that FrameInfo.mat is equal to expected.
function assertSpotsEqualToExpected(testCase, dynamicResultsExperimentPath, testPath, expectedDataSubFolder)
  spotsPath = [dynamicResultsExperimentPath, filesep, 'Spots.mat'];
  testCase.assertTrue(exist(spotsPath, 'file') == 2);
  spotsFile = load(spotsPath);
  
  expectedDataFolder = [testPath, filesep, expectedDataSubFolder, filesep, 'DynamicsResults',...
    filesep, testCase.Prefix];
  expectedSpotsPath = [expectedDataFolder, filesep, 'Spots.mat'];
  expectedSpotsFile = load(expectedSpotsPath);

  testCase.assertEqual(spotsFile.Spots, expectedSpotsFile.Spots);
end
