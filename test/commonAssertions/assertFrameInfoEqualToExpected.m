% Given the path of the dynamic results folder and the expected data path (1st pass or 2nd pass),
% asserts that FrameInfo.mat is equal to expected.
function assertFrameInfoEqualToExpected(testCase, dynamicResultsExperimentPath, testPath, expectedDataSubFolder)
  frameInfoPath = [dynamicResultsExperimentPath, filesep, 'FrameInfo.mat'];
  testCase.assertTrue(exist(frameInfoPath, 'file') == 2);
  frameInfoResult = load(frameInfoPath);

  expectedDataFolder = [testPath, filesep, expectedDataSubFolder, filesep, 'DynamicsResults',...
    filesep, testCase.Prefix];
  expectedFrameInfoPath = [expectedDataFolder, filesep, 'FrameInfo.mat'];
  expectedFrameInfo = load(expectedFrameInfoPath);

  testCase.assertEqual(frameInfoResult.FrameInfo, expectedFrameInfo.FrameInfo);
end
