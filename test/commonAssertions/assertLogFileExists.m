% Given a path, asserts that a log.mat file exists on it.
function  assertLogFileExists(testCase, dynamicResultsExperimentPath)

  logFilePath = [dynamicResultsExperimentPath, filesep, 'log.mat'];
  testCase.assertTrue(exist(logFilePath, 'file') == 2);

end
