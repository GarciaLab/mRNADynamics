function testCase = testTrackmRNADynamics(testCase)
  tic;
  disp(['Running trackmRNADynamics test with prefix ', testCase.Prefix]);
  fprintf('Test run started at %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  
  %Figure out the initial folders. 
  CONFIG_CSV_PATH = ['ComputerFolders.csv'];
  
  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');
  testPath = getConfigValue(configValues, 'TestPath');

  dynamicResultsPath = getConfigValue(configValues, 'DropboxFolder');
  preprocessedDataPath = getConfigValue(configValues, 'PreProcPath');
  processedDataPath = getConfigValue(configValues, 'FISHPath');

  dynamicResultsExperimentPath = [dynamicResultsPath, filesep, testCase.Prefix];
  preprocessedDataExperimentPath = [preprocessedDataPath, filesep, testCase.Prefix];
  % harrypotel: Not sure why the folder on ProcessedData has a _ at the end, is it a bug?
  processedDataExperimentPath = [processedDataPath, filesep, testCase.Prefix, '_'];
  
  % Precondition, copies existing Expected Data to proper folders before running the process
  copyExpectedDataForPrefix(testCase.Prefix, 'TrackNuclei');

  % Executes trackmRNADynamis with known thresholds
  TrackmRNADynamics(testCase.Prefix, testCase.Threshold1, testCase.Threshold1, 'bypassUserPrompt');

  expectedDynamicsResultsFolder = [testPath, filesep, 'TrackmRNADynamics', filesep, 'DynamicsResults', filesep, testCase.Prefix];
  assertStructEqualToExpected(testCase, dynamicResultsExperimentPath, expectedDynamicsResultsFolder, 'Particles.mat');
  assertStructEqualToExpected(testCase, dynamicResultsExperimentPath, expectedDynamicsResultsFolder, 'FrameInfo.mat');

  % Executes tracking a second time to test out the case when Particles.mat already exists
  TrackmRNADynamics(testCase.Prefix, testCase.Threshold1, testCase.Threshold1, 'bypassUserPrompt');

  expectedDynamicsResultsFolder = [testPath, filesep, 'TrackmRNADynamics', filesep, 'secondPass', filesep, 'DynamicsResults',...
    filesep, testCase.Prefix];
  assertStructEqualToExpected(testCase, dynamicResultsExperimentPath, expectedDynamicsResultsFolder, 'Particles.mat');
  assertStructEqualToExpected(testCase, dynamicResultsExperimentPath, expectedDynamicsResultsFolder, 'FrameInfo.mat');

  elapsedTime = toc;
  fprintf('Test run for %s ended successfully at %s\n', testCase.Prefix, datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  fprintf('Elapsed time for test was %d minutes and %f seconds\n', floor(elapsedTime/60), rem(elapsedTime,60));
end
