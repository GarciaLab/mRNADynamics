% Runs the whole pipeline or just segment spots before running the test for trackmRNADynamics
% in order to allow testing for both paths, histone and not histone.
function testCase = testTrackmRNADynamicsBaseTest(testCase, histone)
  tic;
  disp(['Running trackmRNADynamics test with prefix ', testCase.Prefix]);
  fprintf('Test run started at %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  
  %Figure out the initial folders. 
  CONFIG_CSV_PATH = ['ComputerFolders.csv'];
  
  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');
  testPath = getConfigValue(configValues, 'TestPath');

  dynamicResultsPath = getConfigValue(configValues, 'DropboxFolder');

  dynamicResultsExperimentPath = [dynamicResultsPath, filesep, testCase.Prefix];
  
  % Precondition, copies existing Expected Data to proper folders before running the process
  % Depending if the code should use histone or not, it runs track nucleir, or just segment spots
  TestFolder = 'TrackmRNADynamics';
  if (histone)
    copyExpectedDataForPrefix(testCase.Prefix, 'TrackNuclei');
  else 
    copyExpectedDataForPrefix(testCase.Prefix, 'SegmentSpots');
    TestFolder = 'TrackmRNADynamicsNoHistone';
  end

  % Executes trackmRNADynamis with known thresholds
  TrackmRNADynamics(testCase.Prefix);

  expectedDynamicsResultsFolder = [testPath, filesep, TestFolder, filesep, 'DynamicsResults', filesep, testCase.Prefix];
  assertStructEqualToExpected(testCase, dynamicResultsExperimentPath, expectedDynamicsResultsFolder, 'Particles.mat');
  assertStructEqualToExpected(testCase, dynamicResultsExperimentPath, expectedDynamicsResultsFolder, 'FrameInfo.mat');

  % Executes tracking a second time to test out the case when Particles.mat already exists
  TrackmRNADynamics(testCase.Prefix);

  expectedDynamicsResultsFolder = [testPath, filesep, TestFolder, filesep, 'secondPass', filesep, 'DynamicsResults',...
    filesep, testCase.Prefix];
  assertStructEqualToExpected(testCase, dynamicResultsExperimentPath, expectedDynamicsResultsFolder, 'Particles.mat');
  assertStructEqualToExpected(testCase, dynamicResultsExperimentPath, expectedDynamicsResultsFolder, 'FrameInfo.mat');

  elapsedTime = toc;
  fprintf('Test run for %s ended successfully at %s\n', testCase.Prefix, datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  fprintf('Elapsed time for test was %d minutes and %f seconds\n', floor(elapsedTime/60), rem(elapsedTime,60));
end
