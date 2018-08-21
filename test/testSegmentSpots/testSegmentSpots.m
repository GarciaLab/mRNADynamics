function testCase = testSegmentSpots(testCase)
  tic;
  disp(['Running with segment spot test with prefix ', testCase.Prefix, ' and DoG ', testCase.DoG]);
  fprintf('Test run started at %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  
  %Figure out the initial folders. 
  CONFIG_CSV_PATH = ['ComputerFolders.csv'];
  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');
  testPath = getConfigValue(configValues, 'TestPath');
  
  dynamicResultsPath = getConfigValue(configValues, 'DropboxFolder');
  processedDataPath = getConfigValue(configValues, 'FISHPath');
  
  dynamicResultsExperimentPath = [dynamicResultsPath, filesep, testCase.Prefix];
  % harrypotel: Not sure why the folder on ProcessedData has a _ at the end, is it a bug?
  processedDataExperimentPath = [processedDataPath, filesep, testCase.Prefix, '_'];

  % Clean up previous runs
  deleteDirectory(dynamicResultsExperimentPath);
  deleteDirectory(processedDataExperimentPath);

  % Precondition - Run ExportsDataForFISH
  ExportDataForFISH(testCase.Prefix);
  
  % Tests first pass
  % Executes segment spots without specifying DoG
  segmentSpots(testCase.Prefix, []);
  
  % Then verifies log.mat and FrameInfo.mat and expected contents of dogs folder
  assertLogFileExists(testCase, dynamicResultsExperimentPath);
  assertFrameInfoEqualToExpected(testCase, dynamicResultsExperimentPath, testPath, 'SegmentSpots_1stPass');
  assertDogsFolderEqualToExpected(testCase, processedDataExperimentPath, testPath, 'SegmentSpots_1stPass');

  % Tests second pass
  % Executes segment spots with known DoG
  segmentSpots(testCase.Prefix, testCase.DoG);

  assertLogFileExists(testCase, dynamicResultsExperimentPath);
  assertFrameInfoEqualToExpected(testCase, dynamicResultsExperimentPath, testPath, 'SegmentSpots_2ndPass');
  assertDogsFolderEqualToExpected(testCase, processedDataExperimentPath, testPath, 'SegmentSpots_2ndPass');
  assertSpotsEqualToExpected(testCase, dynamicResultsExperimentPath, testPath, 'SegmentSpots_2ndPass');

  elapsedTime = toc;
  fprintf('Test run for %s ended successfully at %s\n', testCase.Prefix, datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  fprintf('Elapsed time for test was %d minutes and %f seconds\n', floor(elapsedTime/60), rem(elapsedTime,60));
end

% Given a path, asserts that a log.mat file exists on it.
function  assertLogFileExists(testCase, dynamicResultsExperimentPath);
  logFilePath = [dynamicResultsExperimentPath, filesep, 'log.mat'];
  testCase.assertTrue(exist(logFilePath, 'file') == 2);
end

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

% Given the path of the processed data folderand the expected data path (1st pass or 2nd pass),
% asserts dogs folder with expected data set.
function assertDogsFolderEqualToExpected(testCase, processedDataExperimentPath, testPath, expectedDataSubFolder);
  dogsPath = [processedDataExperimentPath, filesep, 'dogs'];
  expectedDataFolder = [testPath, filesep, expectedDataSubFolder, filesep, 'ProcessedData',...
    filesep, testCase.Prefix, '_', filesep, 'dogs']; 
  compareExpectedDataDir(testCase, dogsPath, expectedDataFolder);
end

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
