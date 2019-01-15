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

  % Switches to a different directory so the removal does not fail
  cd(testPath);
  % Clean up previous runs
  deleteDirectory(dynamicResultsExperimentPath, testCase.Prefix);
  deleteDirectory(processedDataExperimentPath, testCase.Prefix);

  % Precondition - Run ExportsDataForFISH without deleting TIFs
  ExportDataForLivemRNA(testCase.Prefix, 'keepTifs');
  
  % Tests first pass
  % Generates DoGs
  filterMovie(testCase.Prefix);
  
  % Then verifies log.mat and FrameInfo.mat and expected contents of dogs folder
  assertLogFileExists(testCase, dynamicResultsExperimentPath);
  assertDogsFolderEqualToExpected(testCase, processedDataExperimentPath, testPath, 'SegmentSpots/1stPass');

  % Tests second pass
  % Executes segment spots with known DoG
  segmentSpots(testCase.Prefix, testCase.DoG);

  assertLogFileExists(testCase, dynamicResultsExperimentPath);
  assertSpotsEqualToExpected(testCase, dynamicResultsExperimentPath, testPath, 'SegmentSpots/2ndPass');

  testSegmentSpotsNoThreshold(testCase);
  testSegmentSpotsEmptyThreshold(testCase);

  elapsedTime = toc;
  fprintf('Test run for %s ended successfully at %s\n', testCase.Prefix, datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  fprintf('Elapsed time for test was %d minutes and %f seconds\n', floor(elapsedTime/60), rem(elapsedTime,60));
end

function testCase = testSegmentSpotsNoThreshold(testCase) 
  try
    segmentSpots(testCase.Prefix);
  catch ME
    testCase.assertEqual(ME.message,...
     'Please use filterMovie(Prefix, options) instead of segmentSpots with the argument "[]" to generate DoG images');
  end
end

function testCase = testSegmentSpotsEmptyThreshold(testCase) 
  try
    segmentSpots(testCase.Prefix, []);
  catch ME
    testCase.assertEqual(ME.message,...
     'Please use filterMovie(Prefix, options) instead of segmentSpots with the argument "[]" to generate DoG images');
  end
end
