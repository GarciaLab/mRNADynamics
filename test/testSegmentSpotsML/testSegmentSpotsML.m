function testCase = testSegmentSpotsML(testCase)
  tic;
  disp(['Running segmentSpotsML test with prefix ', testCase.Prefix]);
  fprintf('Test run started at %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  
  %Figure out the initial folders. 
  CONFIG_CSV_PATH = ['ComputerFolders.csv'];
  
  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');
  testPath = getConfigValue(configValues, 'TestPath');

  dynamicResultsPath = getConfigValue(configValues, 'DropboxFolder');

  dynamicResultsExperimentPath = [dynamicResultsPath, filesep, testCase.Prefix];
  
  % Precondition, copies existing Expected Data to proper folders before running the process
  copyExpectedDataForPrefix(testCase.Prefix, 'filterMovieWeka');

  % Executes segment spots with known DoG
  segmentSpotsML(testCase.Prefix, testCase.Threshold);

  assertLogFileExists(testCase, dynamicResultsExperimentPath);
  assertSpotsEqualToExpected(testCase, dynamicResultsExperimentPath, testPath, 'SegmentSpotsML');

  testSegmentSpotsMLTifsUnsupported(testCase);
  testSegmentSpotsMLNoThreshold(testCase);
  testSegmentSpotsMLEmptyThreshold(testCase);

  elapsedTime = toc;
  fprintf('Test run for %s ended successfully at %s\n', testCase.Prefix, datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  fprintf('Elapsed time for test was %d minutes and %f seconds\n', floor(elapsedTime/60), rem(elapsedTime,60));
end

function testCase = testSegmentSpotsMLTifsUnsupported(testCase) 
  try
    segmentSpotsML(testCase.Prefix, testCase.Threshold, 'Tifs');
  catch ME
    testCase.assertEqual(ME.message,...
      'Tifs generation is no longer supported from segmentSpotsML, try filterMovie(Prefix, ''Tifs'') instead.');
  end
end

function testCase = testSegmentSpotsMLNoThreshold(testCase) 
  try
    segmentSpotsML(testCase.Prefix);
  catch ME
    testCase.assertEqual(ME.message,['Not enough input arguments.']);
  end
end

function testCase = testSegmentSpotsMLEmptyThreshold(testCase) 
  try
    segmentSpotsML(testCase.Prefix, []);
  catch ME
    testCase.assertEqual(ME.message,['Please use filterMovie(Prefix, ''Weka'', options) instead ',...
      'of segmentSpots with the argument "[]" to generate DoG images']);
  end
end

function copyExpectedData(Prefix, testPath, dynamicResultsExperimentPath, preprocessedDataExperimentPath, processedDataExperimentPath)
  if 7 ~= exist(testPath, 'dir')
    error('Test path data does not exist, please generate test data before running test case.');
  end

  mkdir(dynamicResultsExperimentPath);
  mkdir(preprocessedDataExperimentPath);
  mkdir(processedDataExperimentPath);

  dynamicsResultsDataExpectedFolder = [testPath, filesep, 'ExportDataForFISH', filesep,...
    'DynamicsResults', filesep, Prefix];
  preProcessedDataExpectedFolder = [testPath, filesep, 'SegmentSpotsML', filesep, 'Tifs', filesep,...
    'PreProcessedData', filesep, Prefix];
  processedDataExpectedFolder = [testPath, filesep, 'SegmentSpotsML', filesep, 'FilterMovie', filesep,...
    'ProcessedData', filesep, Prefix, '_'];

  copyfile([dynamicsResultsDataExpectedFolder, filesep, '*'], dynamicResultsExperimentPath);
  copyfile([preProcessedDataExpectedFolder, filesep, '*'], preprocessedDataExperimentPath);
  copyfile([processedDataExpectedFolder, filesep, '*'], processedDataExperimentPath);
end
