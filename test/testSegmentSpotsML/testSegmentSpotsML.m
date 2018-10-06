function testCase = testSegmentSpotsML(testCase)
  tic;
  disp(['Running segmentSpotsML test with prefix ', testCase.Prefix]);
  fprintf('Test run started at %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  
  %Figure out the initial folders. 
  CONFIG_CSV_PATH = ['ComputerFolders.csv'];
  
  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');
  
  codePath = getConfigValue(configValues, 'MS2CodePath');
  testPath = getConfigValue(configValues, 'TestPath');

  dynamicResultsPath = getConfigValue(configValues, 'DropboxFolder');
  processedDataPath = getConfigValue(configValues, 'FISHPath');

  % Classifier path must end with filesep
  classifiersPath = strcat(codePath, filesep, 'src', filesep, 'classifiers', filesep);
  classifierForTest = ClassifierForTest(classifiersPath, testCase.classifier);

  dynamicResultsExperimentPath = [dynamicResultsPath, filesep, testCase.Prefix];
  % harrypotel: Not sure why the folder on ProcessedData has a _ at the end, is it a bug?
  processedDataExperimentPath = [processedDataPath, filesep, testCase.Prefix, '_'];
  
  % Clean up previous runs
  deleteDirectory(dynamicResultsExperimentPath, testCase.Prefix);
  deleteDirectory(processedDataExperimentPath, testCase.Prefix);
  
  % Precondition - Run ExportsDataForFISH without deleting TIFs
  ExportDataForFISH(testCase.Prefix, 'keepTifs');
  
  % Tests first pass
  % Generates DoGs
  filterMovie(testCase.Prefix, 'Weka', classifierForTest, 'ignoreMemoryCheck');
  
  % Then verifies FrameInfo.mat and expected contents of dogs folder
  expectedPathSubFolderFilter = ['SegmentSpotsML', filesep, 'FilterMovie'];
  assertLogFileExists(testCase, dynamicResultsExperimentPath);
  assertDogsFolderEqualToExpected(testCase, processedDataExperimentPath, testPath, expectedPathSubFolderFilter);
  
  % Tests second pass
  % Executes segment spots with known DoG
  segmentSpotsML(testCase.Prefix, testCase.Threshold);

  expectedPathSubfolderSpots = ['SegmentSpotsML', filesep, 'SegmentSpotsML'];

  assertLogFileExists(testCase, dynamicResultsExperimentPath);
  assertDogsFolderEqualToExpected(testCase, processedDataExperimentPath, testPath, expectedPathSubfolderSpots);
  assertSpotsEqualToExpected(testCase, dynamicResultsExperimentPath, testPath, expectedPathSubfolderSpots);

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

