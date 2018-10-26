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
  segmentSpotsML(testCase.Prefix, [], classifierForTest, 'ignoreMemoryCheck');
  
  % Then verifies FrameInfo.mat and expected contents of dogs folder
  expectedPathSubFolderFilter = ['SegmentSpotsML', filesep, 'FilterMovie'];

  assertDogsFolderEqualToExpected(testCase, processedDataExperimentPath, testPath, expectedPathSubFolderFilter);
  
  % Tests second pass
  % Executes segment spots with known DoG
  segmentSpotsML(testCase.Prefix, testCase.Threshold);

  expectedPathSubfolderSpots = ['SegmentSpotsML', filesep, 'SegmentSpotsML'];

  assertLogFileExists(testCase, dynamicResultsExperimentPath);
  assertDogsFolderEqualToExpected(testCase, processedDataExperimentPath, testPath, expectedPathSubfolderSpots);
  assertSpotsEqualToExpected(testCase, dynamicResultsExperimentPath, testPath, expectedPathSubfolderSpots);

  % This will come in handy after the refactor
  % testSegmentSpotsNoThreshold(testCase);
  % testSegmentSpotsEmptyThreshold(testCase);

  elapsedTime = toc;
  fprintf('Test run for %s ended successfully at %s\n', testCase.Prefix, datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  fprintf('Elapsed time for test was %d minutes and %f seconds\n', floor(elapsedTime/60), rem(elapsedTime,60));
end
