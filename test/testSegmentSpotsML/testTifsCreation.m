function testCase = testTifsCreation(testCase)
  tic;
  disp(['Running segmentSpotsML to generate Tifs test with prefix ', testCase.Prefix]);
  fprintf('Test run started at %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  
  %Figure out the initial folders. 
  CONFIG_CSV_PATH = ['ComputerFolders.csv'];
  
  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');
  testPath = getConfigValue(configValues, 'TestPath');
  
  dynamicResultsPath = getConfigValue(configValues, 'DropboxFolder');
  preprocessedDataPath = getConfigValue(configValues, 'PreProcPath');
  
  expectedDataSubFolder = ['SegmentSpotsML', filesep, 'Tifs'];

  dynamicResultsExperimentPath = [dynamicResultsPath, filesep, testCase.Prefix];
  
  preprocessedDataExperimentPath = [preprocessedDataPath, filesep, testCase.Prefix];
  expectedPreProcessedDataFolder = [testPath, filesep, expectedDataSubFolder, filesep, 'PreProcessedData',...
    filesep, testCase.Prefix]; 

  % Clean up previous runs
  deleteDirectory(dynamicResultsExperimentPath, testCase.Prefix);
  deleteDirectory(preprocessedDataExperimentPath, testCase.Prefix);

  % Precondition - Run ExportsDataForFISH without deleting TIFs
  ExportDataForFISH(testCase.Prefix, 'keepTifs');
  
  % Tests Tifs generation
  segmentSpotsML(testCase.Prefix, [], 'Tifs');

  % Then verifies log.mat and FrameInfo.mat and expected contents of dogs folder
  assertFrameInfoEqualToExpected(testCase, dynamicResultsExperimentPath, testPath, expectedDataSubFolder);
  compareExpectedDataDir(testCase, preprocessedDataExperimentPath, expectedPreProcessedDataFolder);

  elapsedTime = toc;
  fprintf('Test run for %s ended successfully at %s\n', testCase.Prefix, datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  fprintf('Elapsed time for test was %d minutes and %f seconds\n', floor(elapsedTime/60), rem(elapsedTime,60));
end

