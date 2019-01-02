function testCase = testTifsCreation(testCase)
  tic;
  disp(['Running segmentSpotsML to generate Tifs test with prefix ', testCase.Prefix]);
  fprintf('Test run started at %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  
  %Figure out the initial folders. 
  CONFIG_CSV_PATH = ['ComputerFolders.csv'];
  
  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');
  testPath = getConfigValue(configValues, 'TestPath');
  
  preprocessedDataPath = getConfigValue(configValues, 'PreProcPath');

  preprocessedDataExperimentPath = [preprocessedDataPath, filesep, testCase.Prefix];
  expectedPreProcessedDataFolder = [testPath, '/filterMovieTifs/PreProcessedData/', testCase.Prefix]; 

  % Precondition, copies existing Expected Data to proper folders before running the process
  copyExpectedDataForPrefix(testCase.Prefix, 'ExportDataForLivemRNA');
  
  % Tests Tifs generation
  filterMovie(testCase.Prefix, 'Tifs');

  % Then verifies log.mat and FrameInfo.mat and expected contents of dogs folder
  compareExpectedDataDir(testCase, preprocessedDataExperimentPath, expectedPreProcessedDataFolder);

  elapsedTime = toc;
  fprintf('Test run for %s ended successfully at %s\n', testCase.Prefix, datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  fprintf('Elapsed time for test was %d minutes and %f seconds\n', floor(elapsedTime/60), rem(elapsedTime,60));
end

