function testCase = testGenerateDogsWeka(testCase)
  tic;
  disp(['Running genrate Dogs Weka test with prefix ', testCase.Prefix]);
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
  
  % Precondition, copies existing Expected Data to proper folders before running the process
  copyExpectedDataForPrefix(testCase.Prefix, 'filterMovieTifs');

  % Generates DoGs
  filterMovie(testCase.Prefix, 'Weka', classifierForTest, 'ignoreMemoryCheck');
  
  % Then verifies FrameInfo.mat and expected contents of dogs folder
  assertLogFileExists(testCase, dynamicResultsExperimentPath);
  assertDogsFolderEqualToExpected(testCase, processedDataExperimentPath, testPath, 'filterMovieWeka');

  elapsedTime = toc;
  fprintf('Test run for %s ended successfully at %s\n', testCase.Prefix, datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  fprintf('Elapsed time for test was %d minutes and %f seconds\n', floor(elapsedTime/60), rem(elapsedTime,60));
end
