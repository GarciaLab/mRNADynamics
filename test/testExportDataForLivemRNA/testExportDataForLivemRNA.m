function testCase = testExportDataForLivemRNA(testCase)
  tic;
  disp(['Running ExportDataForLivemRNA test with prefix ', testCase.Prefix]);
  fprintf('Test run started at %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  
  %Figure out the initial folders. 
  CONFIG_CSV_PATH = ['ComputerFolders.csv'];
  
  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');
  
  dynamicsResultsPath = getConfigValue(configValues, 'DropboxFolder');
  PreProcPath = getConfigValue(configValues, 'PreProcPath');
  testPath = getConfigValue(configValues, 'TestPath');
  
  %Get file names to compare in preprocessed data folder
  preprocessedDataFolder = strcat(PreProcPath, filesep, testCase.Prefix);
  expectedPreProcFolder = strcat(testPath, filesep, 'ExportDataForLivemRNA', filesep, 'PreProcessedData', filesep,...
    testCase.Prefix);
  
  dynamicResultsPath = strcat(dynamicsResultsPath, filesep, testCase.Prefix);
  
  deleteDirectory(preprocessedDataFolder, testCase.Prefix);
  deleteDirectory(dynamicResultsPath, testCase.Prefix);

  if (~isprop(testCase, 'PreferredFileName')) 
    ExportDataForLivemRNA(testCase.Prefix, 'keepTifs');
  else 
    ExportDataForLivemRNA(testCase.Prefix, testCase.PreferredFileName, 'keepTifs');
  end

  assertFrameInfoEqualToExpected(testCase, dynamicResultsPath, testPath, 'ExportDataForLivemRNA');
  compareExpectedDataDir(testCase, preprocessedDataFolder, expectedPreProcFolder);

  elapsedTime = toc;
  fprintf('Test run for %s ended successfully at %s\n', testCase.Prefix, datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  fprintf('Elapsed time for test was %d minutes and %f seconds\n', floor(elapsedTime/60), rem(elapsedTime,60));
end
