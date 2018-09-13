function testCase = testExportDataForFISH(testCase)
  tic;
  disp(['Running ExportDataForFISH test with prefix ', testCase.Prefix]);
  fprintf('Test run started at %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  
  %Figure out the initial folders. 
  CONFIG_CSV_PATH = ['ComputerFolders.csv'];
  
  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');
  
  PreProcPath = getConfigValue(configValues, 'PreProcPath');
  testPath = getConfigValue(configValues, 'TestPath');
  
  %Get file names to compare in preprocessed data folder
  preprocessedDataFolder = strcat(PreProcPath, filesep, testCase.Prefix);
  expectedDataFolder = strcat(testPath, filesep, 'ExportDataForFISH', filesep, testCase.Prefix);

  deleteDirectory(preprocessedDataFolder);

  if (~isprop(testCase, 'PreferredFileName')) 
    ExportDataForFISH(testCase.Prefix, 'keepTifs');
  else 
    ExportDataForFISH(testCase.Prefix, testCase.PreferredFileName, 'keepTifs');
  end

  compareExpectedDataDir(testCase, preprocessedDataFolder, expectedDataFolder);

  elapsedTime = toc;
  fprintf('Test run for %s ended successfully at %s\n', testCase.Prefix, datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
  fprintf('Elapsed time for test was %d minutes and %f seconds\n', floor(elapsedTime/60), rem(elapsedTime,60));
end
