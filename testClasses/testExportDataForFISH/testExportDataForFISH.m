function testCase = testExportDataForFISH(testCase)
  disp(['Running ExportDataForFISH test with prefix ', testCase.Prefix]);
  
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
    ExportDataForFISH(testCase.Prefix);
  else 
    ExportDataForFISH(testCase.Prefix, testCase.PreferredFileName);
  end

  compareExpectedDataDir(testCase, preprocessedDataFolder, expectedDataFolder);
end
