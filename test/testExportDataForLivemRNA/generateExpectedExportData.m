% Generates expected data for ExportDataForLivemRNA for the given Prefix in the expected data folder,
% so test cases can compare against it.
% The function runs ExportDataForLivemRNA and moves the required files to the expected folder.
% If data exists already in the folder, it's deleted beforehand.
function generateExpectedExportData(testCase)
  disp(['Generating expected data for ExportDataForLivemRNA test case with Prefix', testCase.Prefix]);

  %Figure out the initial folders.
  CONFIG_CSV_PATH = ['ComputerFolders.csv'];

  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');

  dynamicResultsPath = getConfigValue(configValues, 'DropboxFolder');
  PreProcPath = getConfigValue(configValues, 'PreProcPath');
  testPath = getConfigValue(configValues, 'TestPath');
  
  if 7 ~= exist(testPath, 'dir')
    mkdir(testPath);
  end

  %Get file names to compare in preprocessed data folder
  
  preprocessedDataPath = strcat(PreProcPath, filesep, testCase.Prefix);
  dynamicResultsPath = strcat(dynamicResultsPath, filesep, testCase.Prefix);

  deleteDirectory(preprocessedDataPath, testCase.Prefix);
  deleteDirectory(dynamicResultsPath, testCase.Prefix);

  exportAndCopyData(testCase, testPath, dynamicResultsPath, preprocessedDataPath);

end
