% Generates expected data for ExportDataForFISH for the given Prefix in the expected data folder,
% so test cases can compare against it.
% The function runs ExportDataForFISH and moves the required files to the expected folder.
% If data exists already in the folder, it's deleted beforehand.
function generateExpectedExportData(testCase)
  disp(['Generating expected data for ExportDataForFISH test case with Prefix', testCase.Prefix]);

  %Figure out the initial folders.
  CONFIG_CSV_PATH = ['ComputerFolders.csv'];

  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');

  dynamicsResultsPath = getConfigValue(configValues, 'DropboxFolder');
  PreProcPath = getConfigValue(configValues, 'PreProcPath');
  testPath = getConfigValue(configValues, 'TestPath');
  
  if 7 ~= exist(testPath, 'dir')
    mkdir(testPath);
  end

  %Get file names to compare in preprocessed data folder
  experimentTestRootPath = strcat(testPath, filesep, 'ExportDataForFISH', filesep);

  preprocessedDataFolder = strcat(PreProcPath, filesep, testCase.Prefix);
  expectedPreProcFolder = strcat(experimentTestRootPath, 'PreProcessedData', filesep, testCase.Prefix);
  
  dynamicsResultsDataFolder = strcat(dynamicsResultsPath, filesep, testCase.Prefix);
  expectedDynamicsResultsFolder = strcat(experimentTestRootPath, 'DynamicsResults', filesep, testCase.Prefix);

  deleteDirectory(preprocessedDataFolder, testCase.Prefix);
  deleteDirectory(expectedPreProcFolder, testCase.Prefix);
  
  deleteDirectory(dynamicsResultsDataFolder, testCase.Prefix);
  deleteDirectory(expectedDynamicsResultsFolder, testCase.Prefix);

  mkdir(expectedPreProcFolder);
  mkdir(expectedDynamicsResultsFolder);

  if (~isprop(testCase, 'PreferredFileName'))
    ExportDataForFISH(testCase.Prefix, 'keepTifs');
  else
    testCase.initializeTestCase;
    ExportDataForFISH(testCase.Prefix, testCase.PreferredFileName, 'keepTifs');
  end

  disp(['Copying expected data for Prefix ', testCase.Prefix]);
  copyfile([dynamicsResultsDataFolder, filesep, '*'], expectedDynamicsResultsFolder);
  disp(['Expected data copied to folder ', expectedDynamicsResultsFolder]);
  copyfile([preprocessedDataFolder, filesep, '*'], expectedPreProcFolder);
  disp(['Expected data copied to folder ', expectedPreProcFolder]);

end
