% Generates expected data for ExportDataForFISH for the given Prefix in the expected data folder,
% so test cases can compare against it.
% The function runs ExportDataForFISH and moves the required files to the expected folder.
% If data exists already in the folder, it's deleted beforehand.
function generateExpectedData_ExportDataForFISH(testCase)
  disp('Generating expected data for ExportDataForFISH test cases');

  %Figure out the initial folders.
  CONFIG_CSV_PATH = ['ComputerFolders.csv'];

  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');

  PreProcPath = getConfigValue(configValues, 'PreProcPath');
  testPath = getConfigValue(configValues, 'TestPath');
  
  if 7 ~= exist(testPath, 'dir')
    mkdir(testPath);
  end

  %Get file names to compare in preprocessed data folder
  preprocessedDataFolder = strcat(PreProcPath, filesep, testCase.Prefix);
  expectedDataFolder = strcat(testPath, filesep, 'ExportDataForFISH', filesep, testCase.Prefix);

  deleteDirectory(preprocessedDataFolder);
  deleteDirectory(expectedDataFolder);
  mkdir(expectedDataFolder);

  if (~isprop(testCase, 'PreferredFileName'))
    ExportDataForFISH(testCase.Prefix, 'keepTifs');
  else
    testCase.initializeTestCase;
    ExportDataForFISH(testCase.Prefix, testCase.PreferredFileName, 'keepTifs');
  end

  copyfile([preprocessedDataFolder, filesep, '*'], expectedDataFolder);

end
