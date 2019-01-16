function generateExpectedDataTrackmRNADynamicsNotHistone(testCase)
  disp(['Generating test data for trackmRNADynamicsNotHistone test with   Prefix ', testCase.Prefix]);
  
  %Figure out the initial folders. 
  CONFIG_CSV_PATH = ['ComputerFolders.csv'];
  
  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');
  testPath = getConfigValue(configValues, 'TestPath');

  dynamicsResultsPath = getConfigValue(configValues, 'DropboxFolder');
  preProcessedDataPath = getConfigValue(configValues, 'PreProcPath');
  processedDataPath = getConfigValue(configValues, 'FISHPath');

  preProcessedDataExperimentPath = [preProcessedDataPath, filesep, testCase.Prefix];
  dynamicsResultsExperimentPath = [dynamicsResultsPath, filesep, testCase.Prefix];
  processedDataExperimentPath = [processedDataPath, filesep, testCase.Prefix, '_'];

  % Clean up previous runs
  deleteDirectory(dynamicsResultsExperimentPath, testCase.Prefix);
  deleteDirectory(preProcessedDataExperimentPath, testCase.Prefix);
  deleteDirectory(processedDataExperimentPath, testCase.Prefix);

  % Copies required data from SegmentSpotsML step
  copyExpectedDataForPrefix(testCase.Prefix, 'SegmentSpotsML');

  % Creates or cleans existing ExpectedData folder for TrachmRNADynamicsNoHistone
  expectedDataSubFolder = [testPath, filesep, 'TrackmRNADynamicsNoHistone', filesep, 'DynamicsResults', filesep, testCase.Prefix];
  deleteDirectory(expectedDataSubFolder, testCase.Prefix);
  mkdir(expectedDataSubFolder);

  % Executes trackmRNADynamis with known thresholds
  TrackmRNADynamics(testCase.Prefix);
  disp(['Copying expected data for testCase.Prefix ', testCase.Prefix]);
  copyfile([dynamicsResultsExperimentPath, filesep, '*'], expectedDataSubFolder);
  disp(['Expected data copied to folder ', expectedDataSubFolder]);

  % Executes second time to test retracking path  
  expectedDataSubFolder = [testPath, filesep, 'TrackmRNADynamicsNoHistone/secondPass/DynamicsResults', filesep, testCase.Prefix];
  deleteDirectory(expectedDataSubFolder, testCase.Prefix);
  mkdir(expectedDataSubFolder);
  
  TrackmRNADynamics(testCase.Prefix);

  disp(['Copying expected data for testCase.Prefix ', testCase.Prefix]);
  copyfile([dynamicsResultsExperimentPath, filesep, '*'], expectedDataSubFolder);
  disp(['Expected data copied to folder ', expectedDataSubFolder]);
end
