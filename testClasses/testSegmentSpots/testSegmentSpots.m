function testCase = testSegmentSpots(testCase)
  disp(['Running with segment spot test with prefix ', testCase.Prefix, ' and DoG ', testCase.DoG]);
  
  %Figure out the initial folders. 
  CONFIG_CSV_PATH = ['ComputerFolders.csv'];
  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');
  testPath = getConfigValue(configValues, 'TestPath');
  
  dynamicResultsPath = getConfigValue(configValues, 'DropboxFolder');
  processedDataPath = getConfigValue(configValues, 'FISHPath');
  
  dynamicResultsExperimentPath = [dynamicResultsPath, filesep, testCase.Prefix];
  processedDataExperimentPath = [processedDataPath, filesep, testCase.Prefix];

  deleteDirectory(dynamicResultsExperimentPath);
  deleteDirectory(processedDataExperimentPath);

  % Precondition - Run ExportsDataForFISH
  ExportDataForFISH(testCase.Prefix);
  
  % Tests first pass
  % Executes segment spots without specifying DoG
  segmentSpots(testCase.Prefix, []);

  % Then compares expected contents of both full DynamicsResults and ProcessedData folders
  expectedDataFolder = [testPath, filesep, 'SegmentSpots_1stPass', filesep, 'DynamicsResults',...
    filesep, testCase.Prefix];
  compareExpectedDataDir(testCase, dynamicResultsExperimentPath, expectedDataFolder);

  expectedDataFolder = [testPath, filesep, 'SegmentSpots_1stPass', filesep, 'ProcessedData',...
    filesep, testCase.Prefix];
  compareExpectedDataDir(testCase, processedDataExperimentPath, expectedDataFolder);

  % Tests second pass
  % Executes segment spots with known DoG
  segmentSpots(testCase.Prefix, testCase.DoG);

  % Then compares expected contents of both full DynamicsResults and ProcessedData folders
  expectedDataFolder = [testPath, filesep, 'SegmentSpots_2ndPass', filesep, 'DynamicsResults',...
    filesep, testCase.Prefix];
  compareExpectedDataDir(testCase, dynamicResultsExperimentPath, expectedDataFolder);

  expectedDataFolder = [testPath, filesep, 'SegmentSpots_2ndPass', filesep, 'ProcessedData',...
    filesep, testCase.Prefix];
  compareExpectedDataDir(testCase, processedDataExperimentPath, expectedDataFolder);
end
