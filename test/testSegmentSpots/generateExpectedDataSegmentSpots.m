% Generates expected data for Segment Spots for the given Prefix in the expected data folder,
% so test cases can compare against it.
% The function runs ExportDataForFISH, filterMovie and segmentSpots and moves the required files to the 
% expected folders.
% If data exists already in the folder, it's deleted beforehand.
function generateExpectedDataSegmentSpots(testCase)
  disp(['Generating expected data for Segment Spot test case with Prefix', testCase.Prefix]);

  %Figure out the initial folders. 
  CONFIG_CSV_PATH = ['ComputerFolders.csv'];
  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');
  testPath = getConfigValue(configValues, 'TestPath');

  if 7 ~= exist(testPath, 'dir')
    mkdir(testPath);
  end

  dynamicResultsPath = getConfigValue(configValues, 'DropboxFolder');
  processedDataPath = getConfigValue(configValues, 'FISHPath');
  
  dynamicResultsExperimentPath = [dynamicResultsPath, filesep, testCase.Prefix];
  processedDataExperimentPath = [processedDataPath, filesep, testCase.Prefix, '_'];

  % Clean up previous runs
  deleteDirectory(dynamicResultsExperimentPath);
  deleteDirectory(processedDataExperimentPath);
  [expectedDynamicsResults1stPass, expectedProcessedData1stPass, expectedDynamicsResults2ndPass, expectedProcessedData2ndPass] = createExpectedDataStructure(testPath, testCase.Prefix);

  % Precondition - Run ExportsDataForFISH without deleting TIFs
  ExportDataForFISH(testCase.Prefix, 'keepTifs');
  
  % Tests first pass
  % Generates DoGs
  filterMovie(testCase.Prefix);
  
  % Then copy expected data for 1st pass 
  copyDynamicResultsData(dynamicResultsPath, expectedDynamicsResults1stPass, '1st', testCase.Prefix);
  copyProcessedData(processedDataPath, expectedProcessedData1stPass, '1st', testCase.Prefix);

  % Executes segment spots with known DoG
  segmentSpots(testCase.Prefix, testCase.DoG);

  % Then copy expected data for 2nd pass 
  % TODO harrypotel: Since now there are two separate functions, we should rename the 1st pass / 2nd pass logic.
  % I don't want to change to much in a single step, so I'll handle that change later on
  copyDynamicResultsData(dynamicResultsPath, expectedDynamicsResults2ndPass, '2nd', testCase.Prefix);
  copyProcessedData(processedDataPath, expectedProcessedData2ndPass, '2nd', testCase.Prefix);

  disp(['Test data for segment spots completed for Prefix ', testCase.Prefix]);
end

function [expectedDynamicsResults1stPass, expectedProcessedData1stPass, expectedDynamicsResults2ndPass, expectedProcessedData2ndPass] = createExpectedDataStructure(testPath, prefix) 
  % Creates root folder if it does not exist
  if 7 ~= exist(testPath, 'dir')
    mkdir(testPath);
  end

  expectedDynamicsResults1stPass = createOrCleanExpectedDataSubFolder(testPath, '1st', 'DynamicsResults', prefix);
  expectedProcessedData1stPass = createOrCleanExpectedDataSubFolder(testPath, '1st', 'ProcessedData', prefix);
  expectedDynamicsResults2ndPass = createOrCleanExpectedDataSubFolder(testPath, '2nd', 'DynamicsResults', prefix);
  expectedProcessedData2ndPass = createOrCleanExpectedDataSubFolder(testPath, '2nd', 'ProcessedData', prefix);
end

function folder = createOrCleanExpectedDataSubFolder(testPath, pass, subfolder, prefix)
  % Adds underscore to the end of PreprocessedData subfolder
  if strcmpi('ProcessedData', subfolder)
    prefix = [prefix, '_'];
  end

  folder = [testPath, filesep, 'SegmentSpots_', pass, 'Pass', filesep, subfolder, filesep, prefix];
  deleteDirectory(folder);
  mkdir folder;
end

function copyDynamicResultsData(sourceFolder, expectedDataFolder, pass, prefix)
  disp(['Copying DynamicResults data for segment spots ', pass, ' pass of Prefix ', prefix]);
  copyfile([sourceFolder, filesep, prefix, filesep, '*'], expectedDataFolder);
  disp(['Expected data copied to folder ', expectedDataFolder]);
end

function copyProcessedData(sourceFolder, expectedDataFolder, pass, prefix)
  disp(['Copying ProcessedData for segment spots ', pass, ' pass of Prefix ', prefix]);
  copyfile([sourceFolder, filesep, prefix, '_',  filesep, '*'], expectedDataFolder);
  disp(['Expected data copied to folder ', expectedDataFolder]);
end
