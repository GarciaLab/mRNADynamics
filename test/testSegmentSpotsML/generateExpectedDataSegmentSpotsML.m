% Generates expected data for Segment Spots ML for the given Prefix in the expected data folder,
% so test cases can compare against it.
% The function runs ExportDataForFISH, Tifs generation, filterMovie and segmentSpotsML and moves 
% the required files to the expected folders.
% If data exists already in the folder, it's deleted beforehand.
function generateExpectedDataSegmentSpotsML(testCase)
  disp(['Generating expected data for Segment Spot ML test case with Prefix', testCase.Prefix]);

  %Figure out the initial folders. 
  CONFIG_CSV_PATH = ['ComputerFolders.csv'];
  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');
  
  codePath = getConfigValue(configValues, 'MS2CodePath');
  testPath = getConfigValue(configValues, 'TestPath');

  if 7 ~= exist(testPath, 'dir')
    mkdir(testPath);
  end

  dynamicResultsPath = getConfigValue(configValues, 'DropboxFolder');
  preprocessedDataPath = getConfigValue(configValues, 'PreProcPath');
  processedDataPath = getConfigValue(configValues, 'FISHPath');

  preprocessedDataExperimentPath = [preprocessedDataPath, filesep, testCase.Prefix];
  dynamicResultsExperimentPath = [dynamicResultsPath, filesep, testCase.Prefix];
  processedDataExperimentPath = [processedDataPath, filesep, testCase.Prefix, '_'];
  
  deleteDirectory(dynamicResultsExperimentPath, testCase.Prefix);
  deleteDirectory(preprocessedDataExperimentPath, testCase.Prefix);
  deleteDirectory(processedDataExperimentPath, testCase.Prefix);

  % Precondition - Run ExportsDataForFISH without deleting TIFs
  ExportDataForFISH(testCase.Prefix, 'keepTifs');

  segmentSpotsMLAndCopyData(testCase, testPath, codePath, dynamicResultsPath, preprocessedDataPath, processedDataPath);
end
