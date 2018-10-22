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

  % Classifier path must end with filesep
  classifiersPath = strcat(codePath, filesep, 'src', filesep, 'classifiers', filesep);
  classifierForTest = ClassifierForTest(classifiersPath, testCase.classifier);

  % Clean up previous runs
  deleteDirectory(dynamicResultsExperimentPath, testCase.Prefix);
  deleteDirectory(preprocessedDataExperimentPath, testCase.Prefix);
  deleteDirectory(processedDataExperimentPath, testCase.Prefix);
  
  [tifsPreProcessedData, filterMovieProcessedData, segmentMLDynamicsResults, ...
    segmentMLProcessedData] = createExpectedDataStructureML(testPath, testCase.Prefix);

  % Precondition - Run ExportsDataForFISH without deleting TIFs
  ExportDataForFISH(testCase.Prefix, 'keepTifs');

  % Generates Tifs 
  filterMovie(testCase.Prefix, 'Tifs');

  % Then copy expected data for filterMovie pass 
  copyExpectedDataFolder(preprocessedDataPath, tifsPreProcessedData, 'Tifs', testCase.Prefix);

  % Generates expected data for DoGs creation
  filterMovie(testCase.Prefix, classifierForTest, 'ignoreMemoryCheck');
  
  % Then copy expected data for filterMovie pass 
  copyProcessedDataML(processedDataPath, filterMovieProcessedData, 'FilterMovie', testCase.Prefix);

  % Executes segment spots ML with known DoG
  segmentSpotsML(testCase.Prefix, testCase.Threshold);

  % Then copy expected data for 2nd pass 
  % TODO harrypotel: Since now there are two separate functions, we should rename the 1st pass / 2nd pass logic.
  % I don't want to change to much in a single step, so I'll handle that change later on
  copyExpectedDataFolder(dynamicResultsPath, segmentMLDynamicsResults, 'SegmentSpotsML', testCase.Prefix);
  copyProcessedDataML(processedDataPath, segmentMLProcessedData, 'SegmentSpotsML', testCase.Prefix);

  disp(['Test data for segment spots completed for Prefix ', testCase.Prefix]);
end

function [tifsPreProcessedData, filterMovieProcessedData, segmentMLDynamicsResults, segmentMLProcessedData] = createExpectedDataStructureML(testPath, prefix) 
  % Creates root folder if it does not exist
  if 7 ~= exist(testPath, 'dir')
    mkdir(testPath);
  end

  tifsPreProcessedData = createOrCleanExpectedDataSubFolderML(testPath, 'Tifs', 'PreProcessedData', prefix);
  filterMovieProcessedData = createOrCleanExpectedDataSubFolderML(testPath, 'FilterMovie', 'ProcessedData', prefix);
  segmentMLDynamicsResults = createOrCleanExpectedDataSubFolderML(testPath, 'SegmentSpotsML', 'DynamicsResults', prefix);
  segmentMLProcessedData = createOrCleanExpectedDataSubFolderML(testPath, 'SegmentSpotsML', 'ProcessedData', prefix);

end

function expectedDataSubFolder = createOrCleanExpectedDataSubFolderML(testPath, step, subfolder, prefix)
  % Adds underscore to the end of PreprocessedData subfolder
  if strcmpi('ProcessedData', subfolder)
    prefix = [prefix, '_'];
  end

  expectedDataSubFolder = [testPath, filesep, 'SegmentSpotsML', filesep, step, filesep, subfolder, filesep, prefix];
  deleteDirectory(expectedDataSubFolder, prefix);
  mkdir( expectedDataSubFolder);
end

function copyExpectedDataFolder(sourceFolder, expectedDataFolder, step, prefix)
  disp(['Copying DynamicResults data for segment spots ML ', step, ' of Prefix ', prefix]);
  copyfile([sourceFolder, filesep, prefix, filesep, '*'], expectedDataFolder);
  disp(['Expected data copied to folder ', expectedDataFolder]);
end

function copyProcessedDataML(sourceFolder, expectedDataFolder, step, prefix)
  disp(['Copying ProcessedData for segment spots ML ', step, ' pass of Prefix ', prefix]);
  copyfile([sourceFolder, filesep, prefix, '_',  filesep, '*'], expectedDataFolder);
  disp(['Expected data copied to folder ', expectedDataFolder]);
end
