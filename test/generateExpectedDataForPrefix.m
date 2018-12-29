% Generates expected data for Segment Spots ML for the given Prefix in the expected data folder,
% so test cases can compare against it.
% The function runs ExportDataForLivemRNA, Tifs generation, filterMovie and segmentSpotsML and moves
% the required files to the expected folders.
% If data exists already in the folder, it's deleted beforehand.
function generateExpectedDataForPrefix(exportTestCase, segmentSpotsTestCase, trackmRNADynamicsTestCase)
  Prefix = exportTestCase.Prefix;
  disp(['Generating expected data for test case with Prefix', Prefix]);

  %Figure out the initial folders.
  CONFIG_CSV_PATH = ['ComputerFolders.csv'];
  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');

  codePath = getConfigValue(configValues, 'MS2CodePath');
  testPath = getConfigValue(configValues, 'TestPath');

  if 7 ~= exist(testPath, 'dir')
    mkdir(testPath);
  end

  dynamicsResultsPath = getConfigValue(configValues, 'DropboxFolder');
  preProcessedDataPath = getConfigValue(configValues, 'PreProcPath');
  processedDataPath = getConfigValue(configValues, 'FISHPath');

  preProcessedDataExperimentPath = [preProcessedDataPath, filesep, Prefix];
  dynamicsResultsExperimentPath = [dynamicsResultsPath, filesep, Prefix];
  processedDataExperimentPath = [processedDataPath, filesep, Prefix, '_'];

  % Clean up previous runs
  deleteDirectory(dynamicsResultsExperimentPath, Prefix);
  deleteDirectory(preProcessedDataExperimentPath, Prefix);
  deleteDirectory(processedDataExperimentPath, Prefix);

  % Export Data
  exportAndCopyData(exportTestCase, testPath, dynamicsResultsExperimentPath, preProcessedDataExperimentPath);

  % Filter movie to generate Tifs
  filterMovieTifsAndCopyData(segmentSpotsTestCase, testPath, dynamicsResultsExperimentPath, preProcessedDataExperimentPath);

  % Filter movie to generate DoGs with Weka
  filterMovieWekaAndCopyData(segmentSpotsTestCase, testPath, codePath, processedDataExperimentPath);

  % Segment Spots ML
  segmentSpotsMLAndCopyData(segmentSpotsTestCase, testPath, dynamicsResultsExperimentPath, preProcessedDataExperimentPath, processedDataExperimentPath);

  % TrackNuclei
  trackNucleiAndCopyData(Prefix, testPath, dynamicsResultsExperimentPath, preProcessedDataExperimentPath, processedDataExperimentPath);

  % TrackmRNADynamics
  trackmRNADynamicsAndCopyData(trackmRNADynamicsTestCase, testPath, dynamicsResultsExperimentPath, preProcessedDataExperimentPath,...
    processedDataExperimentPath);
end
