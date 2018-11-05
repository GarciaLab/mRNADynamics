function segmentSpotsMLAndCopyData(testCase, testPath, codePath, dynamicResultsPath, preprocessedDataPath, processedDataPath)
  [tifsPreProcessedData, filterMovieProcessedData, segmentMLDynamicsResults] = createExpectedDataStructureML(testPath, testCase.Prefix);
  
  % Classifier path must end with filesep
  classifiersPath = strcat(codePath, filesep, 'src', filesep, 'classifiers', filesep);
  classifierForTest = ClassifierForTest(classifiersPath, testCase.classifier);

  % Generates Tifs
  filterMovie(testCase.Prefix, 'Tifs');

  % Then copy expected data for filterMovie pass
  copyExpectedDataFolderSpotsML(preprocessedDataPath, tifsPreProcessedData, 'Tifs', testCase.Prefix);

  % Generates expected data for DoGs creation
  filterMovie(testCase.Prefix, 'Weka', classifierForTest, 'ignoreMemoryCheck');

  % Then copy expected data for filterMovie pass
  copyProcessedDataSpotsML(processedDataPath, filterMovieProcessedData, 'FilterMovie', testCase.Prefix);

  % Executes segment spots ML with known DoG
  segmentSpotsML(testCase.Prefix, testCase.Threshold);

  % Then copy expected data for 2nd pass
  % TODO harrypotel: Since now there are two separate functions, we should rename the 1st pass / 2nd pass logic.
  % I don't want to change to much in a single step, so I'll handle that change later on
  copyExpectedDataFolderSpotsML(dynamicResultsPath, segmentMLDynamicsResults, 'SegmentSpotsML', testCase.Prefix);

  disp(['Test data for segment spots completed for Prefix ', testCase.Prefix]);
end

function [tifsPreProcessedData, filterMovieProcessedData, segmentMLDynamicsResults] = createExpectedDataStructureML(testPath, Prefix)
  % Creates root folder if it does not exist
  if 7 ~= exist(testPath, 'dir')
    mkdir(testPath);
  end

  tifsPreProcessedData = createOrCleanExpectedDataSubFolderML(testPath, 'Tifs', 'PreProcessedData', Prefix);
  filterMovieProcessedData = createOrCleanExpectedDataSubFolderML(testPath, 'FilterMovie', 'ProcessedData', Prefix);
  segmentMLDynamicsResults = createOrCleanExpectedDataSubFolderML(testPath, 'SegmentSpotsML', 'DynamicsResults', Prefix);

end

function expectedDataSubFolder = createOrCleanExpectedDataSubFolderML(testPath, step, subfolder, Prefix)
  expectedDataSubFolder = getExpectedDataFolderNameSegmentSpotsML(testPath, step, subfolder, Prefix);
  deleteDirectory(expectedDataSubFolder, Prefix);
  mkdir(expectedDataSubFolder);
end

function copyProcessedDataSpotsML(sourceFolder, expectedDataFolder, step, Prefix)
  disp(['Copying ProcessedData for segment spots ML ', step, ' pass of Prefix ', Prefix]);
  copyfile([sourceFolder, filesep, Prefix, '_', filesep, '*'], expectedDataFolder);
  disp(['Expected data copied to folder ', expectedDataFolder]);
end

function copyExpectedDataFolderSpotsML(sourceFolder, expectedDataFolder, step, Prefix)
  disp(['Copying DynamicResults data for segment spots ML ', step, ' of Prefix ', Prefix]);
  copyfile([sourceFolder, filesep, Prefix, filesep, '*'], expectedDataFolder);
  disp(['Expected data copied to folder ', expectedDataFolder]);
end
