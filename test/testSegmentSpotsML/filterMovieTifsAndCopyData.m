function filterMovieTifsAndCopyData(testCase, testPath, dynamicsResultsPath, preProcessedDataPath)
  expectedDataFolder = [testPath, '/filterMovieTifs'];
  if 7 ~= exist(expectedDataFolder, 'dir')
    mkdir(expectedDataFolder);
  end

  dynamicsResultsExpectedPath = [expectedDataFolder, '/DynamicsResults/', testCase.Prefix];
  preProcessedDataExpectedPath = [expectedDataFolder, '/PreProcessedData/', testCase.Prefix];
  
  % Generates Tifs
  filterMovie(testCase.Prefix, 'Tifs');

  % Then copy expected data for filterMovie pass
  copyExpectedDataFolder(dynamicsResultsPath, dynamicsResultsExpectedPath, 'filterMovie Tifs', testCase.Prefix)
  copyExpectedDataFolder(preProcessedDataPath, preProcessedDataExpectedPath, 'filterMovie Tifs', testCase.Prefix)

  disp(['Test data for filterMovie to generate Tifs completed for Prefix ', testCase.Prefix]);
end
