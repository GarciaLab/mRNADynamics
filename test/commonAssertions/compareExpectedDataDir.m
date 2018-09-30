% Compares two folders, one with expected data and the other with real values after execution.
% Takes subfolders into account.
function testCase = compareExpectedDataDir(testCase, dataFolder, expectedDataFolder)
  %Get file names to compare in preprocessed data folder
  disp(dataFolder);
  cd(dataFolder);
  sourceFiles = dir;

  numberOfFiles = length(sourceFiles) - 2;
  filesToCompare = {numberOfFiles};
  foldersToCompare = [];
  filesToCompareIndex = 1;
  foldersToCompareIndex = 1;

  %Traverse the files in folder, but ignore the first two as they are . and ..
  for i = 3:length(sourceFiles)
    fileOrFolder = sourceFiles(i).name;

    if ~ isfolder([dataFolder, filesep, fileOrFolder])
      filesToCompare{filesToCompareIndex} = fileOrFolder;
      filesToCompareIndex = filesToCompareIndex + 1;
    else
      foldersToCompare{foldersToCompareIndex} = fileOrFolder;
      foldersToCompareIndex = foldersToCompareIndex + 1;
    end

  end

  %Run comparison
  expectedCompareResult = 'FC: no differences encountered';

  % Compares all files in the folder
  for i = 1:length(filesToCompare)
    sourceFile = strcat(dataFolder, filesep, filesToCompare{i});
    targetFile = strcat(expectedDataFolder, filesep, filesToCompare{i});
    compareCommand = strcat("fc ", sourceFile, " ", targetFile);
    [status, cmdout] = system(compareCommand, '-echo');
    disp(status);
    testCase.assertFalse(~ contains(cmdout, expectedCompareResult))
  end

  % Now compare folders recursively
  if ~ isempty(foldersToCompare)

    for i = 1:length(foldersToCompare)
      dataSubFolder = strcat(dataFolder, filesep, foldersToCompare{i});
      expectedDataSubFolder = strcat(expectedDataFolder, filesep, foldersToCompare{i});
      compareExpectedDataDir(testCase, dataSubFolder, expectedDataSubFolder);
    end

  end

end
