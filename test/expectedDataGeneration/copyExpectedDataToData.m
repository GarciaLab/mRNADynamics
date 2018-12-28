function copyExpectedDataToData(Prefix, dynamicsResultsPath, dynamicsResultsExpectedPath, preprocessedDataPath, preProcessedDataExpectedPath, processedDataPath, processedDataExpectedPath)
  disp(['Copying expected data for Prefix ', Prefix]);

  copyIfFolderIsNotEmpty(dynamicsResultsExpectedPath, dynamicsResultsPath);
  copyIfFolderIsNotEmpty(preProcessedDataExpectedPath, preprocessedDataPath);
  copyIfFolderIsNotEmpty(processedDataExpectedPath, processedDataPath);
end

function copyIfFolderIsNotEmpty(sourcePath, targetPath)
  if exist(sourcePath, 'dir') == 7
    files = dir(sourcePath);
    if length(files) > 2
      copyfile([sourcePath, filesep, '*'], targetPath);
      disp(['Expected data copied to folder ', targetPath]);
    else
      warning(['Folder ', sourcePath, ' is empty. Data is not copied.']);
    end

  else
    warning(['Folder ', sourcePath, ' does not exist or is not a folder. Data is not copied.']);
  end

end

