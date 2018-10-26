function deleteDirectory(dirPath, expectedSubpath)
  validateDirectory(dirPath, expectedSubpath);

  try 
    if isequal(exist(dirPath, 'dir'), 7)
      dirContents = dir(dirPath);

      for i = 1:length(dirContents)
        fileOrFolder = dirContents(i);
        if isfolder([dirPath, filesep, fileOrFolder.name]) 
          if ~strcmp('.', fileOrFolder.name) && ~strcmp('..', fileOrFolder.name)
            deleteDirectory([dirPath, filesep, fileOrFolder.name]);
          end
        else 
          delete([dirPath, filesep, fileOrFolder.name]);
        end
      end
      
      rehash; 
      rmdir(dirPath);
    end
  catch
    warning('Directory cannot be removed. Try closing Matlab and running again. Proceeding with test.');
  end

end

% Validates that the directory to be deleted is a subfolder of LivemRNA and also a subfolder of an expected subpath.
% These validations are included so nobody can accidentally delete folders in the drive just because ComputerFolders
% file is not properly configured.
% In case that the validations fail, it'll prompt the user about the location and proceed if the user confirms.
function validateDirectory(dirPath, expectedSubpath) 
  jenkinsFolder = 'D:\Data\JuanPabloPicasso\Data';
  directoryOkToDelete = true;  

  % Ignores validation if the dirPath is the Jenkins folder
  if ~contains(dirPath, jenkinsFolder) 
    % Validates that the dirPath contains the 'LivemRNA' string in it
    if ~contains(dirPath, 'LivemRNA')
      directoryOkToDelete = false;
    else
      subPath = extractAfter(dirPath, 'LivemRNA');
      
      % Validates that the remaining path after LivemRNA is not empty, meaning that the whole folder will be deleted
      if length(subPath) < 2
        directoryOkToDelete = false;
      end

      if ~contains(subPath, expectedSubpath)
        directoryOkToDelete = false;
      end

    end

  end

  if ~directoryOkToDelete
    buttonYesText = 'Yes, proceed with directory deletion.';
    buttonCancelText = 'No, cancel process.';

    question = ['You''re about to delete a suspicious directory, please confirm that you''d like to delete ', dirPath];  

    answer = questdlg(question, 'Directory delete confirmation', buttonYesText, buttonCancelText, buttonCancelText);
    if strcmp(answer, buttonCancelText)
      ME = MException('UserAborted:FolderDeletionCancelled', ['User confirmed folder ', dirPath, ' should not be deleted']);
      throw(ME);
    end
    
  end

end
