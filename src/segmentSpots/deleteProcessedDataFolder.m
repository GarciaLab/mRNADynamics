function deleteProcessedDataFolder(folder, Prefix)
  if ~contains(folder, 'ProcessedData') || ~contains(folder, Prefix)
    warning(['ProcessedData folder will not be removed. Given folder does not match prefix or ProcessedData folder, ', folder]);
  else
    disp(['Removing ProcessedData folder for Prefix', Prefix]);
    deleteFolder(folder);
  end
end

function deleteFolder(folder)
  rehash;
  try 
    if isequal(exist(folder, 'dir'), 7)
      dirContents = dir(folder);

      for i = 1:length(dirContents)
        fileOrFolder = dirContents(i);
        if isfolder([folder, filesep, fileOrFolder.name]) 
          if ~strcmp('.', fileOrFolder.name) && ~strcmp('..', fileOrFolder.name)
            subfolderPath = [folder, filesep, fileOrFolder.name];
            deleteFolder(subfolderPath);
          end
        else 
          delete([folder, filesep, fileOrFolder.name]);
        end
      end
      
      rehash; 
      [rmdirResult, rmDirErrorMsg] = rmdir(folder, 's');
      if rmdirResult == 0
        warning(['Folder ', folder, ' cannot be removed. Reason: ', rmDirErrorMsg, ' Try closing Matlab and running again.']);
      else
        disp([folder, ' was deleted.']);
      end
    end
  catch ME
    warning(['Folder ', folder, 'cannot be removed. Try closing Matlab and running again.', ME.identifier, ' - ', ME.message]);
  end

end