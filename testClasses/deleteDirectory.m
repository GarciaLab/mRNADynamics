function deleteDirectory(dirPath)
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
