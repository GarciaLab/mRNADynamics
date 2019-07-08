% Removes all TIF files from the original folder in RawDynamicsData
function removeUnwantedTIFs(rawDataFolder) 

  oldFolder = cd(rawDataFolder);
  
  tifs = dir('*.tif');

  allTifs = {tifs.name};

  if numel(allTifs) > 1
    disp(['Removing TIF files from source folder ', rawDataFolder]);
    parfor i = 2:numel(allTifs)
      delete(allTifs{i});
    end
  end
  
  cd(oldFolder);

end
