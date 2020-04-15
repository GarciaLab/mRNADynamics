% Removes all TIF files from the original folder in RawDynamicsData
function removeUnwantedTIFs(rawDataFolder) 

  oldFolder = cd(rawDataFolder);
  
  tifs = dir('*.tif');
  jpgs = dir('*.jpg');
  tifs = [tifs, jpgs];

  allTifs = {tifs.name};

  if numel(allTifs) > 1
    disp(['Removing TIF files from source folder ', rawDataFolder]);
    for i = 2:numel(allTifs)
      delete(allTifs{i});
    end
  end
  
  cd(oldFolder);

end
