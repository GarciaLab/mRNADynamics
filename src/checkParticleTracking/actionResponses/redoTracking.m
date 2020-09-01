function [Particles, schnitzcells] = redoTracking(...
    NChannels, cptState, DataFolder, FilePrefix, DropboxFolder)
  
  % [Particles, schnitzcells] = redoTracking(...
  %             NChannels, cptState, DataFolder, FilePrefix, DropboxFolder)
  % 
  % DESCRIPTION
  % Saves changes and calls TrackmRNADynamics with retracking option
  %
  % PARAMETERS
  % NChannels: number of channe;s
  % cptState: object that keeps track of current particle characteristics
  % DataFolder: string of data path
  % FilePrefix: [Prefix, '_'];
  % DropboxFolder: Directory containing pipeline output sets
  %
  % OPTIONS
  % N/A
  %
  % OUTPUT
  % useHistone
  % searchRadiusMicrons
  % retrack
  % displayFigures
  %
  %
  % Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
  % Created: 6/30/2020
  % Last Updated: N/A

  Answer=lower(input('Are you sure you want to redo the tracking?  (y/n) ','s'));
  if Answer=='y'
  %     warning('HG: Not clear that this feature will work with the multiple channels')

      %We need to save the data
      saveChanges(NChannels, cptState, DataFolder, FilePrefix, DropboxFolder);

      disp('Calling TrackmRNADynamics...')
      [Particles,schnitzcells]=TrackmRNADynamics(FilePrefix(1:end-1),'retrack'); % Presumably we want this to be a retrack call

  end
end

