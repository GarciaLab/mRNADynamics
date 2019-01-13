% Added PreferredFileName so we can automate testing and bypass the user prompt when there are many files available.
function FrameInfo = processLIFExportMode(Folder, ExperimentType, ProjectionType, Channel1, Channel2, Channel3, Prefix, OutputFolder, PreferredFileNameForTest, keepTifs, nuclearGUI)
  
  %Loads file and metadata
  [XMLFolder, seriesPropertiesXML, seriesXML] = getSeriesFiles(Folder);
  
  [~, ~, LIFImages, LIFMeta] = loadLIFFile(Folder);
  
  %Obtains frames information
  [NSeries, NFrames, NSlices, NPlanes, NChannels, Frame_Times] = getFrames(LIFMeta);
  
  [Frame_Times, First_Time] = obtainFrameTimes(XMLFolder, seriesPropertiesXML, NSeries, NFrames, NSlices, NChannels);
  
  [InitialStackTime, zPosition] = getFirstSliceTimestamp(NSlices, NSeries, NPlanes, NChannels, Frame_Times, XMLFolder, seriesXML);
  
  FrameInfo = recordFrameInfo(NFrames, NSlices, InitialStackTime, LIFMeta, zPosition);
 
  %Find the flat field (FF) information
  LIFExportMode_flatFieldImage(LIFMeta, Folder, OutputFolder, Prefix, PreferredFileNameForTest);
  
  [coatChannel, histoneChannel, fiducialChannel, inputProteinChannel, FrameInfo] =...
    LIFExportMode_interpretChannels(ExperimentType, Channel1, Channel2, Channel3, FrameInfo);

  %Copy the data
  waitbarFigure = waitbar(0, 'Extracting LIFExport images');
  %Create a blank image
  BlankImage = uint16(zeros(size(LIFImages{1}{1,1})));
  
  %Counter for number of frames
  numberOfFrames = 1;        
  %Load the reference histogram for the fake histone channel
  load('ReferenceHist.mat')
  if nuclearGUI
      [Channel1, Channel2, Channel3, ProjectionType] = chooseNuclearChannels(...
        LIFImages, NSeries, NSlices, NChannels, NFrames, ProjectionType, Channel1, Channel2, ...
        Channel3, ReferenceHist);
  end
  for seriesIndex = 1:NSeries
    waitbar(seriesIndex/NSeries, waitbarFigure)
    for framesIndex = 1:NFrames(seriesIndex) 
      processLIFFrame(numberOfFrames, Prefix, BlankImage, OutputFolder, LIFImages, framesIndex, seriesIndex, NChannels, NSlices, ExperimentType, Channel1, Channel2, Channel3, ProjectionType, fiducialChannel, histoneChannel, ReferenceHist, coatChannel, inputProteinChannel);
      numberOfFrames = numberOfFrames + 1;
    end
  end

  if ~keepTifs
    removeUnwantedTIFs(Folder);
  end

  close(waitbarFigure)
end

% Removes all TIF files from the original folder in RawDynamicsData
function removeUnwantedTIFs(Folder) 
  cd(Folder);
  tifs = dir('*.tif');

  allTifs = {tifs.name};

  if numel(allTifs) > 1
    disp(['Removing TIF files from source folder ', Folder]);
    for i = 2:numel(allTifs)
%       disp(['Deleting ', allTifs{i}]);
      delete(allTifs{i});
    end
  end

end
