% Added PreferredFileName so we can automate testing and bypass the user prompt when there are many files available.
function FrameInfo = processLIFExportMode(rawDataFolder, ExperimentType, ProjectionType, Channel1, Channel2, Channel3,...
    Prefix, OutputFolder, PreferredFileNameForTest, keepTifs, nuclearGUI, skipExtraction)
  
  %Loads file and metadata
  [XMLFolder, seriesPropertiesXML, seriesXML] = getSeriesFiles(rawDataFolder);
  
  [~, ~, LIFImages, LIFMeta] = loadLIFFile(rawDataFolder);
  
  %Obtains frames information
  [NSeries, NFrames, NSlices, NPlanes, NChannels, Frame_Times] = getFrames(LIFMeta);
  
  if sum(NFrames)~=0
    [Frame_Times, First_Time] = obtainFrameTimes(XMLFolder, seriesPropertiesXML, NSeries, NFrames, NSlices, NChannels);
    [InitialStackTime, zPosition] = getFirstSliceTimestamp(NSlices, NSeries, NPlanes, NChannels, Frame_Times, XMLFolder, seriesXML);
  else
      InitialStackTime = [];
      zPosition = [];
  end
  FrameInfo = recordFrameInfo(NFrames, NSlices, InitialStackTime, LIFMeta, zPosition);
  
  %Find the flat field (FF) information
  LIFExportMode_flatFieldImage(LIFMeta, rawDataFolder, OutputFolder, Prefix, PreferredFileNameForTest);
  
  [coatChannel, histoneChannel, fiducialChannel, inputProteinChannel, FrameInfo] =...
    LIFExportMode_interpretChannels(ExperimentType, Channel1, Channel2, Channel3, FrameInfo);

  if sum(NFrames) == 0
      NFrames = ~NFrames;
  end
  
  if ~skipExtraction
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
          processLIFFrame(numberOfFrames, Prefix, BlankImage, OutputFolder, LIFImages, framesIndex, seriesIndex,...
            NChannels, NSlices, ExperimentType, Channel1, Channel2, Channel3, ProjectionType, fiducialChannel,...
            histoneChannel, ReferenceHist, coatChannel, inputProteinChannel, false); %JP: hardcode false zPadding for now
          numberOfFrames = numberOfFrames + 1;
        end
      end
  end
  
  close(waitbarFigure)
  
end
