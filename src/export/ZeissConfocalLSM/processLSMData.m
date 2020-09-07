% NL: Updated this drawing heavily from processLIFExportMode
% Last updated: 2020-09-03

function FrameInfo = processLSMData(Folder, D, FrameInfo,...
    Channels, ProjectionType,Prefix, OutputFolder,nuclearGUI,...
    skipExtraction)
    
    disp('Exporting movie file...');
    
    cleanupObj = onCleanup(@myCleanupFun);

    moviePrecision = 'uint16';
    hisPrecision = 'uint16';
    
    %Load the reference histogram for the fake histone channel
    load('ReferenceHist.mat', 'ReferenceHist');    
    
    % get basic info
    NSeries = length(D);      
    
    % This chunk makes FrameInfo                     
    Frame_Times = []; % Store the frame information
    AllLSMImages = cell(1, NSeries);
    for LSMIndex = 1:NSeries        

      %Load the file using BioFormats
      LSMImages = bfopen([Folder, filesep, D(LSMIndex).name]);
      % Extract the metadata for each series
      LSMMeta = LSMImages{:, 4}; % OME Metadata
      LSMMeta2 = LSMImages{:, 2}; % Original Metadata
      AllLSMImages{LSMIndex} = LSMImages{1};

      % Figure out the number of slices in each series
      NSlices(LSMIndex) = str2double(LSMMeta.getPixelsSizeZ(0));
      % Number of channels
      NChannels(LSMIndex) = LSMMeta.getChannelCount(0);
      % Total number of planes acquired
      NPlanes(LSMIndex) = LSMMeta.getPlaneCount(0);
      % Finally, use this information to determine the number of frames in
      % each series
      NFrames(LSMIndex) = NPlanes(LSMIndex) / NSlices(LSMIndex) / NChannels(LSMIndex);

      % Check that the acquisition wasn't stopped before the end of a
      % cycle. If that is the case, some images in the last frame will
      % be blank and we need to remove them.
      if sum(sum(LSMImages{1}{end,1})) == 0
        % Reduce the number of frames by one
        NFrames(LSMIndex) = NFrames(LSMIndex) - 1;
        % Reduce the number of planes by NChannels*NSlices
        NPlanes(LSMIndex) = NPlanes(LSMIndex) - NChannels(LSMIndex) * NSlices(LSMIndex);
      end

      NDigits = getNDigits(NFrames, LSMIndex);

      try
          StartingTime(LSMIndex) = obtainZeissStartingTime(Folder, LSMIndex, LSMMeta2, NDigits);
      catch
          StartingTime(LSMIndex) = obtainZeissStartingTime(Folder, LSMIndex, LSMMeta2, NDigits+1);
      end       
    end
    [ValueField, Frame_Times] = obtainZeissFrameTimes(LSMMeta, NSlices, LSMIndex, NPlanes, NChannels, StartingTime, Frame_Times);
    [~, FrameInfo] = createZeissFrameInfo(LSMIndex, NFrames, NSlices, FrameInfo, LSMMeta, Frame_Times, ValueField);
    
    if nuclearGUI
      warning('NL: GUI option for LSM not implemented. Going with default His export options')
    end
    
    if ~skipExtraction
                                       
      exportTifStacks(FrameInfo, AllLSMImages, 'LSM', min(NChannels), Prefix, ...
          moviePrecision, hisPrecision, nuclearGUI, ProjectionType, Channels, ReferenceHist)
      

      try close(waitbarFigure); catch; end

      [FFPaths, FFToUse, LSMFF] = findFlatFieldInformation(Folder);
      processFlatFieldInformation(Prefix, OutputFolder, FFPaths, FFToUse, LSMFF);
    
    end
end

