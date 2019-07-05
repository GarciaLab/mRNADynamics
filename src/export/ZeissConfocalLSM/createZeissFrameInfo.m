function [FrameRange, FrameInfo] = createZeissFrameInfo(LSMIndex, NFrames, NSlices, FrameInfo, LSMMeta, Frame_Times, ValueField, zslicesPadding)
  %Save the information in FrameInfo
  if LSMIndex == 1
    FrameRange = 1:NFrames(LSMIndex);
  else
    FrameRange = (1:NFrames(LSMIndex)) + length(FrameInfo);
  end
  
  % if zPadding was indicated in the arguments, we round up to the series
  % with more z-slices (because we'll pad with blank images the other series)
  if (zslicesPadding)
    topZSlice = max(NSlices);
  else
    % if no zPadding, we round down to the series with less z-slices
    topZSlice = min(NSlices);
  end

  for i = FrameRange
    FrameInfo(i).LinesPerFrame = str2double(LSMMeta.getPixelsSizeY(0));
    FrameInfo(i).PixelsPerLine = str2double(LSMMeta.getPixelsSizeX(0));
    FrameInfo(i).NumberSlices = topZSlice;
    FrameInfo(i).FileMode = 'LSMExport';

    if ValueField
      FrameInfo(i).PixelSize = str2num(LSMMeta.getPixelsPhysicalSizeX(0).value);
      FrameInfo(i).ZStep = str2double(LSMMeta.getPixelsPhysicalSizeZ(0).value);
    else
      FrameInfo(i).PixelSize = str2num(LSMMeta.getPixelsPhysicalSizeX(0));
      FrameInfo(i).ZStep = str2double(LSMMeta.getPixelsPhysicalSizeZ(0));
    end

    FrameInfo(i).Time = Frame_Times(i); % In seconds
  end
end
