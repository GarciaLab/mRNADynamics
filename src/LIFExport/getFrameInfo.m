function FrameInfo = getFrameInfo(NFrames, NSlices, InitialStackTime, LIFMeta, zGalvo)
  for i = 1:sum(NFrames)
    FrameInfo(i).LinesPerFrame = str2double(LIFMeta.getPixelsSizeY(0));
    FrameInfo(i).PixelsPerLine = str2double(LIFMeta.getPixelsSizeX(0));
    FrameInfo(i).NumberSlices = min(NSlices);
    
    FrameInfo(i).FileMode = 'LIFExport';

    %This is to allow for backwards compatibility with BioFormats
    if ~isempty(str2num(LIFMeta.getPixelsPhysicalSizeX(0)))
        FrameInfo(i).PixelSize = str2num(LIFMeta.getPixelsPhysicalSizeX(0));
        FrameInfo(i).ZStep = str2double(LIFMeta.getPixelsPhysicalSizeZ(0));
    else
        FrameInfo(i).PixelSize = str2num(LIFMeta.getPixelsPhysicalSizeX(0).value);
        FrameInfo(i).ZStep = str2double(LIFMeta.getPixelsPhysicalSizeZ(0).value);
    end

    FrameInfo(i).Time = InitialStackTime(i);

    try
      FrameInfo(i).ZGalvo = zGalvo(i);
    catch
      % This exception should be treated
    end
  end
end
