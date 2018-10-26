function FrameInfo = recordFrameInfo(NFrames, NSlices, InitialStackTime, LIFMeta, zPosition)

  for i = 1:sum(NFrames)
      
    FrameInfo(i).LinesPerFrame = str2double(LIFMeta.getPixelsSizeY(0));
    FrameInfo(i).PixelsPerLine = str2double(LIFMeta.getPixelsSizeX(0));
    FrameInfo(i).NumberSlices = min(NSlices);    
    FrameInfo(i).FileMode = 'LIFExport';
    FrameInfo(i).Time = InitialStackTime(i);
    
    %This is to allow for backwards compatibility with BioFormats
    if ~isempty(str2num(LIFMeta.getPixelsPhysicalSizeX(0)))
        FrameInfo(i).PixelSize = str2num(LIFMeta.getPixelsPhysicalSizeX(0));
        FrameInfo(i).ZStep = str2double(LIFMeta.getPixelsPhysicalSizeZ(0));
    else
        FrameInfo(i).PixelSize = str2num(LIFMeta.getPixelsPhysicalSizeX(0).value);
        FrameInfo(i).ZStep = str2double(LIFMeta.getPixelsPhysicalSizeZ(0).value);
    end

    %currently only correctly records frameinfo for data from the Bateman lab Leica
    try
      FrameInfo(i).zPosition = zPosition(i);
    catch
%       warning('didn''t record zgalvo in frameinfo- getframeinfo')
    end
    
  end
  
end
