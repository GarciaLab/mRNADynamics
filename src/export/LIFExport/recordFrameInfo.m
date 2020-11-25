function FrameInfo = recordFrameInfo(NFrames, NSlices, InitialStackTime, LIFMeta, zPosition)

%initialize frameinfo. this is important for xyz movies that don't have a
%time dimension. 

  FrameInfo(1).LinesPerFrame = str2double(LIFMeta.getPixelsSizeY(0));
  FrameInfo(1).PixelsPerLine = str2double(LIFMeta.getPixelsSizeX(0));
  FrameInfo(1).NumberSlices = min(NSlices);    
  FrameInfo(1).FileMode = 'LIFExport';
  FrameInfo(1).Time = 0;
  FrameInfo(1).Series = 1;
  if ~isempty(str2double(LIFMeta.getPixelsPhysicalSizeX(0))) &...
          ~isnan(str2double(LIFMeta.getPixelsPhysicalSizeX(0)))
      FrameInfo(1).PixelSize = str2double(LIFMeta.getPixelsPhysicalSizeX(0));
      FrameInfo(1).ZStep = str2double(LIFMeta.getPixelsPhysicalSizeZ(0));
  else
      try
          FrameInfo(1).PixelSize = str2double(LIFMeta.getPixelsPhysicalSizeX(0).value);
          FrameInfo(1).ZStep = str2double(LIFMeta.getPixelsPhysicalSizeZ(0).value);
      catch %no idea man
          FrameInfo(1).PixelSize = str2double(LIFMeta.getPixelsPhysicalSizeX(1).value);
          FrameInfo(1).ZStep = str2double(LIFMeta.getPixelsPhysicalSizeZ(1).value);
      end
  end
  if ~isempty(zPosition)
    FrameInfo(1).zPosition = zPosition(1);
  end
  
  for i = 2:sum(NFrames)
      
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
        try
            FrameInfo(i).PixelSize = str2num(LIFMeta.getPixelsPhysicalSizeX(0).value);
            FrameInfo(i).ZStep = str2double(LIFMeta.getPixelsPhysicalSizeZ(0).value);
        catch %no idea man
            FrameInfo(i).PixelSize = str2num(LIFMeta.getPixelsPhysicalSizeX(1).value);
            FrameInfo(i).ZStep = str2double(LIFMeta.getPixelsPhysicalSizeZ(1).value);
        end
    end

    seriesIndex = find(cumsum(NFrames)>=i,1);
    FrameInfo(i).Series = seriesIndex;
    if ~isempty(zPosition)        
        FrameInfo(i).zPosition = zPosition(seriesIndex);
    end
    
  end
  
end
