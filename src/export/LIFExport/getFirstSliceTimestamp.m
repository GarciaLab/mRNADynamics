%Get the time stamp corresponding to the first slice of each Z-stack
function [InitialStackTime, zPosition] = getFirstSliceTimestamp(NSlices, NSeries,...
    NPlanes, NChannels, Frame_Times, XMLFolder, seriesXML)
  
  zPosition = [];
  frameIndex= 1;
 
  for series = 1:NSeries
    
    if series == 1
        startIndex = 1;
    else
        startIndex = sum(NPlanes(1:series-1)) + 1;
    end
        
    frameStep = NSlices(series)*NChannels;
    finalIndex = sum(NPlanes(1:series));
    
    for frame = startIndex:frameStep:finalIndex
      
      InitialStackTime(frameIndex) = Frame_Times(frame);      
      frameIndex = frameIndex + 1;
    
    end
    
  end
  
end
