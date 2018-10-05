%Get the time stamp corresponding to the first slice of each Z-stack
function [InitialStackTime, zPosition] = getFirstSliceTimestamp(NSlices, NSeries, NPlanes, NChannels, Frame_Times, XMLFolder, SeriesFiles3)
  
  m = 1;
  for i = 1:NSeries
    if i == 1
        StartIndex = 1;
    else
        StartIndex = sum(NPlanes(1:i-1)) + 1;
    end

    %Grab the z-galvo position by parsing the XML metadata
    xmltext = fileread([XMLFolder,filesep,SeriesFiles3(i).name]);
    expressionobj = '(?<=" ZPosition=").*?(?=")';
    possiblePositions = regexp(xmltext, expressionobj, 'match'); 
    zPos = str2double(possiblePositions{1}); %AR- I think this is the right element. 
   
    for j = StartIndex:(NSlices(i)*NChannels):sum(NPlanes(1:i))
      InitialStackTime(m) = Frame_Times(j);
%       try
        zPosition(m) = zPos;
%       catch
%         warning('didn''t record zgalvo position')
%       end      
      m = m + 1;
    end
  end
end
