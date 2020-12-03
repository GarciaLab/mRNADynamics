%Get the time stamp corresponding to the first slice of each Z-stack
function [InitialStackTime, zPosition] = getFirstSliceTimestamp(NSlices, NSeries,...
    NPlanes, NChannels, Frame_Times, XMLFolder, seriesXML)

zPosition = [];
frameIndex= 1;

m = 1;

for series = 1:NSeries
    
    if series == 1
        startIndex = 1;
    else
        startIndex = sum(NPlanes(1:series-1)) + 1;
    end
    
    frameStep = NSlices(series)*NChannels;
    finalIndex = sum(NPlanes(1:series));
    
    try
        %Grab the z-galvo position by parsing the XML metadata
        xmltext = fileread([XMLFolder,filesep,seriesXML(series).name]);
        expressionobj = '(?<=" ZPosition=").*?(?=")';
        possiblePositions = regexp(xmltext, expressionobj, 'match');
        zPos = str2double(possiblePositions{1}); %AR- this is sometimes the right element.
    catch
        warning('didn''t record zgalvo position since data hasn''t been exported from LASX');
    end
    
    for frame = startIndex:frameStep:finalIndex
        
        InitialStackTime(frameIndex) = Frame_Times(frame);
        frameIndex = frameIndex + 1;
        
        %currently only correctly records the z-position from data from the
        %Bateman lab Leica
        try
            zPosition(m) = zPos;
        catch
            warning('didn''t record zgalvo position')
        end
        
        m = m + 1;
        
    end
    
end

end
