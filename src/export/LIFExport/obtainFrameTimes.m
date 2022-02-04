%Obtains frames times and first time for all the frames in the series
function Frame_Times = obtainFrameTimes(XMLFolder,...
    seriesPropertiesXML, NSeries, NFrames, NSlices, NChannels)

days2seconds = 86400;
frameTimesIndex = 1;

for seriesIndex = 1:NSeries
    
    xDoc = xmlread([XMLFolder,filesep,seriesPropertiesXML(seriesIndex).name]);
    TimeStampList = xDoc.getElementsByTagName('TimeStamp');
    
    finalIndex = (NFrames(seriesIndex) * NSlices(seriesIndex) * NChannels) - 1;
    
    for framesIndex = 0 : finalIndex
        
        TimeStamp = TimeStampList.item(framesIndex);
        
        if ~isempty(TimeStamp)
            Date = char(TimeStamp.getAttribute('Date'));
            Time = char(TimeStamp.getAttribute('Time'));
            Milli = char(TimeStamp.getAttribute('MiliSeconds'));
        end
        try
        if contains(Time, 'AM')
            time_in_days = datenum(strcat(Date, '-', Time, '-', Milli), 'mm/dd/yyyy-HH:MM:SS AM-FFF');
        elseif contains(Time, 'PM')
            time_in_days = datenum(strcat(Date, '-', Time, '-', Milli), 'mm/dd/yyyy-HH:MM:SS PM-FFF');
        else
            error('something''s wrong with your timestamps. ask AR or HG?')
        end
        
        Frame_Times(frameTimesIndex) = time_in_days;
        catch
            Frame_Times(frameTimesIndex) = datenum('01/01/1900-12:00:00', 'mm/dd/yyyy-HH:MM:SS');
        end
        
        frameTimesIndex = frameTimesIndex + 1;
        
    end
    
end

Frame_Times = days2seconds * (Frame_Times - Frame_Times(1));

end