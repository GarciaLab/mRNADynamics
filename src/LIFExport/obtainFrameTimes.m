%Obtains frames times and first time for all the frames in the series
function [Frame_Times, First_Time] = obtainFrameTimes(XMLFolder, SeriesFiles, NSeries, NFrames, NSlices, NChannels)
  frameTimesIndex = 1;
  for seriesIndex = 1:NSeries
    xDoc = xmlread([XMLFolder,filesep,SeriesFiles(seriesIndex).name]);
    TimeStampList = xDoc.getElementsByTagName('TimeStamp');
    for framesIndex = 0:(NFrames(seriesIndex) * NSlices(seriesIndex) * NChannels) - 1
        TimeStamp = TimeStampList.item(framesIndex);
        Date = char(TimeStamp.getAttribute('Date'));
        Time = char(TimeStamp.getAttribute('Time'));
        Milli = char(TimeStamp.getAttribute('MiliSeconds'));
        if contains(Time, 'AM')
            time_in_days = datenum(strcat(Date, '-', Time, '-', Milli), 'mm/dd/yyyy-HH:MM:SS AM-FFF');
        elseif contains(Time, 'PM')
            time_in_days = datenum(strcat(Date, '-', Time, '-', Milli), 'mm/dd/yyyy-HH:MM:SS PM-FFF');
        else
            error('something''s wrong with your timestamps. ask AR or HG?')
        end
        Frame_Times(frameTimesIndex) = time_in_days * 86400;
        frameTimesIndex = frameTimesIndex + 1;
    end
  end
  
  %Obtains the first time
  First_Time = Frame_Times(1);
  for timeIndex = 1:length(Frame_Times)
    if Frame_Times(timeIndex) == 0
      Frame_Times(timeIndex) = 0;
    else
      Frame_Times(timeIndex) = Frame_Times(timeIndex) - First_Time;
    end
  end
end
