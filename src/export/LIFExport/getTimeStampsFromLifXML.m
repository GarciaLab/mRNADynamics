function stampElapsed = getTimeStampsFromLifXML(file, LIFMeta)

arguments
    file string
    LIFMeta
end

if ~isfile(file)
    error('input must be a valid file path');
end

%testing file. leave here.
% file = 'X:\Armando\scratch_small2series.xml';

fid = fopen(file, 'r');
str = string(fscanf(fid, '%c'));
fclose(fid);

timeStampStrings = extractBetween(str,"NumberOfTimeStamps","/TimeStamp");
NSeries = length(timeStampStrings);
stampList = [];
nPlanes = extractBetween(str,"NumberOfTimeStamps=""",""">");
nChannels = LIFMeta.getChannelCount(0);

nSlices = [];
for seriesIndex = 1:NSeries
  nSlices(seriesIndex) = str2double(LIFMeta.getPixelsSizeZ(seriesIndex-1));
end
  
nFrames = NPlanes./nSlices/NChannels;



for seriesIndex = 1:NSeries
    
    %refine the text string
    cleanStr = extractBetween(timeStampStrings(seriesIndex), ">", "<");
    %divide the string into an array of time stamps
    stampList = [stampList; split(cleanStr)];
    
    
end

%clean the list of empty strings
stampList(stampList=="") = [];

%the time stamps are in hex. let's convert them to decimal
%so we can do numeric operations on them.
stampList_dec = hex2dec(char(stampList));

%get elapsed time in seconds since the first frame.
%for downstream uses, let's also transpose this to
%a row vector.

conversionFactor = 1E-7; %honestly not sure why this number, but it works

stampElapsed = (conversionFactor*stampList_dec -...
    conversionFactor*stampList_dec(1))';

end
