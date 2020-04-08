function [val, substr] = generateLIFMetaDataXML(Prefix, str)

    if nargin < 2
        str = 'Rot';
    end

    [lifFile, metafile] = getFiles(Prefix);
    
    
    cleanasc = getAscii(lifFile);

    [val, substr] = getMetaDataValue(cleanasc, str);

    writeMetaXML(cleanasc, metafile);

end








function writeMetaXML(text, metafile)

xmlStart = min(strfind(text, '<'));
trueXMLEnd = max(strfind(text, '>'));
wholeXML = text(xmlStart:trueXMLEnd);
metafid = fopen(metafile, 'w');
fprintf(metafid, '%c', wholeXML);
fclose(metafid);
dom = xmlread(metafile);
xmlwrite(metafile, dom);

end

function [val, substr] = getMetaDataValue(text, str)

loc = strfind(text, str);
delta = 40; %roughly the number of characters you need to get the whole key value pair
substr = text(loc:loc+delta);
% cleansubstr = strrep(substr, '.', '');
valCell = extractBetween(substr,'"', '"');
val = str2double(valCell{1});

end

function cleanasc = getAscii(lifFile)

s = dir(lifFile);         
filesize = s.bytes;

% xmlEnd = 1E6; %roughly the size of the first xml subfile in the .lif

%seems unlikely the metadata is over a fourth of the file
% xmlEnd = filesize/10;
xmlEnd = filesize;

%wrapping statement w/ evalc to suppress long, annoying output of hexdump
fid = fopen(lifFile, 'r');
[A,count] = fread(fid, xmlEnd, 'uchar');
asc = repmat('.',1, count);
idx = find(double(A)>=32);
asc(idx) = char(A(idx));
% [~, asc] = hexdump(lifFile, xmlEnd)
% evalc('[~, asc] = hexdump(lifFile, xmlEnd)'); 

%every character in the xml subfile is divided
%by 00 bytes, not sure why.
%let's get rid of them. 
cleanasc = strrep(asc, '...', '@');
% cleanasc = strrep(asc, ' ', '');
cleanasc = strrep(cleanasc, '.', '');
cleanasc = strrep(cleanasc, '@', '.');


end

function [lifFile, metafile] = getFiles(Prefix)

thisExperiment = liveExperiment(Prefix);

rawFolder = thisExperiment.rawFolder;
DLIF = dir([rawFolder,filesep,'*.lif']);
lifFile = [rawFolder, filesep, DLIF(1).name];

metaDataFolder = [rawFolder, filesep, 'MetaData'];
mkdir(metaDataFolder)
metafile = [metaDataFolder, filesep, Prefix, '.xml'];

end