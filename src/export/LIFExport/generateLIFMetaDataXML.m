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

xmlEnd = 1E6; %roughly the size of the first xml subfile in the .lif

evalc('[~, asc] = hexdump(lifFile, xmlEnd)');  %wrapping statement w/ evalc to suppress long, annoying output of hexdump

%every character in the xml subfile is divided by 00 bytes, not sure why.
%let's get rid of them. remember to replace this with something
%in the future that won't destroy decimal places (ie, valid '.' characters)
cleanasc = strrep(asc, '...', '@');
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