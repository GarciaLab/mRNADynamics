function [val, substr] = lifReader(Prefix, str)

if nargin < 2
    str = 'Rot';
end

thisExperiment = liveExperiment(Prefix);

rawFolder = thisExperiment.rawFolder;
DLIF = dir([rawFolder,filesep,'*.lif']);
lifFile = [rawFolder, filesep, DLIF(1).name];

xmlEnd = 1E6; %roughly the size of the first xml subfile in the .lif
[~, asc] = hexdump(lifFile, xmlEnd);

%every character in the xml subfile is divided by 00 bytes, not sure why.
%let's get rid of them. remember to replace this with something
%in the future that won't destroy decimal places (ie, valid '.' characters)
cleanasc = strrep(asc, '...', '@');
cleanasc = strrep(cleanasc, '.', '');
cleanasc = strrep(cleanasc, '@', '.');



loc = strfind(cleanasc, str);
delta = 40; %roughly the number of characters you need to get the whole key value pair
substr = cleanasc(loc:loc+delta);
% cleansubstr = strrep(substr, '.', '');
valCell = extractBetween(substr,'"', '"');
val = str2double(valCell{1});


%% make a full xml document

xmlStart = min(strfind(cleanasc, '<'));
trueXMLEnd = max(strfind(cleanasc, '>'));
% wholeXML = extractBetween(cleanasc, '<', '>');
wholeXML = cleanasc(xmlStart:trueXMLEnd);
sketchfile = 'X:\DorsalSynthetics\Data\RawDynamicsData\2020-01-21\1Dg-8D_EfEfEf_9_small\MetaData\sketch.xml';
sketchfid = fopen(sketchfile, 'w');
fprintf(sketchfid, '%c', wholeXML);
fclose(sketchfid);

dom = xmlread(sketchfile);
xmlwrite(dom);


end