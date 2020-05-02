function metadataXML = generateLIFMetaDataXML(in, out)

%input here can be either a Prefix or a direct 
%path to a .lif file.  
%
%the output is optional. if specified, this is the
%path to where the metadata will be saved. else,
%the metadata gets saved in the project's 
%raw metadata folder. 

  arguments
        in char
        out char = '';
  end

    if ~isfile(in)
        [lifFile, metafile] = getFiles(in);
    else
        lifFile = in;
    end
    
    if ~isempty(out)
        metafile = out;
    end
    
    
      
    cleanasc = getAscii(lifFile);

    metadataXML = writeMetaXML(cleanasc, metafile);

end








function wholeXML = writeMetaXML(text, metafile)
  
    arguments
        text char 
        metafile char
    end
  
xmlStart = min(strfind(text, '<'));
trueXMLEnd = max(strfind(text, '>'));
wholeXML = text(xmlStart:trueXMLEnd);

if ~isempty(metafile)
    metafid = fopen(metafile, 'w');
    fprintf(metafid, '%c', wholeXML);
    fclose(metafid);
    dom = xmlread(metafile);
    xmlwrite(metafile, dom);
end

end

function cleanasc = getAscii(lifFile)

s = dir(lifFile);         
filesize = s.bytes;

xmlEnd = min(100E6, filesize); %100MB should cover even very large files

fid = fopen(lifFile, 'r');
[A,count] = fread(fid, xmlEnd, 'uchar');

asc = repmat('.',1, count);
idx = find(double(A)>=32);
asc(idx) = char(A(idx));

%every character in the xml subfile is divided
%by 00 bytes, not sure why.
%let's get rid of them. 
cleanasc = strrep(asc, '...', '@');
cleanasc = strrep(cleanasc, '.', '');
cleanasc = strrep(cleanasc, '@', '.');


end

function [lifFile, metafile] = getFiles(Prefix)

liveExperiment = LiveExperiment(Prefix);

rawFolder = liveExperiment.rawFolder;
DLIF = dir([rawFolder,filesep,'*.lif']);
lifFile = [rawFolder, filesep, DLIF(1).name];

metaDataFolder = [rawFolder, filesep, 'MetaData'];
mkdir(metaDataFolder)
metafile = [metaDataFolder, filesep, Prefix, '.xml'];

end