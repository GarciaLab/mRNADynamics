function metadataXML = generateLIFMetaDataXML(in, out, xmlSize)

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
        xmlSize double = 100E6;
  end
  
  warning('off','MATLAB:MKDIR:DirectoryExists')

    disp('Exporting LIF MetaData as XML ...')
    if ~isfile(in)
        [lifFile, metafile] = getFiles(in);
    else
        lifFile = in;
    end
    
    if ~isempty(out)
        metafile = out;
    end
    
    
      
    cleanasc = getAscii(lifFile, xmlSize);
    
    try
        metadataXML = writeMetaXML(cleanasc, metafile);
    catch
        %occasionally, the original xml file length was too long and
        %produced an error. let's try again with a shorter xml. this is nec
        %essary because it's hard to determine the true offset of the 
        %end of the xml within the .lif binary. 
        shorterXMLSize = xmlSize/2;
        metadataXML = generateLIFMetaDataXML(in, out, shorterXMLSize);
        return;
    end

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

function cleanasc = getAscii(lifFile, optionalXMLSize)

s = dir(lifFile);         
filesize = s.bytes;

xmlSize = 100E6; 
if nargin > 1
    xmlSize = optionalXMLSize;
end
xmlEnd = min(xmlSize, filesize); %100MB should cover even very large files

fid = fopen(lifFile, 'r');
[A,count] = fread(fid, xmlEnd, 'uchar');
fclose(fid);

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