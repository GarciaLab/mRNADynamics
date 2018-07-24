function Stackno = Lattice_FindFromFilename(Filename, pretag, posttag)
% finds the Stack number or other properties from the file naming convention
%
% possible pretags
% 'stack'   stacknumber
% 'ch'      channelnumber
%
% possible posttags
% 'msec'    time of acquisition of this stack
% 'ms'      frame exposure time

pretagpos = strfind(Filename, ['_',pretag])+length(['_',pretag]);
posttagpos = strfind(Filename, [posttag,'_']);
% find the closest distance
tagdistances = zeros(length(posttagpos), length(pretagpos));
for i=1:length(posttagpos)
    tagdistances(i,:) = posttagpos(i)-pretagpos;
end
tagdistances(tagdistances<0) = max(tagdistances(:)); % we do not want negative distances
[val, ind] = min(tagdistances(:));
postind = mod(ind, length(posttagpos));
if postind==0
    postind = length(posttagpos);
end
preind = ceil(ind/length(posttagpos));

    
Stackno = str2double(Filename(pretagpos(preind):posttagpos(postind)-1));