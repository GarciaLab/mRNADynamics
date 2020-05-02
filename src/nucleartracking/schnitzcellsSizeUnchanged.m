function out = schnitzcellsSizeUnchanged(schnitzcellsOld, schnitzcellsNew)
%
%check if routine modified the number of schnitzcells or the number of frames of
%schnitzcells

out = struct;

out.nSchnitzesOld = length(schnitzcellsOld);
out.nFramesOld = nan(nFramesOld, 1);

for schnitzIndex = 1:out.nSchnitzesOld
    
    out.nFramesOld(schnitzIndex) = length(schnitzcellsOld(schnitzIndex).frames);

end


out.nSchnitzesNew = length(schnitzcellsNew);
out.nFramesNew = nan(nFramesNew, 1);

for schnitzIndex = 1:out.nSchnitzesNew
    
    out.nFramesNew(schnitzIndex) = length(schnitzcellsNew(schnitzIndex).frames);

end

if out.nFramesNew ~= out.nFramesOld
    error('something terrible happened. number of frames changed in schnitzcells')
end

if out.nSchnitzesOld ~= out.nSchnitzesNew
    error('something terrible happened. number of schnitzcells changed in Schnitzcells')
end
