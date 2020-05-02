function out = ellipsesSizeUnchanged(EllipsesOld, EllipsesNew)
%
%check if routine modified the number of ellipses or the frames of ellipses
out = struct;

out.nFramesOld = length(EllipsesOld);
out.nEllipsesOld = nan(nFramesOld, 1);

for frameIndex = 1:length(EllipsesOld)
    
    out.nEllipsesOld(frameIndex) = size(EllipsesOld{frameIndex}, 1); 

end

out.nFramesNew = length(EllipsesNew);
out.nEllipsesNew = nan(nFramesNew, 1);

for frameIndex = 1:length(EllipsesNew)
    
    out.nEllipsesNew(frameIndex) = size(EllipsesNew{frameIndex}, 1); 

end

if out.nFramesNew ~= out.nFramesOld
    error('something terrible happened. number of frames changed in Ellipses')
end

if out.nEllipsesNew ~= out.nEllipsesNew
    error('something terrible happened. number of ellipses changed in Ellipses')
end
