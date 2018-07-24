function [ diameters ] = getDiameters( nFrames, indMitosis )
%GETDIAMETERS This function returns a vector with the assumed diameter at
%each frame.
%   Detailed explanation goes here

space_resolution = getDefaultParameters('space resolution');
diameters = zeros(nFrames,1);


mitosis = indMitosis(all(indMitosis ~= 0,2),:);

nCycle = 14;
for j = size(mitosis,1):-1:1
    diameters(mitosis(j,1):mitosis(j,2)) = 0.5*(getDefaultParameters(['d' num2str(nCycle)])+getDefaultParameters(['d' num2str(max(10,nCycle-1))]))/space_resolution;
    if j == size(mitosis,1)
        diameters(mitosis(j,2)+1:end) = getDefaultParameters(['d' num2str(nCycle)])/space_resolution;
    else
        diameters(mitosis(j,2)+1:mitosis(j+1,1)-1) = getDefaultParameters(['d' num2str(nCycle)])/space_resolution;
    end
    nCycle = max(10,nCycle-1);
end
if diameters(1) == 0
    diameters(1:mitosis(1,1)) = getDefaultParameters(['d' num2str(nCycle)])/space_resolution;
end
end

