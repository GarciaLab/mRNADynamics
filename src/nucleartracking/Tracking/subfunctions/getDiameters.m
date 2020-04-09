function [ diameters ] = getDiameters(FrameInfo, nFrames, indMitosis )
%GETDIAMETERS This function returns a vector with the assumed diameter at
%each frame.
%   Detailed explanation goes here

nCycles = 14;
space_resolution = getDefaultParameters(FrameInfo,'space resolution');
diameters = zeros(nFrames,1);
mitosis = indMitosis(all(indMitosis ~= 0,2),:);


try
    
    for j = size(mitosis,1):-1:1
        
        diameters(mitosis(j,1):mitosis(j,2)) = 0.5*(getDefaultParameters(FrameInfo,['d' num2str(nCycles)])+getDefaultParameters(FrameInfo,['d' num2str(max(10,nCycles-1))]))/space_resolution;
        
        if j == size(mitosis,1)
            diameters(mitosis(j,2)+1:end) = getDefaultParameters(FrameInfo,['d' num2str(nCycles)])/space_resolution;
        else
            diameters(mitosis(j,2)+1:mitosis(j+1,1)-1) = getDefaultParameters(FrameInfo,['d' num2str(nCycles)])/space_resolution;
        end
        
        nCycles = max(10,nCycles-1);
        
    end
    
    if diameters(1) == 0
        diameters(1:mitosis(1,1)) = getDefaultParameters(FrameInfo,['d' num2str(nCycles)])/space_resolution;
    end
    
catch
    
    diameters(1:end) =  getDefaultParameters(FrameInfo,['d' num2str(nCycles)])/space_resolution;

end

end

