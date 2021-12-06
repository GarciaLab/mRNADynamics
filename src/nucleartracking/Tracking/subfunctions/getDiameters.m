function [ diameters ] = getDiameters(FrameInfo, nFrames, indMitosis, varargin )
%GETDIAMETERS This function returns a vector with the assumed diameter at
%each frame.
%   Detailed explanation goes here

nCycles = 14;
space_resolution = getDefaultParameters(FrameInfo,'space resolution');
diameters = zeros(nFrames,1);
mitosis = indMitosis(all(indMitosis ~= 0,2),:);

if ~isempty(varargin)
    anaphaseFrames = varargin{1};
end


try
    for j = size(mitosis,1):-1:1
        
        diameters(mitosis(j,1):mitosis(j,2)) =...
            0.5*(getDefaultParameters(FrameInfo,['d' num2str(nCycles)])+ ...
            getDefaultParameters(FrameInfo,['d' num2str(max(10,nCycles-1))]))/space_resolution;
        
        if j == size(mitosis,1)
            diameters(mitosis(j,2)+1:end) = getDefaultParameters(FrameInfo,['d' num2str(nCycles)])/space_resolution;
        else
            diameters(mitosis(j,2)+1:mitosis(j+1,1)-1) =...
                getDefaultParameters(FrameInfo,['d' num2str(nCycles)])/space_resolution;
        end
        
        nCycles = max(10,nCycles-1);
        
    end
    
    if diameters(1) == 0
        diameters(1:mitosis(1,1)) = getDefaultParameters(FrameInfo,['d' num2str(nCycles)])/space_resolution;
    end
catch
    ncFrames = [];
    extendedAnaphaseFrames = [zeros(1,8),anaphaseFrames'];
    for f = 1:nFrames
        nc = find(extendedAnaphaseFrames > f, 1) - 1;
        if isempty(nc)
            % MT, 2021-12-06: Quick, inelegant fix to allow for nuclei 
            % segmentation of single frames
            if nFrames == 1
                ncFrames = nCycles;
            else
                ncFrames(f:nFrames) = max(ncFrames) + 1;
            end
            break;
        else
            ncFrames(f) = nc;
        end
    end
    
    for f = 1:nFrames
        
        diameters(f) = getDefaultParameters(FrameInfo,['d' num2str(ncFrames(f))])/space_resolution;
        
    end
    
end


end