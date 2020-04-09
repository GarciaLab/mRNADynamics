function ellipsesFrameNew = adjustEllipseCentroidsFrame(ellipsesFrame, im)    

    d = struct;
    
    yDim = size(im, 1); xDim = size(im, 2);
    
    ellipsesFrameNew = ellipsesFrame; 
    
    d.hisImage = double(im);
        
    rad0 = median(ellipsesFrame(:, 3), 1);  %3 and 4 are the semimajor/minor axes of the ellipses
    
    d.ellipseMask0 = double(makeNuclearMask(ellipsesFrame, [yDim, xDim])); %make a mask with the initial ellipse configuration
    
    smoothSigma = rad0/2;
    d.hisSmooth = imgaussfilt(d.hisImage, smoothSigma, 'Padding',0); %denoise a little
    
    d.label0 = bwlabel(d.ellipseMask0);
    
    nRegions = max(max(d.label0));
        
    for r = 1:nRegions
        
        maskTemp = (d.label0 == r) .* d.hisSmooth;
        [~,ind] = max(maskTemp,[],'all','linear');
        [y1, x1] = ind2sub([yDim, xDim],ind);
        
             
        ellipsesFrameNew(r, 1) = uint16(x1);
        ellipsesFrameNew(r, 2) = uint16(y1);
        
    end
        
end
    