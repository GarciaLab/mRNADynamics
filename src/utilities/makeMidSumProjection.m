function projection = makeMidSumProjection(imageStack)
    
    nSlices = size(imageStack, 3);
    lowerSlice = round(nSlices * .33);
    upperSlice =round(nSlices * .66);
    projection = sum(imageStack(:, :, lowerSlice:upperSlice), 3);
    
end