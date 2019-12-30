function HisSlices = generateHisSlicesTileScan(images, NSlices, NChannels, fiducialChannel, framesIndex, tileIndex)

% For all 'nuclear' channels, generate HisSlices, and do projection
HisSlices = zeros([size(images{tileIndex}{1, 1}, 1), size(images{tileIndex}{1, 1}, 2), NSlices]);
z = 1;
firstImage = (framesIndex - 1) * NSlices * NChannels + 1 + (fiducialChannel - 1);
lastImage = framesIndex * NSlices * NChannels;

for imagesIndex = firstImage:NChannels:lastImage
    HisSlices(:, :, z) = images{tileIndex}{imagesIndex, 1};
    z = z + 1;
end

end