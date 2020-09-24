function nuclearStack = generateHisSlices(images, NSlices, NChannels,...
    fiducialChannel, framesIndex, seriesIndex)

%Leica SP8 data only

% For all 'nuclear' channels, generate HisSlices
nuclearStack = zeros([size(images{seriesIndex,1}{1, 1}, 1),...
    size(images{seriesIndex}{1, 1}, 2), NSlices(seriesIndex)]);
z = 1;
firstImage = (framesIndex - 1) * NSlices(seriesIndex) * NChannels + 1 + (fiducialChannel - 1);
lastImage = framesIndex * NSlices(seriesIndex) * NChannels;

for imagesIndex = firstImage:NChannels:lastImage
    nuclearStack(:, :, z) = images{seriesIndex}{imagesIndex, 1};
    z = z + 1;
end

end