function processMovieChannel(ExperimentType, channelIndex, numberOfFrames, Prefix, OutputFolder,...
    movieImages, framesIndex, seriesIndex, NChannels,...
    NSlices, coatChannel, inputProteinChannel, zslicesPadding, lowbit)

dataType = 'uint8';

experimentType1 = (strcmpi(ExperimentType,'1spot') || strcmp(ExperimentType,'2spot') || strcmp(ExperimentType,'2spot1color')) && channelIndex == coatChannel;
experimentType2 = strcmpi(ExperimentType,'2spot2color') || strcmpi(ExperimentType, 'inputoutput');
experimentType3 = strcmpi(ExperimentType, 'input') && sum(channelIndex == inputProteinChannel);

%Create a blank image
BlankImage = false(size(movieImages{1}{1,1}));

% if zPadding was indicated in the arguments, we round up to the series
% with more z-slices (because we'll pad with blank images the other series)
if (zslicesPadding)
    topZSlice = max(NSlices);
else
    % if no zPadding, we round down to the series with less z-slices
    topZSlice = min(NSlices);
end

if(experimentType1 || experimentType2 || experimentType3)
    % Save the blank image at the beginning of the stack only for
    % type1,2,or3
    NameSuffix = ['_ch',iIndex(channelIndex,2)];
    NewName = [Prefix, '_', iIndex(numberOfFrames,3), '_z', iIndex(1,2), NameSuffix, '.tif'];
    imwrite(BlankImage, [OutputFolder, filesep, NewName]);
end

%Copy the rest of the images
if(experimentType1 || experimentType2 || experimentType3)
    slicesCounter = 1;
    firstImage = (framesIndex-1) * NSlices(seriesIndex) * NChannels + 1 + (channelIndex - 1);
    lastImage = framesIndex * NSlices(seriesIndex) * NChannels;
    if firstImage == lastImage
        firstImage = 1;
        lastImage = 1;
    end
    for imageIndex = firstImage:NChannels:lastImage
        if slicesCounter <= topZSlice
            % if zPadding, it will process all images (because topZSlice would be max(NSlices)
            % if no zPadding, it will process images rounding down to the series with least
            % zSlices, because topZSlice would be min(NSlices)
            NewName = [Prefix, '_', iIndex(numberOfFrames,3), '_z', iIndex(slicesCounter + 1, 2), NameSuffix, '.tif'];
            im = movieImages{seriesIndex}{imageIndex,1};
            if max(im(:)) < 256 %max uint8 value
                im = uint8(im);
            end
            imwrite(im, [OutputFolder, filesep, NewName]);
            slicesCounter = slicesCounter + 1;
        end
    end
    
    % Save as many blank images at the end of the stack are needed
    % (depending on zPadding being active or not)
    for zPaddingIndex = slicesCounter+1:topZSlice+2
        NewName = [Prefix, '_', iIndex(numberOfFrames,3), '_z', iIndex(zPaddingIndex, 2), NameSuffix, '.tif'];
        imwrite(BlankImage, [OutputFolder, filesep, NewName]);
    end
    
end
end
