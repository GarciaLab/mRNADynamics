% Added PreferredFileName so we can automate testing and bypass the user prompt when there are many files available.
function FrameInfo = processLIFExportMode(rawDataFolder, ProjectionType, Channels,...
    Prefix, PreProcFolder, PreferredFileNameForTest,...
    nuclearGUI, skipExtraction, lowbit)

markandfind = false;

%Loads file and metadata
[XMLFolder, seriesPropertiesXML, seriesXML] = getSeriesFiles(rawDataFolder);

if ~isempty(strfind(seriesPropertiesXML(1).name, 'Mark_and_Find'))
    markandfind = true;
end

[LIFImages, LIFMeta] = loadLIFFile(rawDataFolder);

%Obtains frames information
[NSeries, NFrames, NSlices, NPlanes, NChannels, Frame_Times] = getFrames(LIFMeta);

if sum(NFrames)~=0
    [Frame_Times, First_Time] = obtainFrameTimes(XMLFolder, seriesPropertiesXML, NSeries, NFrames, NSlices, NChannels);
    [InitialStackTime, zPosition] = getFirstSliceTimestamp(NSlices, NSeries, NPlanes, NChannels, Frame_Times, XMLFolder, seriesXML);
else
    InitialStackTime = [];
    zPosition = [];
end

FrameInfo = recordFrameInfo(NFrames, NSlices, InitialStackTime, LIFMeta, zPosition);

if markandfind
    FrameInfo = repmat(FrameInfo, NSeries, 1);
end


%Find the flat field (FF) information
LIFExportMode_flatFieldImage(LIFMeta, rawDataFolder, PreProcFolder, Prefix, PreferredFileNameForTest);

if sum(NFrames) == 0
    NFrames = ~NFrames;
end

if ~skipExtraction
    %Copy the data
    waitbarFigure = waitbar(0, 'Extracting LIFExport images');
    
    
    %Counter for number of frames
    numberOfFrames = 1;
    %Load the reference histogram for the fake histone channel
    load('ReferenceHist.mat')
    if nuclearGUI
        [Channels, ProjectionType] = chooseNuclearChannels(...
            LIFImages, NSeries, NSlices, NChannels, NFrames, ProjectionType, Channels, ReferenceHist);
    end
    
    ySize = size(LIFImages{1}{1,1}, 1);
    xSize = size(LIFImages{1}{1,1}, 2);
%     BlankImage = uint16(zeros(ySize, xSize));
    
    nPadding = 2;
    movieMat = zeros(NChannels, max(NSlices)+nPadding, sum(NFrames), ySize, xSize, 'uint16'); % ch z t x y
   hisMat = zeros(sum(NFrames), ySize, xSize, 'uint16'); % f x y

    
%     zslicesPadding = false;
    
    for seriesIndex = 1:NSeries
        waitbar(seriesIndex/NSeries, waitbarFigure)
        for framesIndex = 1:NFrames(seriesIndex)
            
            for channelIndex = 1:NChannels
                
               
                topZSlice = min(NSlices);
                
                
%                 NameSuffix = ['_ch',iIndex(channelIndex,2)];
%                 NewName = [Prefix, '_', iIndex(numberOfFrames,3), '_z', iIndex(1,2), NameSuffix, '.tif'];
%                 imwrite(BlankImage, [PreProcFolder, filesep, NewName]);
                
                %Copy the rest of the images
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
                        movieMat(channelIndex, slicesCounter + 1, numberOfFrames, :, :) = LIFImages{seriesIndex}{imageIndex,1};
%                         NewName = [Prefix, '_', iIndex(numberOfFrames,3), '_z', iIndex(slicesCounter + 1, 2), NameSuffix, '.tif'];
%                         imwrite(LIFImages{seriesIndex}{imageIndex,1}, [PreProcFolder, filesep, NewName]);
                        slicesCounter = slicesCounter + 1;
                    end
                end
                
                % Save as many blank images at the end of the stack are needed
                % (depending on zPadding being active or not)
%                 for zPaddingIndex = slicesCounter+1:topZSlice+2
%                     NewName = [Prefix, '_', iIndex(numberOfFrames,3), '_z', iIndex(zPaddingIndex, 2), NameSuffix, '.tif'];
%                     imwrite(BlankImage, [PreProcFolder, filesep, NewName]);
%                 end
                %
                %                 processMovieChannel(channelIndex, numberOfFrames, Prefix, OutputFolder,...
                %                     LIFImages, framesIndex, seriesIndex, NChannels, NSlices,...
                %                     zslicesPadding, lowbit);
            end
            
            %Now copy nuclear tracking images
            hisMat(numberOfFrames, :, :) = generateNuclearChannel(numberOfFrames, LIFImages,...
                framesIndex, seriesIndex, NSlices, NChannels,ProjectionType,...
                Channels, ReferenceHist, PreProcFolder, Prefix, lowbit);
            
            numberOfFrames = numberOfFrames + 1;
        end
    end
    
    save([PreProcFolder,filesep, Prefix, '_movieMat.mat'],'movieMat', '-v7.3', '-nocompression');
    save([PreProcFolder, filesep, Prefix, '_hisMat.mat'],'hisMat', '-v7.3', '-nocompression');

    close(waitbarFigure)
    
end

end
