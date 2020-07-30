% Added PreferredFileName so we can automate testing and bypass the user prompt when there are many files available.
function FrameInfo = processLIFExportMode(rawDataFolder, ProjectionType, Channels,...
    Prefix, PreProcFolder, PreferredFileNameForTest,...
    nuclearGUI, skipExtraction,...
    shouldExportNuclearProjections, shouldExportMovieFiles)

disp('Exporting movie file...');

cleanupObj = onCleanup(@myCleanupFun);

liveExperiment = LiveExperiment(Prefix);

resultsFolder = liveExperiment.resultsFolder;

if ~shouldExportMovieFiles
    FrameInfo = [];
end
moviePrecision = 'uint16';
hisPrecision = 'uint16';

%Load the reference histogram for the fake histone channel
load('ReferenceHist.mat', 'ReferenceHist');

shouldMakeMovieMat = shouldExportMovieFiles ||...
    (shouldExportNuclearProjections &&...
    ~exist([PreProcFolder,filesep, Prefix, '_movieMat.mat'],'file'));

if shouldMakeMovieMat
    
    markandfind = false;
    
    %%
    %this section is being deprecated
    try
        %Loads file and metadata
        [XMLFolder, seriesPropertiesXML, seriesXML] = getSeriesFiles(rawDataFolder);
    catch
        XMLFolder = '';
        seriesPropertiesXML = '';
        seriesXML = '';
    end
    
    try
        if contains(seriesPropertiesXML(1).name, 'Mark_and_Find')
            markandfind = true;
        end
    catch % do nothing
    end
    %%
    [LIFImages, LIFMeta] = loadLIFFile(rawDataFolder);
    
    %obtain precision from meta data
%     moviePrecision = char(LIFMeta.getPixelsType(0).getValue());
%     hisPrecision = char(LIFMeta.getPixelsType(0).getValue());
%     
    %Obtains frames information
    [NSeries, NFrames, NSlices,...
        NPlanes, NChannels, Frame_Times] = getFrames(LIFMeta);
    
    %use the old method(exported from lasx) if the files are exported
    %already. if they're not, just use bioformats. the lasx method is being
    %deprecated.
    if ~isempty(XMLFolder)
        timeStampRetrievalMethod = 'lasx';
    else
        timeStampRetrievalMethod = 'bioformats';
    end
    
    if sum(NFrames)~=0
        
        switch timeStampRetrievalMethod
            
            case 'manual'
                
                xml_file = [liveExperiment.rawFolder, filesep, 'lifMeta.xml'];
                
                if ~exist(xml_file, 'file')
                    generateLIFMetaDataXML(Prefix, xml_file);
                end
                
                Frame_Times = getTimeStampsFromLifXML(xml_file);
                
            case 'lasx'
                
                Frame_Times = obtainFrameTimes(XMLFolder, seriesPropertiesXML,...
                    NSeries, NFrames, NSlices, NChannels);
                
            case 'bioformats'
                
                Frame_Times = getFrameTimesFromBioFormats(LIFMeta, NSlices);
                
            otherwise, error('what?');
                
        end
    end
    
    [InitialStackTime, zPosition] = getFirstSliceTimestamp(NSlices,...
        NSeries, NPlanes, NChannels, Frame_Times, XMLFolder, seriesXML);
    
    FrameInfo = recordFrameInfo(NFrames, NSlices, InitialStackTime, LIFMeta, zPosition);
    
    if markandfind
        FrameInfo = repmat(FrameInfo, NSeries, 1);
    end
    
    save([resultsFolder, filesep, 'FrameInfo.mat'], 'FrameInfo', '-v6');
    
    %Find the flat field (FF) information
    LIFExportMode_flatFieldImage(LIFMeta,...
        rawDataFolder, PreProcFolder, Prefix, PreferredFileNameForTest);
    
    if sum(NFrames) == 0
        NFrames = ~NFrames;
    end
    
end

if ~skipExtraction
    
    if shouldMakeMovieMat
        %Copy the data
        waitbarFigure = waitbar(0, 'Extracting LIFExport images');
        
        
        %Counter for number of frames
        numberOfFrames = 1;
        
        
        ySize = size(LIFImages{1}{1,1}, 1);
        xSize = size(LIFImages{1}{1,1}, 2);
        BlankImage = cast(zeros(ySize, xSize),moviePrecision);
        
        nPadding = 2;
        hisMat = zeros(ySize, xSize, sum(NFrames), hisPrecision);
    end
    
    
    if shouldExportMovieFiles
        
        topZSlice = min(NSlices);
        
        movieMat = zeros(ySize, xSize,...
            max(NSlices)+nPadding, sum(NFrames),NChannels, moviePrecision);
        
        
        for seriesIndex = 1:NSeries
            waitbar(seriesIndex/NSeries, waitbarFigure)
            for framesIndex = 1:NFrames(seriesIndex)
                
                for channelIndex = 1:NChannels
                    
                    NameSuffix = ['_ch',iIndex(channelIndex,2)];
                    
                    NewName = [Prefix, '_', iIndex(numberOfFrames,3),...
                        NameSuffix, '.tif'];
                    
                    imwrite(BlankImage, [PreProcFolder, filesep, NewName]);
                    %
                    %Copy the rest of the images
                    slicesCounter = 1;
                    firstImageIndex = (framesIndex-1) * NSlices(seriesIndex) * NChannels +...
                        1 + (channelIndex - 1);
                    lastImageIndex = framesIndex * NSlices(seriesIndex) * NChannels;
                    if firstImageIndex == lastImageIndex
                        firstImageIndex = 1;
                        lastImageIndex = 1;
                    end
                    for imageIndex = firstImageIndex:NChannels:lastImageIndex
                        if slicesCounter <= topZSlice
                            % if zPadding, it will process all images (because topZSlice would be max(NSlices)
                            % if no zPadding, it will process images rounding down to the series with least
                            % zSlices, because topZSlice would be min(NSlices)
                            imSlice = cast(LIFImages{seriesIndex}{imageIndex,1},moviePrecision);
                            movieMat(:, :,slicesCounter + 1,  numberOfFrames,...
                                channelIndex) = imSlice;
                            imwrite(imSlice,...uint16(LIFImages{seriesIndex}{imageIndex,1}),...
                                [PreProcFolder, filesep, NewName], 'WriteMode', 'append');
                            slicesCounter = slicesCounter + 1;
                        end
                    end
                    
                    %Save as many blank images at the end of the stack are needed
                    %(depending on zPadding being active or not)
                    for zPaddingIndex = slicesCounter+1:topZSlice+2
                        imwrite(BlankImage, [PreProcFolder, filesep, NewName], 'WriteMode', 'append');
                    end
                end
                
                
                
                %Now create nuclear projection movies
                if ~nuclearGUI
                    
                    hisMat(:, :, numberOfFrames) = generateNuclearChannel(...
                        numberOfFrames, LIFImages,...
                        framesIndex, seriesIndex, NSlices, NChannels,ProjectionType,...
                        Channels, ReferenceHist, PreProcFolder, Prefix);
                    
                    saveNuclearProjection(hisMat, [PreProcFolder, filesep, Prefix, '-His.tif']);
                    
                end
                
                numberOfFrames = numberOfFrames + 1;
            end
        end
        
        
    end
    
    
    if nuclearGUI && shouldExportNuclearProjections
        
        if ~exist('movieMat', 'var')
            movieMat = getMovieMat(LiveExperiment(Prefix));
        end
        
        chooseAnaphaseFrames(...
            Prefix, 'ProjectionType', ProjectionType,...
            'ReferenceHist', ReferenceHist, 'movieMat', movieMat);
        
    end
    
    try close(waitbarFigure); catch; end
    
end


disp('Movie files exported.');


end