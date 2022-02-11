function tile_array = NewTileArrayFromMetadata(Prefix, ID)
    if ~exist('Prefix')
        FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
        Dashes=strfind(FolderTemp,filesep);
        Prefix=FolderTemp((Dashes(end)+1):end);
    end    
    [SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
        DetermineLocalFolders(Prefix);
    %Find out the date it was taken
    Dashes=findstr(Prefix,'-');
    Date=Prefix(1:Dashes(3)-1);
    EmbryoName=Prefix(Dashes(3)+1:end);
    if ~isempty(strfind(lower(ID), 'mid'))
        filename = 'MidTile';
    elseif ~isempty(strfind(lower(ID), 'surf'))
        filename = 'SurfTile';
    else
        filename = ID;
    end
    LIFPath = [SourcePath, filesep, Date, filesep, EmbryoName,filesep,...
    'FullEmbryo\', filename,'.lif'];

    r = bfGetReader();
    % Decorate the reader with the Memoizer wrapper
    r = loci.formats.Memoizer(r);
    r.setId(LIFPath);
    LIFImages = bfopen(LIFPath);
    LIFMeta = LIFImages{:, 4};
    r.close();
    
    [Prefix, SkipFrames, ProjectionType, PreferredFileNameForTest, keepTifs,...
    generateTifStacks, nuclearGUI, skipExtraction, rootFolder, zslicesPadding,...
    lowbit] = exportDataForLivemRNA_processInputParameters(Prefix);
    [rawDataPath, ~, DropboxFolder, ~, PreProcPath, rawDataFolder, Prefix, ExperimentType, Channel1, Channel2, ~,...
        Channel3] = readMovieDatabase(Prefix,'rootFolder', rootFolder);
    [NSeries, NFrames, NSlices, NPlanes, NChannels, Frame_Times] = getFrames(LIFMeta);
    if sum(NFrames)~=0
        [Frame_Times, First_Time] = obtainFrameTimes(XMLFolder, seriesPropertiesXML, NSeries, NFrames, NSlices, NChannels);
        [InitialStackTime, zPosition] = getFirstSliceTimestamp(NSlices, NSeries, NPlanes, NChannels, Frame_Times, XMLFolder, seriesXML);
    else
        InitialStackTime = [];
        zPosition = [];
    end
    framesIndex  = 1;
    NTiles = NSeries;
    FrameInfo = recordFrameInfo(NFrames, NSlices, InitialStackTime, LIFMeta, zPosition);

    [coatChannel, histoneChannel, fiducialChannel, inputProteinChannel, FrameInfo] =...
    LIFExportMode_interpretChannels(ExperimentType, Channel1, Channel2, Channel3, FrameInfo);

    


    framesIndex = 1;
    tiles = {};
    zstacks = {};
    for n=1:NTiles
        ImageSlices = generateHisSlicesTileScan(LIFImages, NSlices(n),...
            NChannels, fiducialChannel, framesIndex, n);
        temp = max(ImageSlices, [], 3);
        zstacks{n} = uint16(ImageSlices);
        tiles{n} = uint16(temp);
    end


    PixelSize = double(LIFMeta.getPixelsPhysicalSizeX(1).value);% units: microns
    PixelSize_m = double(PixelSize)*10^(-6);
    ypos = [];
    xpos = [];
    for i=0:(NTiles-1)
        xpos(length(xpos)+1) = -1*double(LIFMeta.getPlanePositionX(i,0).value);
        ypos(length(ypos)+1) = double(LIFMeta.getPlanePositionY(i,0).value);
    end
    uxpos = sort(unique(xpos), 'ascend');
    uypos = sort(unique(ypos), 'ascend');
    xdim = length(uxpos); ydim = length(uypos);
    dx = round(abs((uxpos(xdim)-uxpos(xdim-1))/PixelSize_m), 0);
    dy = round(abs((uypos(ydim)-uypos(ydim-1))/PixelSize_m), 0);
    tile_array.tiles = tiles;
    tile_array.zstacks = zstacks;
    tile_array.rows = {};
    tile_array.cols = {};
    tile_array.grid_positions = {};
    tile_array.imgs = {};
    tile_array.use_tiles = [];
    sigma = .6/PixelSize;
    for i=1:NTiles
        ri = round((xpos(i)-uxpos(1))/PixelSize_m, 0)+1;
        ci = round((ypos(i)-uypos(1))/PixelSize_m, 0)+1;
        gri = find(xpos(i) == uxpos);
        gci = find(ypos(i) == uypos);
        tile_array.rows{i} = ri;
        tile_array.cols{i} = ci;
        tile_array.grid_positions{i} = [gri, gci];
        tile_array.imgs{i} = imfilter(tiles{i},... 
            fspecial('gaussian', 2, 1), 'symmetric');
        tile_array.use_tiles{i} = true;
    end

    tile_array.prevrows = {};
    tile_array.prevcols = {};
    
       %Save the information
       
    %Create the output folder if it doesn't exist 
    outputFolder = [DropboxFolder,filesep,Prefix,filesep,'FullEmbryoStitching'];
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    

    saveVars = {};
    saveVars = [saveVars, 'tile_array'];
    outputDatafile = [upper(ID(1)), ID(2:end), 'TileArray.mat'];
    save([outputFolder, filesep, outputDatafile],saveVars{:});
    GenerateStitchedData(Prefix, ID);
end