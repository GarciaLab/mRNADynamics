function Spots = fit3DGaussiansToAllSpots(prefix)

[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
    Folder, Prefix, ExperimentType,Channel1,Channel2,OutputFolder, Channel3] = readMovieDatabase(prefix);
DataFolder=[DropboxFolder,filesep,prefix];

load([DataFolder,filesep,'Spots.mat'])
load([DataFolder,filesep,'FrameInfo.mat'])

xSize = FrameInfo(1).PixelsPerLine;
ySize = FrameInfo(1).LinesPerFrame;
pixelSize = FrameInfo(1).PixelSize*1000; %nm

nFrames = length(Spots);
ch = 1;
zMax = FrameInfo(1).NumberSlices+2;

snips3D = {};

% parpool(12);

parfor frame = 1:nFrames %frames
    
    nSpotsPerFrame = length(Spots(frame).Fits);
    
    for spot = 1:nSpotsPerFrame
        
        s = Spots(frame).Fits(spot);
        bZ = s.brightestZ;
        zInd = s.z==bZ;
        xSpot = s.xDoG(zInd);
        ySpot = s.yDoG(zInd);
        
        if isfield(s, 'snippet_size') && ~isempty(s.snippet_size)
            snippet_size = s.Fits(spot).snippet_size;
        else
            snippet_size = 13; %pixels
        end
        
        snipDepth = 2;
        zBot = bZ - snipDepth;
        zTop = bZ + snipDepth;
        width = 200/pixelSize; %nm. empirically determined and seems to work width of spot psf
        offsetGuess = mean(s.Offset);
        snip3D = [];
        try
            k = 1;
            for z = zBot:zTop
                
                FullSlice=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(frame,3)...
                    ,'_z' iIndex(z,2) '_ch' iIndex(ch,2) '.tif']);
                
                snip3D(:,:,k) = double(FullSlice(max(1,ySpot-snippet_size):min(ySize,ySpot+snippet_size),...
                    max(1,xSpot-snippet_size):min(xSize,xSpot+snippet_size))); %#ok<*SAGROW>
                k = k + 1;
                
            end
        catch
            snipDepth = zMax;
            for z = 1:zMax
                FullSlice=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(frame,3)...
                    ,'_z' iIndex(z,2) '_ch' iIndex(ch,2) '.tif']);
                
                snip3D(:,:,z) = double(FullSlice(max(1,ySpot-snippet_size):min(ySize,ySpot+snippet_size),...
                    max(1,xSpot-snippet_size):min(xSize,xSpot+snippet_size))); %#ok<*SAGROW>
            end
        end
        %         snips3D = [snips3D, currentSnippet3D]; %#ok<*AGROW>
        initial_params = [max(max(max(snip3D))), NaN,NaN, snipDepth + 1, width,offsetGuess];
        [fits, intensity] = fitGaussian3D(snip3D, initial_params);
        Spots(frame).Fits(spot).fits3D = fits;
        Spots(frame).Fits(spot).gauss3DIntensity = intensity;
    end
    
end

save([DataFolder,filesep,'Spots.mat'])

end

