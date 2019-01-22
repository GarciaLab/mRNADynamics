function Spots = fit3DGaussiansToAllSpots(prefix, varargin)

[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
    Folder, Prefix, ExperimentType,Channel1,Channel2,OutputFolder, Channel3] = readMovieDatabase(prefix);
DataFolder=[DropboxFolder,filesep,prefix];

segmentSpots = 0;
displayFigures = 0;
for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'displayFigures')
        displayFigures = 1;
    elseif strcmpi(varargin{i}, 'segmentSpots')
        Spots = varargin{i+1};
        segmentSpots = 1;
    else
    end
end

if ~segmentSpots
    load([DataFolder,filesep,'Spots.mat'], 'Spots');
end

load([DataFolder,filesep,'FrameInfo.mat'], 'FrameInfo');

xSize = FrameInfo(1).PixelsPerLine;
ySize = FrameInfo(1).LinesPerFrame;
pixelSize = FrameInfo(1).PixelSize*1000; %nm
zstep = FrameInfo(1).ZStep;

nCh = 1;
if iscell(Spots)
    nCh = length(Spots);
else
    Spots = {Spots};
end

zMax = FrameInfo(1).NumberSlices+2;


for ch = 1:nCh
    
    nFrames = length(Spots{ch});
    SpotsCh = Spots{ch};
    for frame = 1:nFrames %frames
        nSpotsPerFrame = length(SpotsCh(frame).Fits);
        SpotsFrame = SpotsCh(frame).Fits;
            
       for spot = 1:nSpotsPerFrame

            s = SpotsFrame(spot);
            bZ = s.brightestZ;
            zInd = s.z==bZ;
            xSpot = s.xDoG(zInd);
            ySpot = s.yDoG(zInd);

            if isfield(s, 'snippet_size') && ~isempty(s.snippet_size)
                snippet_size = s.snippet_size;
            else
                snippet_size = 13; %pixels
            end

            snipDepth = 2;
            zBot = bZ - snipDepth;
            zTop = bZ + snipDepth;
            width = 200/pixelSize; %nm. empirically determined and seems to work width of spot psf
            offsetGuess = nanmean(s.Offset);
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
            if displayFigures
                initial_params = [max(max(max(snip3D))), NaN,NaN, snipDepth + 1, width,offsetGuess];
                [fits, intensity, ci95] = fitGaussian3D(snip3D, initial_params, zstep,'displayFigures');
            else
                initial_params = [max(max(max(snip3D))), NaN,NaN, snipDepth + 1, width,offsetGuess];
                    [fits, intensity, ci95] = fitGaussian3D(snip3D, initial_params, zstep);
            end

            x = fits(2) - snippet_size + xSpot;
            y = fits(3) - snippet_size + ySpot;
            z = fits(4) - snipDepth + bZ;
            
            dxLow = ci95(2, 1) - snippet_size + xSpot;
            dxHigh = ci95(2, 2) - snippet_size + xSpot;
            dyLow = ci95(3, 1) - snippet_size + ySpot;
            dyHigh = ci95(3, 2) - snippet_size + ySpot;
            dzLow = ci95(4, 1) - snipDepth + bZ;
            dzHigh = ci95(4, 2) - snipDepth + bZ;
            
            SpotsCh(frame).Fits(spot).GaussPos = [x,y,z];
            SpotsCh(frame).Fits(spot).GaussPosCI95 = [dxLow,dxHigh; dyLow, dyHigh; dzLow, dzHigh];
            SpotsCh(frame).Fits(spot).fits3D = fits;
            SpotsCh(frame).Fits(spot).gauss3DIntensity = intensity;
            SpotsCh(frame).Fits(spot).fits3DCI95 = ci95;

       end
    end
    
    Spots{ch} = SpotsCh;
    
end

if length(Spots) < 2
    Spots = Spots{1};
end

save([DataFolder,filesep,'Spots.mat'],'Spots', '-v7.3');
disp('Fitting done.')

end

