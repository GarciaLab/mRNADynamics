function Spots = fit3DGaussiansToAllSpots(prefix, varargin)

optionalResults = '';

segmentSpots = 0;
displayFigures = 0;
nWorkers = 8;

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'displayFigures')
        displayFigures = 1;
    elseif strcmpi(varargin{i}, 'segmentSpots')
        Spots = varargin{i+1};
        segmentSpots = 1;
    elseif strcmpi(varargin{i}, 'optionalResults')
        optionalResults = varargin{i+1};
    elseif strcmpi(varargin{i}, 'nWorkers')
        nWorkers = varargin{i+1};
    else
    end
end

[~,ProcPath,DropboxFolder,~, PreProcPath,...
    ~, Prefix, ~,Channel1,Channel2,~, Channel3] = readMovieDatabase(prefix, optionalResults);

spotChannels = getCoatChannel(Channel1, Channel2, Channel3);

DataFolder=[DropboxFolder,filesep,prefix];

if ~segmentSpots
    load([DataFolder,filesep,'Spots.mat'], 'Spots');
end

load([DataFolder,filesep,'FrameInfo.mat'], 'FrameInfo');

xSize = FrameInfo(1).PixelsPerLine;
ySize = FrameInfo(1).LinesPerFrame;
pixelSize = FrameInfo(1).PixelSize*1000; %nm
zStep = FrameInfo(1).ZStep;

zMax = FrameInfo(1).NumberSlices+2;

startParallelPool(nWorkers, displayFigures);    

if ~iscell(Spots)
    Spots = {Spots};
end

for ch = spotChannels

    SpotsCh = Spots{ch};
    nFrames = length(SpotsCh);

        
    parfor frame = 1:nFrames %frames
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
            
            zoomFactor = 1; %replace this with zStep from FrameInfo later- AR 1/31/2019
            snipDepth = round(3*zoomFactor);
            zBot = bZ - snipDepth;
            zTop = bZ + snipDepth;
            width = 200/pixelSize; %nm. empirically determined and seems to work width of spot psf
            offsetGuess = nanmean(s.Offset);
            snip3D = [];
            dogSnip3D = [];
%           try
                k = 1;
                for z = zBot:zTop
                    if z > 1 && z < zMax
                        FullSlice=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(frame,3)...
                            ,'_z' iIndex(z,2) '_ch' iIndex(ch,2) '.tif']);
                        nameSuffix = ['_ch', iIndex(ch, 2)];
                        dog_name = ['DOG_', Prefix, '_', iIndex(frame, 3), '_z', iIndex(z, 2), nameSuffix, '.tif'];
                        dog_full_path = [ProcPath, filesep,Prefix,'_',filesep,'dogs',filesep,dog_name];
                         FullDoGSlice= double(imread(dog_full_path));
                        snip3D(:,:,k) = double(FullSlice(max(1,ySpot-snippet_size):min(ySize,ySpot+snippet_size),...
                            max(1,xSpot-snippet_size):min(xSize,xSpot+snippet_size))); %#ok<*SAGROW>
                        dogSnip3D(:, :, k) = double(FullDoGSlice(max(1,ySpot-snippet_size):min(ySize,ySpot+snippet_size),...
                            max(1,xSpot-snippet_size):min(xSize,xSpot+snippet_size))); %#ok<*SAGROW>
                        k = k + 1;
                    end
                end
                
%                 dogSnip3D = (dogSnip3D/10) - 100; 
                
%                 catch
%                     snipDepth = zMax;
%                     for z = 1:zMax
%                         FullSlice=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(frame,3)...
%                             ,'_z' iIndex(z,2) '_ch' iIndex(ch,2) '.tif']);
% 
%                         snip3D(:,:,z) = double(FullSlice(max(1,ySpot-snippet_size):min(ySize,ySpot+snippet_size),...
%                             max(1,xSpot-snippet_size):min(xSize,xSpot+snippet_size))); %#ok<*SAGROW>
%                     end
%                 end
            %         snips3D = [snips3D, currentSnippet3D]; %#ok<*AGROW>
            
            initial_params = [max(snip3D(:)), snipDepth + 1, width,offsetGuess];
            fitOptions = {};
            if displayFigures
                fitOptions = [fitOptions,'displayFigures'];
            end
            
            [fits, intensity, ci95, intensityError95] = fitGaussian3D(snip3D, initial_params, zStep, pixelSize, fitOptions{:});

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
            SpotsCh(frame).Fits(spot).gauss3DIntensityCI95 = intensityError95;
            
            midind = ceil(size(dogSnip3D,3)/2);
            if size(dogSnip3D,3) > 2
                SpotsCh(frame).Fits(spot).ampdog3 = sum(sum(sum(dogSnip3D(:,:,midind-1:midind+1))) );
                SpotsCh(frame).Fits(spot).ampdog3Max = max(max(max(dogSnip3D(:,:,midind-1:midind+1))) );
            else
                SpotsCh(frame).Fits(spot).ampdog3 = NaN;
                SpotsCh(frame).Fits(spot).ampdog3Max = max(max(max(dogSnip3D(:,:,:))));
            end
            
            %this is a flag that the fit was done over few z-frames so the
            %user can decide if they want to keep the fit or not
            if k < 3
                SpotsCh(frame).Fits(spot).weeFit = 1;
            else
                SpotsCh(frame).Fits(spot).weeFit = 0;
            end
            
       end
    end
    
    Spots{ch} = SpotsCh;
    
end

if length(Spots) < 2
    Spots = Spots{1};
end


save([DataFolder,filesep,'Spots.mat'],'Spots', '-v7.3');
disp('3D fitting done on all spots.')

end