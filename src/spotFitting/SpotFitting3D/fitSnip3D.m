function [spotsFrame, fitInfo] = fitSnip3D(spotsFrame, spotIndex, liveExperiment, imStack, nSpots)
%
% DESCRIPTION
% Sub-function called by spotFittingLoop that generates a 3D snip 
% containing a single spot and then fits 3D Gaussian(s) to that particular 
% spot. The 3D snippet is a ~1.3-1.5um neighborhood centered around the
% spot of interest.
%
% INPUT ARGUMENTS
% spotsFrame: 2D fits for all spots detected in a single frame
% spotsIndex: index of the spot within the given frame (used to index into
%             the spotsFrame structure) 
% liveExperiment: LiveExperiment instance for this particular dataset
% imStack: 3D array containing the image data for this single frame
% nSpots: number of Gaussians to fit. Should be 2 for MS2 spots (to
%         account for sister chromatids) and 1 for transcription factor
%         clusters
% 
% OPTIONS
% N/A
%
% OUPUT
% spotsFrame: 
% fitInfo: Data structure containing key fit parameter results for the
%          1 or 2 3D Gaussians fit to the spot.
%          Parameter identity is as follows: 
%               (1) amplitude of gaussian (spot 1)    
%               (2-3) xy, and z sigma values (both spots) 
%               (4-6) y,x,and z center positions (spot 1)         
%               (7) amplitude of gaussian (spot 2)
%               (8-10) y,x,and z center positions (spot 2)        
%               (11-17) inferred background gradient
%
% Author (contact): Nicholas Lammers (nlammers@berkeley.edu)
% Created: 2019-2020ish
%
% Documented by: Meghan Turner, 2022-05-31
%

FrameInfo = getFrameInfo(liveExperiment);

% extract basic fit parameters
spot = spotsFrame.Fits(spotIndex);
xSize = FrameInfo(1).PixelsPerLine;
ySize = FrameInfo(1).LinesPerFrame;
pixelSize_nm = FrameInfo(1).PixelSize*1000; %nm
zStep_nm = FrameInfo(1).ZStep*1000; %nm
zMax = FrameInfo(1).NumberSlices+1;
snipDepth = ceil(3500/zStep_nm); % sets Z depth above and below brightest plane that is extracted for fitting

% NL: Ideally we would make this independent of 2D fit info 
brightestZPlane = spot.brightestZ;
xSpot = spot.xDoG(spot.z==brightestZPlane);
ySpot = spot.yDoG(spot.z==brightestZPlane);

if isfield(spot, 'snippet_size') && ~isempty(spot.snippet_size)
    % spot.snippet_size is hardcoded in the main segmentSpots function to 
    % be 1.3um, converted to pixels
    snippet_size = spot.snippet_size; % in pixels
else
    neighboorhood_nm = 1500; % in nm, set to be around 1.5 um 
    snippet_size = round(neighboorhood_nm/pixelSize_nm); %in pixels
end
snippet_size = uint16(snippet_size(1));

%% %%%%%%%%%%%%%%%%%%%%%% generate 3D snip %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zBot = max([2,brightestZPlane - snipDepth]);
zTop = min([zMax-1, brightestZPlane + snipDepth]);
zRange = zBot:zTop;
xRange = max([1,xSpot-snippet_size]):min([xSize,xSpot+snippet_size]);
yRange = max([1,ySpot-snippet_size]):min([ySize,ySpot+snippet_size]);

snip3D = imStack(yRange, xRange, zRange);
  
%% %%%%%%%%%%%%%%%%%%%%%% perform fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get stack origins
xMin = single(min(xRange));
yMin = single(min(yRange));
zMin = single(min(zRange));

% Perform fits
fitInfo = fit3DGaussians(snip3D, pixelSize_nm, zStep_nm, nSpots, []);

% add fit info
spotsFrame.Fits(spotIndex).SpotFitInfo3D = fitInfo;%single(GaussParams1);

% combined metrics
spotsFrame.Fits(spotIndex).GaussIntensity3D = fitInfo.GaussIntegralTot;
spotsFrame.Fits(spotIndex).GaussIntensity3DSE = fitInfo.GaussIntegralSETot;
spotsFrame.Fits(spotIndex).Offset3D = fitInfo.offset;
spotsFrame.Fits(spotIndex).GaussIntensity3DRaw = fitInfo.GaussIntegralRaw;
spotsFrame.Fits(spotIndex).GaussIntensity3DRawSE = fitInfo.GaussIntegralRawSE;

% spotsFrame.Fits(spotIndex).Offset3DSE = offsetSE;   
spotsFrame.Fits(spotIndex).GaussPos3D = fitInfo.SpotCentroid + [yMin xMin zMin] - 1.0;
spotsFrame.Fits(spotIndex).GaussPos3DSE = fitInfo.SpotCentroidSE;


end