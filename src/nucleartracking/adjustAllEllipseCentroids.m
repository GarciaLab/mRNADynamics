function Ellipses = adjustAllEllipseCentroids(Prefix, varargin)


%function description- we want to use existing ellipses, but mildly adjust the position and
%diameter to better model the nuclei

nWorkers = 1;
displayFigures = false;
min_rad_um = 1; % set min and max acceptable area for nucleus segmentation
max_rad_um = 7; %this needs to be 6um for nc12. 4um for nc14
bouncingFig = [];
T = [];
tileFig = [];




for i = 1:(numel(varargin)-1)
    if i ~= numel(varargin)
        if ~ischar(varargin{i+1})
            eval([varargin{i} '=varargin{i+1};']);
        end
    end
end

disp('Adjusting ellipses...');

[~,~,DropboxFolder,~,PreProcPath]=...
    DetermineLocalFolders(Prefix);

load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo');

[xDim, yDim, pixelSize_nm, ~, ~,...
    nFrames, ~, ~] = getFrameInfoParams(FrameInfo);

PixelSize_um = pixelSize_nm / 1000;
minRad_px =  min_rad_um / PixelSize_um;
maxRad_px = max_rad_um / PixelSize_um;
minArea_px = round(pi*minRad_px^2);
maxArea_px = round(pi*maxRad_px^2);

load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], 'Ellipses');

hisFile = [PreProcPath, filesep, Prefix, filesep, Prefix, '_hisMat.mat'];
hisMat = double(loadHisMat(hisFile));


if displayFigures
    tileFig = figure('Units', 'normalized', 'Position',[0.6441 0.0744 0.3184 0.3844]);
    T = tiledlayout('flow', 'TileSpacing', 'none', 'Padding', 'none')
    bouncingFig = figure;
end

EllipsesNew = Ellipses;

parfor f = 1:nFrames
    
    disp(['Adjusting frame: ', num2str(f), '...']);
    
    hisImage = squeeze(hisMat(:, :, f));
        
    rad0 = nanmedian(Ellipses{f}(:, 3), 1);  %3 and 4 are the semimajor/minor axes of the ellipses
        
    %make sure rad0 is set to a reasonable value
    rad0 = min(maxRad_px, rad0);
    rad0 = max(minRad_px, rad0);
    
    %correct any mistakes in Ellipses
    
    Ellipses{f}(:, 3) = min(maxRad_px, Ellipses{f}(:, 3));
    Ellipses{f}(:, 3) = min(maxRad_px, Ellipses{f}(:, 3));
    

    
    ellipseMask0 = double(makeNuclearMask(Ellipses{f}, [yDim, xDim], 'radiusScale', 1.3)); %make a mask with the initial ellipse configuration
    
    hisImage = wiener2( hisImage);
    
    smoothSigma = rad0/2;
    hisSmooth = imgaussfilt( hisImage, smoothSigma, 'Padding',0); %denoise a little
    
    smoothMasked =  hisSmooth.* ellipseMask0;
    
    hisMasked =  hisImage.* ellipseMask0;
    
    
    label0 = bwlabel( ellipseMask0);
    
    nRegions = max(max( label0));
    
    for r = 1:nRegions
        
        maskTemp = ( label0 == r) .*  hisSmooth;
        %find the centroids as local maxima
        [~,ind] = max(maskTemp,[],'all','linear');
        %get the subscripts from the linear indices of the maxima
        [y1, x1] = ind2sub([yDim, xDim],ind);
        
        
        b = regionprops( label0 == r, 'BoundingBox', 'Area');
        
        if b.Area > minArea_px && b.Area < maxArea_px
            
            br = round(b.BoundingBox);
            c = imcrop(imgaussfilt( hisImage,8), br);
            do = imcrop( label0 == r, br);
            snakesFun = @(b, s, sigma) activecontour(c, do, 'Chan-Vese', 'ContractionBias', b, 'SmoothFactor', s);
            snakey = snakesFun(.4, 0, 0);
            boundaryCell = bwboundaries(snakey, 8, 'noholes');
            
            
            if ~isempty(boundaryCell)
                xs = boundaryCell{1}(:, 1);
                ys = boundaryCell{1}(:, 2);
                [~,~, Rfit]= circfit(xs,ys);
                if Rfit > minRad_px && Rfit < maxRad_px && ~isnan(Rfit)
                    EllipsesNew{f}(r, 3) = Rfit;
                    EllipsesNew{f}(r, 4) = Rfit;
                end
            end
        end
        
        
        
        
        if displayFigures
            set(0, 'currentfigure', bouncingFig);
            imagesc(maskTemp, [min(min(maskTemp(maskTemp>0))), max(maskTemp(:))]);
            hold on;
            plot(x1, y1, 'xk', 'LineWidth', 2);
            hold off;
            drawnow;
        end
        
        EllipsesNew{f}(r, 1) = x1;
        EllipsesNew{f}(r, 2) = y1;
        
        
        
    end
    if displayFigures, imagescUpdate(nexttile(T), smoothMasked, []); drawnow; end
    
    if displayFigures, imagescUpdate(nexttile(T), hisMasked, []); drawnow; end
    
    ellipseNewMask = double(makeNuclearMask(EllipsesNew{f}, [yDim, xDim], 'radiusScale', 1.3)); %make a mask with the initial ellipse configuration
    if displayFigures, imagescUpdate(nexttile(T), ellipseNewMask, []); drawnow; end
    
    smoothNewMasked =  hisSmooth.* ellipseNewMask;
    if displayFigures, imagescUpdate(nexttile(T), smoothNewMasked, []); drawnow; end
    
    
    hisNewMasked =  hisImage.* ellipseNewMask;
    if displayFigures, imagescUpdate(nexttile(T), hisNewMasked, []); drawnow; end
    
    
end

Ellipses = EllipsesNew;

save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], 'Ellipses', '-v6');

end