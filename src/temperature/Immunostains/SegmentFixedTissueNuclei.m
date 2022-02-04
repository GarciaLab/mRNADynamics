function Ellipses = SegmentFixedTissueNuclei(Prefix, UseCustom)
if ~exist('UseCustom', 'var')
    UseCustom = false;
end
ShowPlots = false;
liveExperiment = LiveExperiment(Prefix);
if UseCustom
    RotatedHisFile = [liveExperiment.preFolder, filesep, Prefix, '-CustomHis_Rotated.tif'];
else
    RotatedHisFile = [liveExperiment.preFolder, filesep, Prefix, '-His_Rotated.tif'];
end

RotatedHisMat = imreadStack2(RotatedHisFile, liveExperiment.yDim,...
    liveExperiment.xDim, liveExperiment.nFrames);
CompiledEmbryoPath = [liveExperiment.resultsFolder, filesep, 'CompiledEmbryos.Mat'];
load(CompiledEmbryoPath);

PixelSize_um =liveExperiment.pixelSize_um;
xSize = size(RotatedHisMat,2);
ySize =  size(RotatedHisMat,1);
NEmbryos = size(RotatedHisMat,3);

ellipse = cell(NEmbryos,1);


%%
embryoIndex = 10; % NC < 14
embryoIndex = 18; % NC probably 14?
embryoIndex = 28; % NC = 14
embryoIndex = 37;% NC probably 14?
embryoIndex = 66;% NC probably 14?
embryoIndex = 67; % older NC14 with focus issues
embryoIndex = 75; % CF but a little rumpled
embryoIndex = 90; % older NC14

%%
for embryoIndex =1:NEmbryos
    disp(['Embryo Index: ', num2str(embryoIndex)])
    
    if CompiledEmbryos.Approved(embryoIndex)
        
        EmbryoImage = RotatedHisMat(:,:,embryoIndex);
        % EmbryoImage= histeq(mat2gray(EmbryoImage), ReferenceHist);
        %             EmbryoImage =EmbryoImage *255;
        
        
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
        
        %
        % EmbryoFigure=figure;
        % set(EmbryoFigure,'units', 'normalized', 'position',[0.01, .1, .7, .7]);
        %
        % embryoAxes = axes(EmbryoFigure,'Units', 'normalized', 'Position', [0 0 1 1]);
        % cc=1;
        %
        % % Show the first image
        % imEmbryo = imshow(EmbryoImage,DisplayRange,'Parent',embryoAxes);
        pixelValues = histc(EmbryoImage(:), [0:256]);
        ObservedPixelValues = find(pixelValues);
        ThreshImage = EmbryoImage > ObservedPixelValues(3)-1;
        
        sigma = round(2.5/PixelSize_um); % 2.5 micron sigma
        
        imGauss = imgaussfilt(EmbryoImage,sigma);
        
        %
        % GaussEmbryoFigure=figure(2);
        % set(GaussEmbryoFigure,'units', 'normalized', 'position',[0.01, .1, .7, .7]);
        %
        % gaussAxes = axes(GaussEmbryoFigure,'Units', 'normalized', 'Position', [0 0 1 1]);
        % cc=1;
        %
        % gaussRange=[min(min(imGauss)),max(max(imGauss))];
        % % Show the first image
        % imEmbryo = imshow(imGauss,gaussRange,'Parent',gaussAxes);
        %
        % %%
        
        NeighborhoodSize_um = 4; % in microns
        NeighborhoodRad = floor(round(NeighborhoodSize_um/PixelSize_um)/2);
        NeighborhoodSize = NeighborhoodRad*2+1;
        
        meanIntensityMap = zeros(ySize, xSize, 'double');
        varIntensityMap = zeros(ySize, xSize, 'double');
        normedIntensityMap = zeros(ySize, xSize, 'double');
        
        for i = 1+NeighborhoodRad:ySize-NeighborhoodRad-1
            for j =1+NeighborhoodRad:xSize-NeighborhoodRad-1
                meanIntensityMap(i,j) = mean(EmbryoImage(i-NeighborhoodRad:i+NeighborhoodRad,...
                    j-NeighborhoodRad:j+NeighborhoodRad), 'all');
                varIntensityMap(i,j) = var(EmbryoImage(i-NeighborhoodRad:i+NeighborhoodRad,...
                    j-NeighborhoodRad:j+NeighborhoodRad),0, 'all');
                
                normedIntensityMap(i,j)  = (EmbryoImage(i,j)-meanIntensityMap(i,j))/sqrt(varIntensityMap(i,j));
            end
        end
        
        %
        % MeanMapFigure=figure(3);
        % set(MeanMapFigure,'units', 'normalized', 'position',[0.01, .1, .7, .7]);
        %
        % meanMapAxes = axes(MeanMapFigure,'Units', 'normalized', 'Position', [0 0 1 1]);
        %
        %
        % MeanMapRange=[min(min(meanIntensityMap)),max(max(meanIntensityMap))];
        % % Show the first image
        % imMeanMap = imshow(meanIntensityMap,MeanMapRange,'Parent',meanMapAxes);
        %
        %
        % NormedMapFigure=figure(4);
        % set(NormedMapFigure,'units', 'normalized', 'position',[0.01, .1, .7, .7]);
        %
        % normedMapAxes = axes(NormedMapFigure,'Units', 'normalized', 'Position', [0 0 1 1]);
        % cc=1;
        %
        % NormedMapRange=[min(min(normedIntensityMap)),max(max(normedIntensityMap))];
        % % Show the first image
        % imNormedMap = imshow(normedIntensityMap,NormedMapRange,'Parent',normedMapAxes);
        
        %%
        sigma2 = round(1/PixelSize_um);
        normedGaussMap = imgaussfilt(normedIntensityMap,sigma2);
        normedGaussMap(isnan(normedGaussMap)) = min(min(normedGaussMap));
        % NormedGaussMapFigure=figure(5);
        % set(NormedGaussMapFigure,'units', 'normalized', 'position',[0.01, .1, .7, .7]);
        %
        % normedGaussMapAxes = axes(NormedGaussMapFigure,'Units', 'normalized', 'Position', [0 0 1 1]);
        % cc=1;
        %
        % % Show the first image
        normedGaussMap = mat2gray(normedGaussMap);
        
        NormedGaussMapRange=[min(min(normedGaussMap)),max(max(normedGaussMap))];
        % normedGaussMap(normedGaussMap < prctile(normedGaussMap(:),90)) = 0;
        % imNormedGauss = imshow(normedGaussMap,NormedGaussMapRange,'Parent',normedGaussMapAxes);
        %
        
        %%
        %h = fspecial('average',NeighborhoodSize);
        hisImage = normedGaussMap;
        FrameInfo = getFrameInfo(liveExperiment);
        diameters = [getDefaultParameters(FrameInfo,'d10'),getDefaultParameters(FrameInfo,'d11'),...
            getDefaultParameters(FrameInfo,'d12'),getDefaultParameters(FrameInfo,'d13'),...
            getDefaultParameters(FrameInfo,'d14')];
        
        nucleusDiameter = diameters(5);
        LoGratio = getDefaultParameters(FrameInfo,'LoGratio');
        space_resolution = getDefaultParameters(FrameInfo,'space resolution');
        localMaximumRadius = LoGratio*nucleusDiameter/space_resolution;
        LoGradius = nucleusDiameter/space_resolution*LoGratio;
        edgeClearance = getDefaultParameters(FrameInfo,'edge clearance')*nucleusDiameter/space_resolution;
        
        % Edited by GM on 1/7/2019
        xDim_um = FrameInfo(1).PixelsPerLine * FrameInfo(1).PixelSize;
        yDim_um = FrameInfo(1).LinesPerFrame * FrameInfo(1).PixelSize;
        
        
        
        %%
        localMaxMask = fspecial('disk',localMaximumRadius);
        localMaxMask = imbinarize(mat2gray(localMaxMask),graythresh(mat2gray(localMaxMask)));
        localMaxMask(round(length(localMaxMask)/2),round(length(localMaxMask)/2))  = 0;
        
        % Filter the image
        % filteredImg = fourierFilterWithSymmetricBoundaryConditions(...
        %     hisImage,-fspecial('log',round(10*LoGradius),LoGradius));
        
        dilatedImage =  imdilate(hisImage,localMaxMask);
        %         maximaMask = (hisImage > dilatedImage);
        %         G = imfilter(hisImage,fspecial('disk',3),'symmetric');
        %         [y_max2, x_max2] =ind2sub(size(maximaMask),find(maximaMask));
        %         maximaLinearIndices = find(maximaMask>0);
        %         try
        %             [v,~] = voronoin([y_max2 x_max2]);
        %         catch
        %             %not sure when this exception occurs-
        %             %it adds additional options (http://www.qhull.org/html/qh-optq.htm) -AR
        %             %Qbb - scale last coordinate to [0,m] for Delaunay (this is the
        %             %default for voronoin. it doesn't need to be included)
        %             %Qz -add a point-at-infinity for Delaunay triangulations
        %             [v,~] = voronoin([y_max2 x_max2],{'Qbb','Qz'});
        %         end
        %
        %         %get the voronoi edge pixel indices(linear) that are within the image
        %         ind_v = v(:,1) > 0.5 &...
        %             v(:,2) > 0.5...
        %             & v(:,1) < size(maximaMask,1)...
        %             & v(:,2) < size(maximaMask,2);
        %
        %         %get the corresponding subscript indices from the above linear indices
        %         background_ind = sub2ind(size(maximaMask),round(v(ind_v,1)),round(v(ind_v,2)));
        %
        %         %append the maxima to the voronoi edges. this corresponds to an image
        %         %with a bunch of voronoi edges with dots somewhere inside them
        %         sample_ind = [maximaLinearIndices(:) ; background_ind(:)];
        %
        %         %these are fluorescence values in the (smoothed) nuclear image at the
        %         %locations of sample_ind (the voronoi vertices and maxima)
        %         nuc1 = G(sample_ind);
        %
        %         %autothresholding with Otsu's method
        %         thresh = graythresh(mat2gray(nuc1));
        %         indNuclei = imbinarize(mat2gray(nuc1),thresh);
        %         indNuclei = indNuclei(1:numel(maximaLinearIndices));
        %
        %         % Get rid of all nuclei below thresh
        %         maximaLinearIndices = maximaLinearIndices(indNuclei);
        %
        %         % Convert into coordinates
        %         [y_max2,x_max2] = ind2sub(size(hisImage),maximaLinearIndices);
        %
        [hValues,~] = histcounts(dilatedImage(dilatedImage > 0),[0:0.01:1.01]);
        [y1, x1] =max(hValues);
        x2 = 101;
        y2 = hValues(x2);
        Slope = (y2-y1)/(x2-x1);
        Intercept = y2-Slope*x2;
        xrange = x1:x2;
        height_diffs = Slope*xrange+Intercept-hValues(xrange);
        [~,auto_thresh] = max(height_diffs);
        auto_thresh = (auto_thresh + x1-2)/100;
        
        dilatedMask = dilatedImage > 0.65;
        % Find local maxima
        % Commented out GM on 9/3/20 Seems to no longer work after updates & embryoMask;
        ThreshImage = hisImage;
        ThreshImage(~dilatedMask) = 0;
        maxima = (ThreshImage > dilatedImage);
        % Smooth the image to get more robust maxima values:
        
        maximaLinearIndices = find(maxima>0);
        [y_maxima,x_maxima] = ind2sub(size(maxima),maximaLinearIndices);
        
        if ShowPlots
            close all
            figure(1)
            imshow(ThreshImage)
            hold on
            scatter(x_maxima,y_maxima,'r.')
            %scatter(x_max2,y_max2,'y.')
            hold off
            
            figure(2)
            imshow(EmbryoImage, DisplayRange)
            hold on
            scatter(x_maxima,y_maxima,'r.')
            %scatter(x_max2,y_max2,'y.')
            hold off
            
        end
        %Check whether anything was found in this image. If nothing is found,
        %this is usually the result of an empty frame
        % if isempty(xm) || isempty(ym)
        %    error(['No nuclei found in frame ',num2str(frameNumber),'. Check that that frame is not blank.'])
        % end
        
        
        
        
        
        %(y, x, major axis, minor axis, orientation angle, maxcontourvalue, time,
        %particle_id %optionally, schnitz id
        %only the first 4 columns are used, so we'll just set the rest to 0.
        radiusScale = 1.0;
        centroid_orig = TransformToOriginalImageCoordinates([x_maxima,y_maxima], Prefix, CompiledEmbryos, liveExperiment, embryoIndex);
        
        nEllipses = length(x_maxima);
        ellipse{embryoIndex} = zeros(nEllipses,8);
        
        for ellipseIndex = 1:nEllipses
            
            centroid = [x_maxima(ellipseIndex), y_maxima(ellipseIndex)];
            majorAxis = nucleusDiameter*radiusScale;
            minorAxis = nucleusDiameter*radiusScale;
            centroid2 = centroid_orig(ellipseIndex,:);
            orientationAngle = 0;
            maxContourValue = 0;
            ellipse{embryoIndex}(ellipseIndex,:) = [centroid, majorAxis, minorAxis,...
                centroid2, orientationAngle, maxContourValue];
        end
    else
        ellipse{embryoIndex} = zeros(0,8);
        
    end
    
    
end

%%
Ellipses = ellipse;
ellipsesFile = [liveExperiment.resultsFolder,filesep,'Ellipses.mat'];
save2(ellipsesFile, Ellipses);


