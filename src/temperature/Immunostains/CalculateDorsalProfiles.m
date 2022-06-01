function CalculateDorsalProfiles(Prefix, varargin)


ShowPlots = false;

liveExperiment= LiveExperiment(Prefix);

CompiledEmbryoPath = [liveExperiment.resultsFolder, filesep, 'CompiledEmbryos.Mat'];
load(CompiledEmbryoPath);

Ellipses = getEllipses(liveExperiment);
FrameInfo = getFrameInfo(liveExperiment);
try
    clrmp = single(hsv(20));
    clrmp = clrmp(randperm(length(clrmp)), :);
catch
    %in case the user doesn't have this colormap, just keep going.
end
nSlices=FrameInfo(1).NumberSlices;
%%
ProcPath = liveExperiment.userProcFolder;
DropboxFolder = liveExperiment.userResultsFolder;

Channel1 = liveExperiment.Channel1;
Channel2 = liveExperiment.Channel2;
Channel3 = liveExperiment.Channel3;
Channel4 = liveExperiment.Channel4;
Channel5 = liveExperiment.Channel5;


xSize = liveExperiment.xDim;
ySize = liveExperiment.yDim;
PixelSize_um = liveExperiment.pixelSize_um;

%Get the nuclei segmentation data
Ellipses = getEllipses(liveExperiment);
FrameInfo = getFrameInfo(liveExperiment);

try
    clrmp = single(hsv(20));
    clrmp = clrmp(randperm(length(clrmp)), :);
catch
    %in case the user doesn't have this colormap, just keep going.
end
Channels = {Channel1, Channel2, Channel3, Channel4, Channel5};





InputChannelIndexes = find(contains(Channels, 'input', 'IgnoreCase', true));
HisChannelIndexes = find(contains(Channels, 'his', 'IgnoreCase', true));
ChannelsToIntegrate = unique([HisChannelIndexes InputChannelIndexes]);
if isempty(InputChannelIndexes)
    warning(['No input channel found. Check correct definition in MovieDatabase.',...
        ' Input channels should use the :input notation.'])
    return;
end
%%

close all
nEmbryos = length(FrameInfo);
%Get information about the image size
% HisImage=imread([PreProcPath,filesep,Prefix,filesep,D(1).name]);
GoodEmbryos = 1:nEmbryos;
GoodEmbryos = GoodEmbryos(CompiledEmbryos.Approved);
EllipsesFluoInfo = cell(size(Ellipses));
% %%
% close all
% if sum(InputChannelIndexes)
%     for CurrentEmbryo =1:nEmbryos
%         if ismember(CurrentEmbryo, GoodEmbryos)
%             
%             ellipseFrame = double(Ellipses{CurrentEmbryo});
%             ellipseFrame = ellipseFrame(ellipseFrame(:,1) > 0, :);
%             ellipseFrame = ellipseFrame(ellipseFrame(:,2) > 0, :);
%             Ellipses{CurrentEmbryo} = ellipseFrame;
%             IntegrationRadius=ellipseFrame(1,3)/2;       %Radius of the integration region in um
%             x = 1:floor(ceil(2*IntegrationRadius)/2)*2+1;
%             y = 1:floor(ceil(2*IntegrationRadius)/2)*2+1;
%             [X, Y] = meshgrid(x, y);
%             FluoFilter = double(sqrt((X-median(x)).^2 + (Y-median(y)).^2) < ceil(2*IntegrationRadius)/2);
%             %Initialize fields
%             EmbryoSchnitzFluo = NaN(size(ellipseFrame, 1),8+nSlices, length(Channels));
%             
%             CoordA = CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,:);
%             CoordP = CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,:);
%             CoordD = CompiledEmbryos.RotatedCoordDs(CurrentEmbryo,:);
%             CoordV = CompiledEmbryos.RotatedCoordVs(CurrentEmbryo,:);
%             %bIndex = boundary(ellipseFrame(:,1), ellipseFrame(:,2), 1);
%             APLength = CoordP(1)-CoordA(1);
%             DVLength = CoordV(2)-CoordD(2);
%             
%             for channelIndex = 1:length(Channels)
%                 
%                 EmbryoSchnitzFluo(:,1:3, channelIndex) = ellipseFrame(:,1:3);
%                 EmbryoSchnitzFluo(:,4:5, channelIndex) = ellipseFrame(:,5:6);
%                 EmbryoSchnitzFluo(:,6, channelIndex) = (EmbryoSchnitzFluo(:,1, channelIndex)-CoordA(1))/APLength;
%                 EmbryoSchnitzFluo(:,7, channelIndex) = (EmbryoSchnitzFluo(:,2, channelIndex)-CoordD(2))/DVLength;
%                 if ismember(channelIndex, ChannelsToIntegrate)
%                     EmbryoChStackFile = [liveExperiment.preFolder, filesep, Prefix,...
%                         '_Position',iIndex(CurrentEmbryo,3),'_001_ch', iIndex(channelIndex,2),'.tif'];
%                     
%                     
%                     EmbryoChStack = imreadStack2(EmbryoChStackFile, liveExperiment.yDim, liveExperiment.xDim,...
%                         liveExperiment.nFrames);
%                     zIndex = 1;
%                     if ShowPlots
%                         close all
%                         EmbryoFig = figure(1);
%                         fullAxes = axes(EmbryoFig);
%                         EmbryoImage = imshow(EmbryoChStack(:,:,zIndex), [min(min(EmbryoChStack(:,:,zIndex))), ...
%                             max(max(EmbryoChStack(:,:, zIndex)))], 'Border','Tight','Parent',fullAxes);
%                         FilterEmbryoImage = imgaussfilt(EmbryoChStack(:,:,zIndex), 2/PixelSize_um);
%                         hold(fullAxes, 'on')
%                         
%                         
%                         for k=1:size(ellipseFrame, 1)
%                             n =k;
%                             %         PlotHandle{k} = drawellipse('Center',[ellipseFrame(n, 1) ellipseFrame(n, 2)],...
%                             %             'SemiAxes',[ellipseFrame(n, 3) ellipseFrame(n, 4)], ...
%                             %             'RotationAngle',ellipseFrame(n, 5) * (360/(2*pi)), 'FaceAlpha', 0,...
%                             %             'InteractionsAllowed', 'none', 'LabelVisible', 'hover', 'Label', num2str(ellipseFrame(n, 9)));
%                             colorhash = uint8(mod(round(ellipseFrame(n, 1)+ellipseFrame(n, 2)),20)+1);
%                             PlotHandle{k} = ellipse(ellipseFrame(n, 3), ellipseFrame(n, 4),...
%                                 ellipseFrame(n, 5) * (360/(2*pi)), ellipseFrame(n, 5),...
%                                 ellipseFrame(n, 6), clrmp(colorhash,:), 10, fullAxes, 0.5);
%                             
%                         end
%                         hold(fullAxes, 'off')
%                     end
%                     
%                     
%                     
%                     %%
%                     
%                     
%                     %Create the circle that we'll use as the mask
%                     
%                     
%                     
%                     
%                     %Extract the fluroescence of each schnitz, for each channel,
%                     
%                     
%                     
%                     %         nameSuffix=['_ch',iIndex(ChN,2)];
%                     
%                     %
%                     convImage = imfilter(EmbryoChStack, FluoFilter, 'same');
%                     
%                     for j=1:length(EmbryoSchnitzFluo)
%                         for zIndex = 1:nSlices
%                             x_coord = uint16(round(EmbryoSchnitzFluo(j, 4, channelIndex)));
%                             y_coord = uint16(round(EmbryoSchnitzFluo(j, 5, channelIndex)));
%                             EmbryoSchnitzFluo(j, 8+zIndex, channelIndex) = convImage(y_coord, x_coord, zIndex)/sum(FluoFilter(:));
%                         end
%                         
%                     end %loop of nuclei in a frame
%                     EmbryoSchnitzFluo(:, 8, channelIndex) = max(EmbryoSchnitzFluo(:, 9:9+nSlices-1, channelIndex), [], 2);
%                 end
%             end
%             EllipsesFluoInfo{CurrentEmbryo} = EmbryoSchnitzFluo;
%         else
%             EllipsesFluoInfo{CurrentEmbryo} = NaN(size(0,8+nSlices, length(Channels)));
%         end
%     end
%     
%     
%     
% else
%     
%     error('Input channel not recognized. Check correct definition in MovieDatabase.Input channels should use the :input notation.');
%     
% end

%%
% EmbryoSchnitzPath = [liveExperiment.resultsFolder, filesep, 'EllipsesFluoInfo.mat'];
%save(EmbryoSchnitzPath, 'EllipsesFluoInfo')
%%

RotatedHisFile = [liveExperiment.preFolder, filesep, Prefix, '-His_Rotated.tif'];



rotatedHisMat = imreadStack2(RotatedHisFile, liveExperiment.yDim, liveExperiment.xDim,...
    liveExperiment.nFrames);
CurrentEmbryo = 28;
IntegrationRadius=1;      %Radius of the integration region in um
x = 1:floor(ceil(2*IntegrationRadius)/2)*2+1;
y = 1:floor(ceil(2*IntegrationRadius)/2)*2+1;
[X, Y] = meshgrid(x, y);
FluoFilter = double(sqrt((X-median(x)).^2 + (Y-median(y)).^2) < ceil(2*IntegrationRadius)/2);
%Initialize fields

CoordA = CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,:);
CoordP = CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,:);
CoordD = CompiledEmbryos.RotatedCoordDs(CurrentEmbryo,:);
CoordV = CompiledEmbryos.RotatedCoordVs(CurrentEmbryo,:);
%bIndex = boundary(ellipseFrame(:,1), ellipseFrame(:,2), 1);
APLength = CoordP(1)-CoordA(1);
DVLength = CoordV(2)-CoordD(2);

channelIndex = 2;


EmbryoChStackFile = [liveExperiment.preFolder, filesep, Prefix,...
    '_Rotated_Position',iIndex(CurrentEmbryo,3),'_001_ch', iIndex(channelIndex,2),'.tif'];

for i=1:nSlices
EmbryoChStack(:,:,i) = imread(EmbryoChStackFile, liveExperiment.yDim, liveExperiment.xDim,...
    i);
end
%%
MaxStack = max(EmbryoChStack,[],3);
close all
EmbryoFig = figure(1);
fullAxes = axes(EmbryoFig);
EmbryoImage = imshow(MaxStack, [min(min(MaxStack)), ...
    max(max(MaxStack))], 'Border','Tight','Parent',fullAxes);
EmbryoFig2 = figure(2);
fullAxes2 = axes(EmbryoFig2);
EmbryoImage2 = imshow(EmbryoChStack(:,:,3), [min(EmbryoChStack(:)), ...
    max(EmbryoChStack(:))], 'Border','Tight','Parent',fullAxes2);


h = fspecial3('ellipsoid', [4/PixelSize_um, 12/PixelSize_um, 1]);
EllipsoidFilteredImage = imfilter(MaxStack, h);
DorsalImage = double(EllipsoidFilteredImage);
DorsalImage(round(CoordA(2))+1:end,:,:) = 0;
BWDorsalImage = (DorsalImage-min(DorsalImage(:)))/(max(DorsalImage(:))-min(DorsalImage(:)));
APbins = .1:.0125:.9;
DVbinStarts = 0:0.01:round(0.4-10/PixelSize_um/DVLength, 2);
DVbinheight = 2.5/PixelSize_um/DVLength;
DVwindows = zeros(length(APbins)-1, 2);
DVFilter = zeros(size(MaxStack));
for i = 1:length(APbins)-1
    PixelBoundaries_AP = [round(CoordA(1)+APbins(i)*APLength), round(CoordA(1)+APbins(i+1)*APLength)];
    DorsalAPWindowImage = DorsalImage(round(CoordD(2)):round(CoordD(2)+.4*DVLength),PixelBoundaries_AP(1):PixelBoundaries_AP(2));
    DorsalAPWindowSum = sum(DorsalAPWindowImage, 2);
    SmoothedDVinfo = zeros(1, length(DVbinStarts));
    for j = 1:length(DVbinStarts)
        try
            SmoothedDVinfo(j) = ...
                sum(DorsalAPWindowSum(1+round(DVbinStarts(j)*DVLength):round(DVbinheight*DVLength+DVbinStarts(j)*DVLength)));
        catch
           SmoothedDVinfo(j) = NaN;
        end
    end
    [~, maxind ] = nanmax(SmoothedDVinfo);
    DVwindows(i, 1) = round(DVbinStarts(maxind)*DVLength)+round(CoordD(2));
    DVwindows(i,2) = round(DVbinheight*DVLength+DVbinStarts(maxind)*DVLength)+round(CoordD(2))-1;
    DVFilter(DVwindows(i, 1):DVwindows(i,2), PixelBoundaries_AP(1):PixelBoundaries_AP(2)) = 1;
end
ThresholdedImage = MaxStack;
ThresholdedImage(DVFilter == 0) = 0;

nbd_size = floor(ceil(2.5/PixelSize_um)/2)*2+mod(ceil(ceil(2.5/PixelSize_um)/2), 2); 
h = fspecial('average', nbd_size);
MeanImage = imfilter(MaxStack, h);
StdImage = stdfilt(MaxStack, ones(nbd_size, nbd_size, 'logical'));
VarImage = StdImage.^2;
NormedImage = (double(MaxStack)-double(MeanImage))./StdImage;

GaussImage = imgaussfilt(NormedImage, 2);
GaussImage(DVFilter == 0) = min(GaussImage(:));
% % Show the first image
        normedGaussMap = mat2gray(GaussImage);
        
        NormedGaussMapRange=[min(min(normedGaussMap)),max(max(normedGaussMap))];
     
%         if ShowPlots
%         NormedGaussMapFigure=figure(5);
%         set(NormedGaussMapFigure,'units', 'normalized', 'position',[0.01, .1, .7, .7]);
%         
%         normedGaussMapAxes = axes(NormedGaussMapFigure,'Units', 'normalized', 'Position', [0 0 1 1]);
%         
%         imNormedGauss = imshow(normedGaussMap,NormedGaussMapRange,'Parent',normedGaussMapAxes);
%         
%         end
        %%
        %h = fspecial('average',NeighborhoodSize);
        hisImage = normedGaussMap;
   
        
        FrameInfo = getFrameInfo(liveExperiment);
        diameters = [getDefaultParameters(FrameInfo,'d10'),getDefaultParameters(FrameInfo,'d11'),...
            getDefaultParameters(FrameInfo,'d12'),getDefaultParameters(FrameInfo,'d13'),...
            getDefaultParameters(FrameInfo,'d14')];
        
        nucleusDiameter = diameters(5)/2;
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
        [hValues,~] = histcounts(hisImage(hisImage > 0),[0:0.05:1.05]);
        [y1, x1] =max(hValues);
        x2 = 21;
        y2 = hValues(x2);
        Slope = (y2-y1)/(x2-x1);
        Intercept = y2-Slope*x2;
        xrange = x1:x2;
        height_diffs = Slope*xrange+Intercept-hValues(xrange);
        [~,auto_thresh] = max(height_diffs);
        auto_thresh = (auto_thresh + x1-2)/20;
        
        hisMask = hisImage > 0.65;
        %dilatedMask = dilatedImage > auto_thresh;
        % Find local maxima
        % Commented out GM on 9/3/20 Seems to no longer work after updates & embryoMask;
        ThreshImage2 = hisImage;
        ThreshImage2(~hisMask) = 0;
        maxima = (ThreshImage2 > hisImage);
        % Smooth the image to get more robust maxima values:
        
        maximaLinearIndices = find(maxima>0);
        [y_maxima,x_maxima] = ind2sub(size(maxima),maximaLinearIndices);
        RotatedCoordA = CompiledEmbryos.RotatedCoordAs(embryoIndex,:);
        RotatedCoordP = CompiledEmbryos.RotatedCoordPs(embryoIndex,:);
        RotatedCoordD = CompiledEmbryos.RotatedCoordDs(embryoIndex,:);
        RotatedCoordV = CompiledEmbryos.RotatedCoordVs(embryoIndex,:);
        keep_index = x_maxima >= RotatedCoordA(1) & ...
            x_maxima <= RotatedCoordP(1) & ...
            y_maxima >= RotatedCoordD(2) & ...
            y_maxima <= RotatedCoordV(2);
%       
        y_maxima = y_maxima(keep_index);
        x_maxima = x_maxima(keep_index);
        DVlength = abs(RotatedCoordV(2)-RotatedCoordD(2));
        APlength = abs(RotatedCoordA(1)-RotatedCoordP(1));
        Mid_y = RotatedCoordA(2);
        Mid_x = RotatedCoordA(1) + APlength/2;
        y_dorsal = y_maxima(y_maxima < Mid_y);
        x_dorsal = x_maxima(y_maxima < Mid_y);
        y_ventral = y_maxima(y_maxima >= Mid_y);
        x_ventral = x_maxima(y_maxima >= Mid_y);
        b_index = boundary(x_maxima, y_maxima);
        D = squareform(pdist([x_maxima y_maxima]));
        min_distances = zeros(size(x_maxima), 'double');
        APcoords = (x_maxima - RotatedCoordA(1))/APlength;
        DVcoords = (y_maxima - RotatedCoordD(2))/DVlength;
        for i = 1:length(x_maxima)
            min_distances(i) = min(D(i,D(i,:) > 0));
        end
        keep_index = ones(size(x_maxima), 'logical');
        max_dist = 2*median(min_distances);
        for i = 1:length(x_maxima)
            if min_distances > max_dist
                keep_index(i) = 0;  
            end
            if APcoords(i) > 0.1 & APcoords(i) < 0.9 & DVcoords(i) > 0.4 & DVcoords(i) < 0.6
                keep_index(i) = 0;
            end
            if APcoords(i) > 0.2 & APcoords(i) < 0.8 & DVcoords(i) > 0.25 & DVcoords(i) < 0.75
                keep_index(i) = 0;
            end
        end
        x_maxima2 = x_maxima(keep_index);
        y_maxima2 = y_maxima(keep_index);
        boundary_index = unique(boundary(x_maxima2, y_maxima2));
        x_boundary = x_maxima2(ismember(1:length(x_maxima2), boundary_index)); 
        y_boundary = y_maxima2(ismember(1:length(x_maxima2), boundary_index)); 
        x_interior = x_maxima2(~ismember(1:length(x_maxima2), boundary_index));
        y_interior = y_maxima2(~ismember(1:length(x_maxima2), boundary_index));
        all_x = [x_boundary ;  x_interior];
        all_y = [y_boundary ; y_interior];
        D2 = squareform(pdist([double(all_x)  double(all_y)]));
        D2 = D2(length(x_boundary)+1:length(all_x), 1:length(x_boundary))*PixelSize_um;
        min_distances2 = min(D2, [], 2);
        x_maxima3 = [x_boundary; x_interior(min_distances2 < 10)];
        y_maxima3 = [y_boundary; y_interior(min_distances2 < 10)];
%%
zIndex = 1;
if ShowPlots
    close all
    EmbryoFig = figure(1);
    fullAxes = axes(EmbryoFig);
    EmbryoImage = imshow(EmbryoChStack(:,:,zIndex), [min(min(EmbryoChStack(:,:,zIndex))), ...
        max(max(EmbryoChStack(:,:, zIndex)))], 'Border','Tight','Parent',fullAxes);
    FilterEmbryoImage = imgaussfilt(EmbryoChStack(:,:,zIndex), 2/PixelSize_um);
    hold(fullAxes, 'on')
    
    
    for k=1:size(ellipseFrame, 1)
        n =k;
        %         PlotHandle{k} = drawellipse('Center',[ellipseFrame(n, 1) ellipseFrame(n, 2)],...
        %             'SemiAxes',[ellipseFrame(n, 3) ellipseFrame(n, 4)], ...
        %             'RotationAngle',ellipseFrame(n, 5) * (360/(2*pi)), 'FaceAlpha', 0,...
        %             'InteractionsAllowed', 'none', 'LabelVisible', 'hover', 'Label', num2str(ellipseFrame(n, 9)));
        colorhash = uint8(mod(round(ellipseFrame(n, 1)+ellipseFrame(n, 2)),20)+1);
        PlotHandle{k} = ellipse(ellipseFrame(n, 3), ellipseFrame(n, 4),...
            ellipseFrame(n, 5) * (360/(2*pi)), ellipseFrame(n, 5),...
            ellipseFrame(n, 6), clrmp(colorhash,:), 10, fullAxes, 0.5);
        
    end
    hold(fullAxes, 'off')
end



%%
