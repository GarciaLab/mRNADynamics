function nucMask=GetNuclearMask(I,smoothing,imopenParam)

%The following has been canibalized from Michael Tikhonov's FISH code.
%It creates a quick nuclear mask so we cna check alignments


I_mask = getEmbryoMaskLive(I, max(10*smoothing,10));

% Restrict area of consideration to the inside of the embryo
if isempty(I_mask)
    I_mask = getEmbryoMaskLive(I, max(10*smoothing,10));
end
if sum(I_mask(:))/numel(I_mask)<0.1
    I_mask = true(size(I_mask));
end
% Use watershed algorithm to partition the image; prevent oversegmentation by first
% gaussian filtering the image 
    segmentationOk = false;
    while ~segmentationOk
        if smoothing>0
            iSmooth = imfilter(I,...
                fspecial('gaussian',smoothing*2, smoothing),'symmetric');
        else
            iSmooth = I;
        end
        imSegmented = single(watershed(imcomplement(iSmooth))).*single(I_mask);
        %imSegmented = watershed(imcomplement(mexican_hat_filter(I,smoothing))).*single(I_mask);

        % Areas that are too large or too small do not correspond to single nuclei; remove
        % them.
        cc=bwconncomp(imSegmented);
        if cc.NumObjects<5000
            segmentationOk = true;
        else 
            smoothing = smoothing+5;
        end
    end
    
    P=regionprops(imSegmented,'Area','PixelIdxList');
    areas = [P.Area];
    bounds = mean(areas)*[0.5, 2.5];
    
%     % show a diagnostic plot
%     f = figure;
%     hist(areas);
%     a = axis;
%     hold on;
%     plot(bounds([1,1]), a(3:4),'r-');
%     plot(bounds([2,2]), a(3:4),'r-');
%     pause;
%     close(f);
    
    imSegmented(imSegmented==0)=length(areas)+1;
    goodSegments = areas > bounds(1); % & areas < bounds(2);
    
    % For each good segment, choose a local threshold value using graythresh on the
    % original image. Then use the list of local thresholds to turn the original image
    % into a nuclear mask 
    %
    % In case of overpartitioning (a likely event for early nuclear cycles), some
    % thresholds may be chosen too low and the telltale sign of this is that the points
    % that are white after thresholding are right next to the edge of the region of
    % partitionining (as opposed to forming a blob in the center of the region). If this
    % occurs, mark the corresponding threshold as possibly too low and replace it by the
    % mean of good thresholds.
    
    threshList = zeros([1,length(areas)+1]);
    % catches all pixels that were not inside any segment
    threshList(length(areas)+1) = Inf;
    for i=1:length(areas)
        if goodSegments(i)
            % grayThresh only uses 256 levels, so make sure we use full dynamic range.
            offset = min(I(P(i).PixelIdxList));
            im = I(P(i).PixelIdxList) - offset;
            scale = 255/max(double(im(:)));
            im = uint8(immultiply(im, scale));
            threshList(i) = offset + 255 * graythresh(im) / scale;
        else
            threshList(i) = Inf;
        end
    end    
    % Remove abnormally low thresholds (usually caused by proximity to the boundary)
    % For this, declare that [mean - sigma] is the lower bound and replace
    % thresholds values that are too low by this value.
    
%     % Threshold the image for the first time
%     nucMask = imSegmented .* (I > threshList(imSegmented));
%     % This nucMask is imSegmented seen through our first tentative mask which may have
%     % some regions under-thresholded. To find under-thresholded regions, check for
%     % regions that are, after application of the mask, immediately adjacent to the
%     % border of the partitioning region. For a correctly chosen threshold this should
%     % not happen. Therefore, dilate this nucMask by 1 pixel and check if any region
%     % starts overlapping with the borders (pixels in imSegmented equal to the dummy
%     % segment number length(areas)+1) 
%     nucMask = imdilate(nucMask,strel('square',2));
%     tooLow = unique(nucMask(imSegmented(:)==(length(areas)+1)));
%     if tooLow(1)==0
%         % Exclude the zero
%         tooLow = tooLow(2:end);
%     end
%     %threshList(tooLow) = Inf;
%    %threshList(tooLow) = minThresh;
%    %threshList(threshList == Inf) = mean(validThresh);
    validThresh = threshList(threshList<Inf);
    minThresh = mean(validThresh)-std(validThresh);
    threshList(threshList < minThresh) = minThresh;
    
    nucMask = (I > threshList(imSegmented));
    
    % fill holes
    nucMask = imfill(nucMask,'holes');
    
    if imopenParam>0
        % remove small specks    
        nucMask = imopen(nucMask,strel('disk',imopenParam));
    end
    
%     figure(1)
%     imshow(nucMask)
    
    
    
    
    
%%    
    
% nucMaskResize=imresize(nucMask, ResizeFactor);
% 
% 
% FigureOverlay=figure;
% cc=1;
% 
% %Plot the image and overlay with the different particles found
% 
% while (cc~=13)
% 
% 
%     SurfMaskResizeZoom=...
%         nucMaskResize(RowsResized/2-RowsZoom/2+ShiftRow*ZoomRatio:RowsResized/2+RowsZoom/2-1+ShiftRow*ZoomRatio,...
%         ColumnsResized/2-ColumnsZoom/2+ShiftColumn*ZoomRatio:ColumnsResized/2+ColumnsZoom/2-1+ShiftColumn*ZoomRatio);
% 
%                     ImOverlay=cat(3,mat2gray(SurfMaskResizeZoom),...
%                         +mat2gray(nucMask),zeros(size(SurfMaskResizeZoom)));
%     imshow(ImOverlay)
%     
%     set(gcf,'name',['ShiftRow: ',num2str(ShiftRow),'. ShiftColumn:',num2str(ShiftColumn),'.'])
%     
%     figure(FigureOverlay)
%     ct=waitforbuttonpress;
%     cc=get(FigureOverlay,'currentcharacter');
%     cm=get(gca,'CurrentPoint');
%     
%     
% 
% 
%     if (ct~=0)&(cc=='.')
%         ShiftColumn=ShiftColumn+1;
%     elseif (ct~=0)&(cc==',')
%         ShiftColumn=ShiftColumn-1;
%     elseif (ct~=0)&(cc=='a')
%         ShiftRow=ShiftRow-1;
%     elseif (ct~=0)&(cc=='z')
%         ShiftRow=ShiftRow+1;
%     end
% end

