function [nucMask]=getNuclearMaskTile(I, nucleusDiameter, PixelSize, ImposeEmbryoMask)

    f_sigma = round(nucleusDiameter /(6*PixelSize));
    imopenParam= round(f_sigma/4);
    if ~ImposeEmbryoMask


        I_blurred = imfilter(I,...
             fspecial('gaussian',2*f_sigma,f_sigma),'symmetric','conv');
        I_mask = imbinarize(I_blurred);
    else
        I_mask = ones(size(I));
    end
    segmentationOk = false;
    %smoothing,imopenParam These are the second and third arguments to the
    %GetNuclearMask function
    while ~segmentationOk
        if f_sigma > 0
           I_smooth = imfilter(I,...
                    fspecial('gaussian',max(f_sigma*2, 1), f_sigma),'symmetric'); 
        else
            I_smooth = I;
        end
        I_segmented = single(watershed(imcomplement(I_smooth))).*single(I_mask);
        imshow(I_segmented, [])
        cc=bwconncomp(I_segmented);
        if cc.NumObjects<5000
            segmentationOk = true;
        else 
            f_sigma = f_sigma+5;
        end
    end

        



%% 




    
    P=regionprops(I_segmented,'Area','PixelIdxList');
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
    
    I_segmented(I_segmented==0)=length(areas)+1;
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
    
    nucMask = (I > threshList(I_segmented));
    
    % fill holes
    nucMask = imfill(nucMask,'holes');
    
    if imopenParam>0
        % remove small specks    
        nucMask = imopen(nucMask,strel('disk',imopenParam));
    end
    
end