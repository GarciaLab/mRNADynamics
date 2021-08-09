function [I_mask, blurSigma] =getEmbryoMaskLive(I,pixelSize_um)

    % NL: rewrote this on 2021-08-07
    max_dist = round(50/pixelSize_um); % maximum distance mask can be from center
    blurSigma = round(5/pixelSize_um);
    % blur
    I_blurred = imgaussfilt(double(I),blurSigma);
    % find threshold
    level = multithresh(I_blurred);
    % generate mask
    I_mask_temp = I_blurred>level;
    
    % find mask that is closest to image center
    stats = regionprops(I_mask_temp,'ConvexArea','Centroid','ConvexHull');
    
    % calculate distances
    CentroidArray = vertcat(stats.Centroid);        
    imageCenter = [size(I,2)/2+1 size(I,1)/2+1];
    distVec = sqrt((CentroidArray(:,1)-imageCenter(1)).^2 + (CentroidArray(:,2)-imageCenter(2)).^2);        
    
    % get convex areas
    AreaVec = vertcat(stats.ConvexArea);
    
    % find which masks (if any) meet our criteria
    elligible_mask_flags = distVec <= max_dist & AreaVec >= (pi*max_dist^2 / 10);    
    
    % if the closest mask is close enough, take it    
    if any(elligible_mask_flags)
        distVec(~elligible_mask_flags) = Inf;
        [~, nearest_id] = min(distVec);
        hull_points = stats(nearest_id).ConvexHull;
        I_mask = poly2mask(hull_points(:,1),hull_points(:,2),size(I,1),size(I,2));
    % otherwise, let's see what happens if we re-threshold, now excluding
    % the original regions. Sometimes bright objects will obscure other
    % structure, so this is worth a shot
    else
        I_blurred2 = I_blurred;
        se = strel('disk',blurSigma);
        I_mask_temp_d = imdilate(I_mask_temp,se);
        I_blurred2(I_mask_temp_d) = median(I_blurred(~I_mask_temp_d));
        level2 = multithresh(I_blurred2);
        I_mask_temp2 = I_blurred2>level2;
        if any(I_mask_temp2(:))
            % find mask that is closest to image center
            stats = regionprops(I_mask_temp2,'ConvexArea','Centroid','ConvexHull');

            % calculate distances
            CentroidArray = vertcat(stats.Centroid);
            imageCenter = [size(I,2)/2 size(I,1)/2];
            distVec = sqrt((CentroidArray(:,1)-imageCenter(1)).^2 + (CentroidArray(:,2)-imageCenter(2)).^2);                

            % get convex areas
            AreaVec = vertcat(stats.ConvexArea);

            % find which masks (if any) meet our criteria
            elligible_mask_flags = distVec <= max_dist & AreaVec >= pi*max_dist^2; 

            % if the closest mask is close enough, take it    
            if any(elligible_mask_flags)
                distVec(~elligible_mask_flags) = Inf;
                [~, nearest_id] = min(distVec);
                % extract convex hull points
                hull_points = stats(nearest_id).ConvexHull;
                I_mask = poly2mask(hull_points(:,1),hull_points(:,2),size(I,1),size(I,2));
            else
                I_mask = true(size(I));
                warning('Failed to calculate embryo mask')
            end
        else
            I_mask = true(size(I));
            warning('Failed to calculate embryo mask')
        end
    end
        
end
