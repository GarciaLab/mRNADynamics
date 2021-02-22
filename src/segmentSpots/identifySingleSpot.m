

function Fits = identifySingleSpot(particle_index, image, image_label,...
    dog_image, searchRadius, snippet_size, ...
    pixelSize, show_status, graphicsHandles, microscope, addition,...
    forced_centroid, ml_string, currentFrame, spotIndex, zIndex)
% identifySingleSpot(awholelot)
%
% DESCRIPTION
% This is a subfunction for 'segmentSpots' that locates a transcriptional locus (the k'th locus in an image)
% and assigns a Gaussian to it. In a future release, this will be replaced
% by identifySpot.
%
% ARGUMENTS
% This is a massive list that, honestly, just needs to be refactored. (AR)
%
% OUTPUT
% tFits:  Returns a structure containing properties of the
%                 transcriptional locus such as the intensity, size, position, etc.
%
% Author (contact): Armando Reimer (areimer@berkeley.edu)
% Created: 01/01/2016
% Last Updated: 12/31/2016
%
% Documented by: Armando Reimer (areimer@berkeley.edu)


Fits = [];

ML = false;
if strcmp(ml_string, 'ML')
    ML = true;
end

doCyl = false;
if iscell(image)
    doCyl = true;
    imageAbove = image{2};
    imageBelow = image{3};
    image = image{1};
end

%Find spot centroids in the actual image by hunting for global maxima in
%neighborhoods around spots that were just located

possible_centroid_intensity = [];
possible_centroid_location = {};
if ~addition(1) %this gets flagged if we're not manually adding a particle in checkparticletracking
    [k_row, k_column] = find(image_label == particle_index); %the position of the k'th locus in the image
    row = k_row(1);
    col = k_column(1);
else %if we're in checkparticletracking
    row = addition(3);
    col = addition(2);
    k_row = row;
    k_column = col;
end


for y = 1:2*searchRadius
    for x = 1:2*searchRadius
        if row - searchRadius + y > 0 && col - searchRadius + x > 0 ...
                && row - searchRadius + y < size(image,1)  && col - searchRadius + x < size(image,2)
            if ML
                possible_centroid_intensity(y,x) = sum(sum(image(row-searchRadius+y, col-searchRadius+x)));
            else
                if addition(1) || round(pixelSize)~=212
                                possible_centroid_intensity(y,x) = sum(sum(image(row-searchRadius+y, col-searchRadius+x)));
                else
                    possible_centroid_intensity(y,x) = sum(sum(image(row-searchRadius+y, col-searchRadius+x)));
                end
            end
            possible_centroid_location{y,x} = [row-searchRadius+y, col-searchRadius+x];
        end
    end
end
clear row;
clear col;

%Compute some preliminary properties of the located spots

if ~isempty(possible_centroid_intensity) && sum(sum(possible_centroid_intensity)) ~= 0
    
    [intensity, centroid_index] = max(possible_centroid_intensity(:));
    [row, col] = ind2sub(size(possible_centroid_intensity),centroid_index);
    centroid_y = possible_centroid_location{row,col}(1);
    centroid_x = possible_centroid_location{row,col}(2);
    
    %if the predicted centroid is outside of the image while manually
    %adding spots in checkparticletracking, we'll try to use the
    %position where the user actually clicked
    if addition(1) && ~(centroid_y - snippet_size > 1 && centroid_x - snippet_size > 1 ...
            && centroid_y + snippet_size < size(image, 1) && centroid_x + snippet_size < size(image,2))
        centroid_y = k_row;
        centroid_x = k_column;
    end
    
    if ~isempty(forced_centroid)
        centroid_y = forced_centroid(2);
        centroid_x = forced_centroid(1);
    end
    
    %Now, we'll calculate Gaussian fits
    if centroid_y - snippet_size > 1 && centroid_x - snippet_size > 1 && ...
            centroid_y + snippet_size < size(image, 1) && centroid_x + snippet_size < size(image,2)
        
        snippet = image(centroid_y-snippet_size:centroid_y+snippet_size, centroid_x-snippet_size:centroid_x+snippet_size);
        dogsnip = dog_image(centroid_y-snippet_size:centroid_y+snippet_size, centroid_x-snippet_size:centroid_x+snippet_size);
        
        if doCyl
            snippetAbove = imageAbove(centroid_y-snippet_size:centroid_y+snippet_size, centroid_x-snippet_size:centroid_x+snippet_size);
            snippetBelow = imageBelow(centroid_y-snippet_size:centroid_y+snippet_size, centroid_x-snippet_size:centroid_x+snippet_size);
        end
        
        % Set parameters to use as initial guess in the fits.
        if strcmp(microscope, 'LAT')
            neighborhood_Size = 2000/pixelSize; %nm
            maxThreshold = 2000; %intensity
            widthGuess = 500 / pixelSize; %nm
            offsetGuess = 1000; %intensity
            
            % Confocal Leica SP8 Data
        else
            neighborhood_Size = 2000/pixelSize; %nm
            maxThreshold = 30; %intensity
            widthGuess = 200 / pixelSize; %nm
            offsetGuess = median(snippet(:)); %intensity
        end
          
        [fits, relative_errors, ~, confidence_intervals, gaussianIntensity, gaussian, mesh] =  ...
            fitSingleGaussian(snippet, neighborhood_Size, maxThreshold, ...
            widthGuess, offsetGuess, show_status, graphicsHandles);
        %fits: [amplitude, x position, x width, y position, y width, offset, angle]
        
        % @(A, x0, y0, rho, sigma_x, sigma_y, offset, offset_x, offset_y)
        
        %         sigma_x = fits(3);
        %         sigma_y = fits(5);
        %         offset = fits(6);
        
        sigma_x = fits(5);
        sigma_y = fits(6);
        offset = fits(7);
        
        
        gaussianArea = pi*sigma_x*sigma_y; %in pixels. this is one width away from peak
        integration_radius = 6*ceil(sqrt(212/pixelSize)); %integrate 109 pixels around the spot with 212nm pixel size
        spot_x = fits(2) - snippet_size + centroid_x; %final reported spot position
        %         spot_y = fits(4) - snippet_size + centroid_y;
        spot_y = fits(3) - snippet_size + centroid_y;
        
        if show_status && ~isempty(graphicsHandles)
            dogAx = graphicsHandles(2);
            ellipse(searchRadius/2,searchRadius/2,0,centroid_x, centroid_y,'r',[],dogAx);
            drawnow;
        end
        
        %disp(rel_errors);
        % Quality control.
        % TODO: make some quality control using the errors in
        % the fits. Using any(rel_errors > 0.3) (for
        % example) is not the best thing to do, because
        % sometimes the second gaussian doesn't get a good fit
        % but the first one does, and the second one is good
        % enough to position its center.
        
        snippet_mask = double(snippet);
        dog_mask = double(dogsnip);
        if doCyl
            snippet_mask_above = double(snippetAbove);
            snippet_mask_below = double(snippetBelow);
        end
        
        if length(fits)>7
            linearOffset = true;
            off_x = fits(8);
            off_y = fits(9);
        end
        
        maskArea = 0;
        for y = 1:size(snippet, 1)
            for x = 1:size(snippet,2)
                d = sqrt( (x - (size(snippet,1)+1)/2)^2 + (y - (size(snippet,2)+1)/2)^2) ;
                if d >= integration_radius
                    snippet_mask(y, x) = 0;
                    snippet_mask_above(y,x) = 0;
                    snippet_mask_below(y,x) = 0;
                    dog_mask(y,x) = 0;
                else
                    maskArea = maskArea+1;
                    if linearOffset
                        snippet_mask(y,x) = snippet_mask(y,x) - off_x*x - off_y*y;
                        snippet_mask_above(y,x) = snippet_mask_above(y,x) - off_x*x - off_y*y;
                        snippet_mask_below(y,x) = snippet_mask_below(y,x) - off_x*x - off_y*y;
                    end
                end
            end
        end
        
        sigma_x2 = 0;
        sigma_y2 = 0;
        fixedAreaIntensity = sum(sum(snippet_mask)) - (offset*maskArea); %corrected AR 7/13/2018
        
        dogFixedAreaIntensity = sum(dog_mask(:));
        
        if doCyl
            fixedAreaIntensityCyl3 =  sum(sum(snippet_mask)) + sum(sum(snippet_mask_above))...
                + sum(sum(snippet_mask_below)) - 3*offset*maskArea;
        end
        
        if  .1<sigma_x && sigma_x<(600/pixelSize) && .1<sigma_y && sigma_y<(600/pixelSize)...
                || addition(1) %here is the place to introduce quality control
            try
                max_dog = max(max(dog_image(k_row,k_column)));
            catch exceptionMaxDOG
                exceptionMessage = exceptionMaxDOG.message;
                
                if contains(exceptionMessage, 'array exceeds maximum array size preference')
                    fprintf('***Error Suggestions: Please read!***');
                    fprintf('\nYou have requested memory for an array that exceeds maximum array size preferences set by MATLAB.');
                    fprintf('\nThis could be due to one of the following reasons:');
                    fprintf('\n(1) Your threshold is set too low, thus finding too many false positives. Try a slightly higher threshold.');
                    fprintf('\n(2) Your ML classifier is not good enough and is finding too many false positives. Go back and retrain your classifier to find fewer false positives.')
                    fprintf('\n(3) In Weka, you made class 1 (red) "not spots" and class 2 (green) "spots". Train a new classifier where class 1 is "spots" and class 2 is "not spots".\n');
                    rethrow(exceptionMaxDOG);
                else
                    rethrow(exceptionMaxDOG);
                end
            end
                        
            Fits.FixedAreaIntensity = single(fixedAreaIntensity);
            Fits.xFit = single(spot_x);
            Fits.yFit = single(spot_y);
            Fits.Offset = single(offset);
            Fits.xFitWidth = single(sigma_x);
            Fits.yFitWidth = single(sigma_y);
            Fits.yDoG = uint16(centroid_y);
            Fits.xDoG = uint16(centroid_x);
            Fits.GaussianIntensity = single(gaussianIntensity);
            Fits.CentralIntensity = single(intensity);
            Fits.DOGIntensity = single(max_dog);
            Fits.ConfidenceIntervals = confidence_intervals;
            Fits.gaussParams = {fits};
            Fits.dogFixedAreaIntensity = single(dogFixedAreaIntensity);
            Fits.intArea = uint16(maskArea);
            Fits.z = uint8(zIndex);
            Fits.discardThis = false;
            Fits.r = false;
            Fits.FixedAreaIntensity3  = [];
            Fits.brightestZ =[];
            Fits.snippet_size = uint8(snippet_size);
            Fits.Approved = 0; 
        else
            Fits = [];
        end
        
    end
    
end
end
