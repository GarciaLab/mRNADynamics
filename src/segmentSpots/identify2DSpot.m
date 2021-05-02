

function Fits = identify2DSpot(particle_index, image, image_label,...
    dog_image, searchRadius, snippet_size, ...
    pixelSize, show_status, graphicsHandles, addition,...
    spot_props, ml_string, currentFrame, spotIndex, zIndex)
%%
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
% Author (contact): Gabriella Martini (martini@berkeley.edu)
% Created: 03/27/2021
% Last Updated: 03/27/2021
%
% Documented by: Gabriella Martini (martini@berkeley.edu)

forced_centroid = spot_props{1} ;
forced_diameter = spot_props{2} ;
snippet_size = uint16(snippet_size);

Fits = [];

ML = false;
if strcmp(ml_string, 'ML')
    ML = true;
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

bwLabelIntensity = sum(sum(image(image_label == particle_index)));


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
        centroid_y = uint16(forced_centroid(2));
        centroid_x = uint16(forced_centroid(1));
    end
    
    %Now, we'll calculate Gaussian fits
    if centroid_y - snippet_size > 1 && centroid_x - snippet_size > 1 && ...
            centroid_y + snippet_size < size(image, 1) && centroid_x + snippet_size < size(image,2)
        
        snippet = image(centroid_y-snippet_size:centroid_y+snippet_size, centroid_x-snippet_size:centroid_x+snippet_size);
        dogsnip = dog_image(centroid_y-snippet_size:centroid_y+snippet_size, centroid_x-snippet_size:centroid_x+snippet_size);


        neighborhood_Size = forced_diameter*2; %nm
        maxThreshold = 2000; %intensity
        widthGuess = forced_diameter/2; %nm
        offsetGuess = 1000; %intensity
        
        
          
        [fits, relative_errors, residuals, confidence_intervals, gaussianIntensity, gaussian, mesh] =  ...
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
        spot_x = fits(2) - double(snippet_size + centroid_x); %final reported spot position
        %         spot_y = fits(4) - snippet_size + centroid_y;
        spot_y = fits(3) - double(snippet_size + centroid_y);
        
        if show_status && ~isempty(graphicsHandles)
            dogAx = graphicsHandles(2);
            ellipse(double(searchRadius/2),double(searchRadius/2),0,double(centroid_x), double(centroid_y),'r',[],dogAx);
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
                    dog_mask(y,x) = 0;
                else
                    maskArea = maskArea+1;
                    if linearOffset
                        snippet_mask(y,x) = snippet_mask(y,x) - off_x*x - off_y*y;
                    end
                end
            end
        end
        
        sigma_x2 = 0;
        sigma_y2 = 0;
        fixedAreaIntensity = sum(sum(snippet_mask)) - (offset*maskArea); %corrected AR 7/13/2018
        
        dogFixedAreaIntensity = sum(dog_mask(:));
        
      
        
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
            Fits.GaussianInfo = {};
            Fits.GaussianInfo.A = single(fits(1));
            Fits.GaussianInfo.x0 = single(fits(2));
            Fits.GaussianInfo.y0 = single(fits(3));
            Fits.GaussianInfo.rho = single(fits(4));
            Fits.GaussianInfo.sigma_x = single(fits(5));
            Fits.GaussianInfo.sigma_y = single(fits(6));
            Fits.GaussianInfo.offset = single(fits(7));
            Fits.GaussianInfo.offset_x = single(fits(8));
            Fits.GaussianInfo.offset_y = single(fits(9));
            Fits.GaussianError.A = single(relative_errors(1));
            Fits.GaussianError.x0 = single(relative_errors(2));
            Fits.GaussianError.y0 = single(relative_errors(3));
            Fits.GaussianError.rho = single(relative_errors(4));
            Fits.GaussianError.sigma_x = single(relative_errors(5));
            Fits.GaussianError.sigma_y = single(relative_errors(6));
            Fits.GaussianError.offset = single(relative_errors(7));
            Fits.GaussianError.offset_x = single(relative_errors(8));
            Fits.GaussianError.offset_y = single(relative_errors(9));
            Fits.GaussianResiduals = residuals;
            Fits.GaussianFitValues = gaussian;
            Fits.CentralIntensity = single(intensity);
            Fits.DOGIntensity = single(max_dog);
            Fits.ConfidenceIntervals = confidence_intervals;
            Fits.gaussParams = {fits};
            Fits.dogFixedAreaIntensity = single(dogFixedAreaIntensity);
            Fits.intArea = uint16(maskArea);
            Fits.z = uint8(zIndex);
            Fits.bwIntensity = bwLabelIntensity;
            Fits.bwArea = spot_props{3};
            Fits.bwDiameter = forced_diameter;
            Fits.bwCircularity = spot_props{4};
            Fits.bwEccentricity = spot_props{5};
            Fits.bwMajorAxisLength = spot_props{6};
            Fits.bwMinorAxisLength = spot_props{7};
            Fits.bwPixelList= spot_props{8};
            Fits.discardThis = false;
            Fits.r = false;
            Fits.snippet_size = uint8(snippet_size);
            Fits.Approved = 0; 
            
            
            
            

        else
            Fits = [];
        end
        
    end
    
end
end
