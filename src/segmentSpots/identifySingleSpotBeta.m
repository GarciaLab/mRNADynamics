function temp_particles = identifySingleSpotBeta(particle_index, image, image_label, dog_image, searchRadius, snippet_size, ...
    pixelSize, show_status, fg, microscope, addition, forced_centroid, ml_string, intScale)
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
% tempParticles:  Returns a structure containing properties of the
%                 transcriptional locus such as the intensity, size, position, etc.
%
% Author (contact): Armando Reimer (areimer@berkeley.edu)
% Created: 01/01/2016
% Last Updated: 12/31/2016
%
% Documented by: Armando Reimer (areimer@berkeley.edu)

    ML = 0;
    if strcmp(ml_string, 'ML')
        ML = 1;
    end

    doCyl = 0;
    if iscell(image)
        doCyl = 1;
        imageAbove = image{2};
        imageBelow = image{3};
        image = image{1};
    end
    
    %Find spot centroids in the actual image by hunting for global maxima in
    %neighborhoods around spots that were just located
    
    possible_centroid_intensity = [];
    possible_centroid_location = {};
    if addition(1) %this gets flagged if we're not manually adding a particle in checkparticletracking
        centroid_y = addition(3);
        centroid_x = addition(2);
        intensity = image(centroid_y,centroid_x);
    elseif ~isempty(forced_centroid) %if we're in checkparticletracking 
        centroid_y = forced_centroid(2);
        centroid_x = forced_centroid(1);
        intensity = image(centroid_y,centroid_x);       
    else
        [k_row, k_column] = find(image_label == particle_index); %the position of the k'th locus in the image
        mask = image;
        mask(image_label~=particle_index) = nan;
        [intensity, mi] = nanmax(mask(:));
        [centroid_y, centroid_x] = ind2sub(size(mask),mi);
        row = k_row(1);
        col = k_column(1);
    end    
    %Compute some preliminary properties of the located spots
    temp_particles = {[]};
    
    
    tp = struct('fixedAreaIntensity', [], 'xFit', [], 'yFit', [], 'Offset', [],...
       'Snippet', [], 'Area', [], 'xFitWidth', [], 'yFitWidth', [], 'centroidY',...
       [], 'centroidX', [], 'GaussianIntensity', [], 'CentralIntensity', [],...
       'DOGIntensity', [], 'snippet_mask', [], 'sistersXWidth', [], 'sistersYWidth',...
       [], 'sisterSeparation', [], 'relative_errors', [], 'ConfidenceIntervals',...
       [], 'gaussSpot', [], 'mesh', [], 'sistersXPos', [], 'sistersYPos', [],...
       'eccentricity', [],'sisterSeparation2', [],'sistersXPos2', [], 'sistersYPos2', [],...
       'sistersXWidth2', [], 'sistersYWidth2',[], 'sisterAmps', []);
   
   %AR 7192018- i want to switch the tempparticles over to this structure
   %in the future for reading clarity. currently doesn't have the right
   %fields and is unused. 
      
    if ~isnan(intensity)
        
       %Now, we'll calculate Gaussian fits 
       if centroid_y - snippet_size > 1 && centroid_x - snippet_size > 1 && ...
               centroid_y + snippet_size < size(image, 1) && centroid_x + snippet_size < size(image,2)
           
            snippet = image(centroid_y-snippet_size:centroid_y+snippet_size, centroid_x-snippet_size:centroid_x+snippet_size);
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
                offsetGuess = 0; %intensity
            end
            
            chisq = [];
            fitstruct = struct();
            %steps = round(1/pixelSize):round(50/pixelSize):round(400/pixelSize);
            steps = 200/pixelSize:200/pixelSize;
            for i = 1:length(steps)
                [fits, relative_errors, residual, confidence_intervals, gaussianIntensity, gaussian, mesh] =  ...
                    fitSingleGaussian(snippet, neighborhood_Size, maxThreshold, ...
                    steps(i), offsetGuess, show_status);
                fitstruct(i).fits = fits;
                fitstruct(i).relative_errors = relative_errors;
                fitstruct(i).residual = residual;
                fitstruct(i).confidence_intervals = confidence_intervals;
                fitstruct(i).gaussianIntensity = gaussianIntensity;
                fitstruct(i).gaussian = gaussian;
                fitstruct(i).mesh = mesh;
                chisq(i) = sum(sum(residual.^2));
            end     
            [~, optindex] = min(chisq);
            fits = fitstruct(optindex).fits;
            relative_errors = fitstruct(optindex).relative_errors;
            residual = fitstruct(optindex).residual;
            confidence_intervals = fitstruct(optindex).confidence_intervals;
            gaussianIntensity = fitstruct(optindex).gaussianIntensity;
            mesh = fitstruct(optindex).mesh;
            
            
            sigma_x = fits(3);
            sigma_y = fits(5);

            gaussianArea = pi*sigma_x*sigma_y; %in pixels. this is one width away from peak
            integration_radius = 6*intScale; %integrate 109 pixels around the spot or more optionally
            spot_x = fits(2) - snippet_size + centroid_x; %final reported spot position
            spot_y = fits(4) - snippet_size + centroid_y;    
            
            if show_status && ~isempty(fg)
                figure(fg)
                ellipse(searchRadius/2,searchRadius/2,0,centroid_x, centroid_y,'r',[],gca);
                pause(.1) %Ellipses won't be plotted correctly without this pause.
                %figure(5)
                %imshow(image,[])
                %ellipse(searchRadius/2,searchRadius/2,0,spot_x,spot_y,'r');
                %pause(.1) %Ellipses won't be plotted correctly without this pause.
                
            end
            
            %disp(rel_errors);
            % Quality control.
            % TODO: make some quality control using the errors in
            % the fits. Using any(rel_errors > 0.3) (for
            % example) is not the best thing to do, because
            % sometimes the second gaussian doesn't get a good fit
            % but the first one does, and the second one is good
            % enough to position its center.
            
            snippet_mask = snippet;
            if doCyl
                snippet_mask_above = snippetAbove;
                snippet_mask_below = snippetBelow;
            end
            maskArea = 0;
            tic
            for i = 1:size(snippet, 1)
                for j = 1:size(snippet,2)
                    d = sqrt( (j - (size(snippet,1)+1)/2)^2 + (i - (size(snippet,2)+1)/2)^2) ;
                    if d >= integration_radius
                        snippet_mask(i, j) = 0;
                        snippet_mask_above(i,j) = 0;
                        snippet_mask_below(i,j) = 0;
                    else
                        maskArea = maskArea+1;
                    end 
                end
            end
            toc
            sigma_x2 = 0;
            sigma_y2 = 0;
            sister_chromatid_distance = fits(end);
            fixedAreaIntensity = sum(sum(snippet_mask)) - fits(end-1)*maskArea; %corrected AR 7/13/2018
            fixedAreaIntensityCyl3 = NaN;
            if doCyl
                fixedAreaIntensityCyl3 =  sum(sum(snippet_mask)) + sum(sum(snippet_mask_above))...
                    + sum(sum(snippet_mask_below)) - 3*fits(end-1)*maskArea;
            end
            
            if  .1<sigma_x && sigma_x<(600/pixelSize) && .1<sigma_y && sigma_y<(600/pixelSize)...
                    || addition(1) %here is the place to introduce quality control
                try
                    max_dog = max(dog_image(image_label == particle_index));
                catch exceptionMaxDOG
                    exceptionMessage = exceptionMaxDOG.message;
                
                    if contains(exceptionMessage, 'array exceeds maximum array size preference')
                        fprintf('***Error Suggestions: Please read!***');
                        fprintf('\nYou have requested memory for an array that exceeds maximum array size preferences set by MATLAB.');
                        fprintf('\nThis could be due to one of the following reasons:');
                        fprintf('\n(1) Your ML classifier is not good enough and is finding too many false positives. Go back and retrain your classifier to find less false positives (i.e. Do better next time.).')
                        fprintf('\n(2) In Weka, you made class 1 (red) "not spots" and class 2 (green) "spots". Train a new classifier where class 1 is "spots" and class 2 is "not spots".\n');
                        rethrow(exceptionMaxDOG);
                    else
                        rethrow(exceptionMaxDOG);
                    end
                end
                temp_particles = {{fixedAreaIntensity, spot_x, spot_y, fits(end-1), snippet, ...
                    gaussianArea, sigma_x, sigma_y, centroid_y, centroid_x, gaussianIntensity,intensity,...
                    max_dog, snippet_mask, sigma_x2, sigma_y2, sister_chromatid_distance, relative_errors, confidence_intervals, gaussian, mesh,fits, maskArea, fixedAreaIntensityCyl3}};
            else                
                temp_particles = {{}};   
            end

       end
    end
end
