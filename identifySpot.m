function tp = identifySpot(particle_index, image, image_label, dog_image, distance_to_neighbor, snippet_size, ...
    pixelSize, show_status, fg, microscope, addition)
% identifySpot(awholelot)
%
% DESCRIPTION
% This is a subfunction for 'segmentSpots' that locates a transcriptional locus (the k'th locus in an image)
% and assigns a Gaussian to it.
% This will likely replace identifySingleSpot in a future release.
%
% ARGUMENTS
% This is a massive list that, honestly, just needs to be refactored. (AR)
%
% OUTPUT
% tp:  Returns a structure containing properties of the
%                 transcriptional locus such as the intensity, size, position, etc.
%
% Author (contact): Armando Reimer (areimer@berkeley.edu)
% Created: 01/01/2016
% Last Updated: 12/31/2016
%
% Documented by: Armando Reimer (areimer@berkeley.edu)

    %Find spot centroids in the actual image by hunting for global maxima in
    %neighborhoods around spots that were just located
   tp = struct('fixedAreaIntensity', [], 'xFit', [], 'yFit', [], 'Offset', [],...
       'Snippet', [], 'Area', [], 'xFitWidth', [], 'yFitWidth', [], 'centroidY',...
       [], 'centroidX', [], 'GaussianIntensity', [], 'CentralIntensity', [],...
       'DOGIntensity', [], 'snippet_mask', [], 'sistersXWidth', [], 'sistersYWidth',...
       [], 'sisterSeparation', [], 'relative_errors', [], 'ConfidenceIntervals',...
       [], 'gaussSpot', [], 'mesh', [], 'sistersXPos', [], 'sistersYPos', [],...
       'eccentricity', [],'sisterSeparation2', [],'sistersXPos2', [], 'sistersYPos2', [],...
       'sistersXWidth2', [], 'sistersYWidth2',[], 'sisterAmps', []);
    possible_centroid = [];
    possible_centroid_location = {};
    if ~addition(1)
        [k_row, k_column] = find(image_label == particle_index); %the position of the k'th locus in the image
        row = k_row(1);
        col = k_column(1);
    else 
        row = addition(3);
        col = addition(2);
        k_row = row;
        k_column = col;
        
    end
    for i = 1:2*distance_to_neighbor
        for j = 1:2*distance_to_neighbor
            if row - distance_to_neighbor + i > 0 && col - distance_to_neighbor + j > 0 ... 
                    && row - distance_to_neighbor + i < size(image,1)  && col - distance_to_neighbor + j < size(image,2)
                possible_centroid(i,j) = image(row-distance_to_neighbor+i, col-distance_to_neighbor+j);
                possible_centroid_location{i,j} = [row-distance_to_neighbor+i, col-distance_to_neighbor+j];
            end
        end
    end
    clear row;
    clear col;
    
    %Compute some preliminary properties of the located spots
    if ~isempty(possible_centroid)
        [intensity, centroid_index] = max(possible_centroid(:));
        [row, col] = ind2sub(size(possible_centroid),centroid_index);
        centroid_y = possible_centroid_location{row,col}(1); 
        centroid_x = possible_centroid_location{row,col}(2);
       
        if addition(1) && ~(centroid_y - snippet_size > 1 && centroid_x - snippet_size > 1 && centroid_y + snippet_size < size(image, 1) && centroid_x + snippet_size < size(image,2))    
            centroid_y = k_row;
            centroid_x = k_column;
        end
        
       %Now, we'll calculate Gaussian fits 
       if centroid_y - snippet_size > 1 && centroid_x - snippet_size > 1 && centroid_y + snippet_size < size(image, 1) && centroid_x + snippet_size < size(image,2)
           
            snippet = image(centroid_y-snippet_size:centroid_y+snippet_size, centroid_x-snippet_size:centroid_x+snippet_size);

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
                    fitGaussians(snippet, neighborhood_Size, maxThreshold, ...
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
            
            if .1<fits(1).sigmaX && fits(1).sigmaX<(600/pixelSize) && .1<fits(1).sigmaY && fits(1).sigmaY<(600/pixelSize)...
                || addition(1) %some quality control based on widths
            
                area = pi*fits(1).sigmaX*fits(1).sigmaY; %in pixels. this is one width away from peak
                integration_radius = 6; %integrate 109 pixels around the spot
                %Shift the positions to be relative to the larger image
                fits(1).XPosition = fits(1).XPosition - snippet_size + centroid_x; 
                fits(1).YPosition = fits(1).YPosition - snippet_size + centroid_y;   
                fits(2).XPosition = fits(2).XPosition - snippet_size + centroid_x; 
                fits(2).YPosition = fits(2).YPosition - snippet_size + centroid_y;   
                fits(3).XPosition = fits(3).XPosition - snippet_size + centroid_x; 
                fits(3).YPosition = fits(3).YPosition - snippet_size + centroid_y;
                fits(4).XPosition = fits(4).XPosition - snippet_size + centroid_x; 
                fits(4).YPosition = fits(4).YPosition - snippet_size + centroid_y;
                fits(5).XPosition = fits(5).XPosition - snippet_size + centroid_x; 
                fits(5).YPosition = fits(5).YPosition - snippet_size + centroid_y; 

                if show_status && ~isempty(fg)
                    figure(fg)
                    ellipse(distance_to_neighbor/2,distance_to_neighbor/2,0,fits(1).XPosition,fits(1).YPosition,'r');
                    pause(.1) %Ellipses won't be plotted correctly without this pause.
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
                for i = 1:size(snippet, 1)
                    for j = 1:size(snippet,2)
                        d = sqrt( (j - (size(snippet,1)+1)/2)^2 + (i - (size(snippet,2)+1)/2)^2) ;
                        if d >= integration_radius
                            snippet_mask(i, j) = 0;
                        end 
                    end
                end
                if fits(1).sigmaX > fits(1).sigmaY
                    semimaj = 2*fits(1).sigmaX;
                    semimin = 2*fits(1).sigmaY;
                else
                    semimaj = 2*fits(1).sigmaY;
                    semimin = 2*fits(1).sigmaX;
                end
                ecc = sqrt(1 - (semimin^2/semimaj^2)); %this measure probably isn't indicative of sister separation
                sister_separation = sqrt((fits(2).XPosition-fits(3).XPosition)^2 + (fits(2).YPosition-fits(3).YPosition)^2); % in pixels
                fixedAreaIntensity = sum(sum(snippet_mask)) - fits(1).Offset*sum(sum(snippet_mask~=0));
                max_dog = max(max(dog_image(k_row,k_column)));
                %return a temporary particle structure that will be
                %reshaped by segmentSpots
                tp.fixedAreaIntensity = fixedAreaIntensity;
                tp.xFit = fits(2).XPosition;
                tp.yFit = fits(2).YPosition;
                tp.Offset = fits(1).Offset;
                tp.Snippet = snippet;
                tp.Area = area;
                tp.xFitWidth = fits(1).sigmaX;
                tp.yFitWidth = fits(1).sigmaY;
                tp.centroidY = centroid_y;
                tp.centroidX = centroid_x;
                tp.GaussianIntensity = gaussianIntensity;
                tp.CentralIntensity = intensity;
                tp.DOGIntensity = max_dog;
                tp.snippet_mask = snippet_mask;
                tp.sistersXWidth = [fits(2).sigmaX, fits(3).sigmaX];
                tp.sistersYWidth = [fits(2).sigmaY, fits(3).sigmaY];
                tp.sisterSeparation = sister_separation;
                tp.relative_errors = relative_errors;
                tp.ConfidenceIntervals = confidence_intervals;
                tp.gaussSpot = gaussian;
                tp.mesh = mesh;
                tp.sistersXPos = [fits(2).XPosition,fits(3).XPosition];
                tp.sistersYPos = [fits(2).YPosition,fits(3).YPosition];
                tp.eccentricity = ecc;
                tp.sistersXPos2 = [fits(4).XPosition,fits(5).XPosition];
                tp.sistersYPos2 = [fits(4).YPosition,fits(5).YPosition];
                tp.sisterSeparation2 = sqrt( (fits(4).XPosition-fits(5).XPosition).^2 ...
                + (fits(4).YPosition-fits(5).YPosition).^2 );
                tp.sistersXWidth2 = [fits(4).sigmaX, fits(5).sigmaX];
                tp.sistersYWidth2 = [fits(4).sigmaY, fits(5).sigmaY];
                tp.sisterAmps = [fits(4).Amplitude, fits(5).Amplitude];               
            end
        end
    end
end
