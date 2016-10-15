function temp_particles = identifySpot(particle_index, image, image_label, dog_image, distance_to_neighbor, snippet_size, ...
    pixelSize, show_status, figure, microscope, threshold)

    %This function locates a transcriptional locus (the k'th locus in an image)
    %and assigns a Gaussian
    %to it. It returns a structure containing properties of the
    %transcriptional locus such as the intensity, size, position, etc.
    
    %Arguments: This requires an image, its difference of gaussians image,
    %as well as some properties of the recording. 

    
    
    %Find spot centroids in the actual image by hunting for global maxima in
    %neighborhoods around spots that were just located

    possible_centroid = [];
    possible_centroid_location = {};
    [k_row, k_column] = find(image_label == particle_index); %the position of the k'th locus in the image
    row = k_row(1);
    col = k_column(1);
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
    temp_particles = {[]};
    if ~isempty(possible_centroid)
        [intensity, centroid_index] = max(possible_centroid(:));
        [row, col] = ind2sub(size(possible_centroid),centroid_index);
        centroid_y = possible_centroid_location{row,col}(1); 
        centroid_x = possible_centroid_location{row,col}(2);
       
        if show_status && ~isempty(figure)
            set(0,'CurrentFigure', figure);...
            ellipse(distance_to_neighbor/2,distance_to_neighbor/2,0,centroid_x,centroid_y,'r');
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
            steps = 1:50:400;
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
            
            
            sigma_x = fits(3);
            sigma_y = fits(3);
            sigma_x2 = fits(7);
            sigma_y2 = fits(7);
            area = pi*(2*sigma_x^2)^2; %in pixels. this is two widths away from peak
            integration_radius = 6; %integrate 109 pixels around the spot
            spot_x = fits(2) - snippet_size + centroid_x; %AR 7/14/16: same deal as line above
            spot_y = fits(4) - snippet_size + centroid_y;    

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
                                  
            if 1   %here is the place to introduce quality control
                fixedAreaIntensity = sum(sum(snippet_mask)) - fits(end-1)*sum(sum(snippet_mask~=0));
                max_dog = max(max(dog_image(k_row,k_column)));
                temp_particles = {{fixedAreaIntensity, spot_x, spot_y, fits(end-1), snippet, ...
                    area, sigma_x, sigma_y, centroid_y, centroid_x, gaussianIntensity,intensity,...
                    max_dog, snippet_mask, sigma_x2, sigma_y2, fits(end), relative_errors, confidence_intervals, gaussian, mesh}};
            end
       end
    end
end
