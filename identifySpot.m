function temp_particles = identifySpot(k, im, im_label, dog, neighb, rad, ...
    pixelSize, show_status, f,microscope)

    [r,c] = find(im_label == k);

    % Keep track of the maximum dog value

    max_dog = max(max(dog(r,c)));

    %Find spot centroids in the actual image by hunting for global maxima in
    %neighborhoods around spots that were just located

    possible_cent = [];
    pcentloc = {};

    for o = 1:2*neighb
        for p = 1:2*neighb
            if r(1) - neighb + o > 0 && c(1) - neighb + p > 0 ... 
                    && r(1) - neighb + o < size(im,1)  && c(1) - neighb + p < size(im,2)
                possible_cent(o,p) = im(r(1)-neighb+o, c(1)-neighb+p);
                pcentloc{o,p} = [r(1)-neighb+o, c(1)-neighb+p];
            end
        end
    end
    if ~isempty(possible_cent)
        [inten, index] = max(possible_cent(:));
        [row, col] = ind2sub(size(possible_cent),index);
        cent_y = pcentloc{row,col}(1); 
        cent_x = pcentloc{row,col}(2);
       % temp_particles = [temp_particles,[0, cent_x, cent_y, 0, 0]];
    %    temp_particles = {};
       if show_status && ~isempty(f)
            set(0,'CurrentFigure', f);...
            ellipse(neighb/2,neighb/2,0,cent_x,cent_y,'r');
       end

       if cent_y - rad > 1 && cent_x - rad > 1 && cent_y + rad < size(im, 1) && cent_x + rad < size(im,2)
           snip = im(cent_y-rad:cent_y+rad, cent_x-rad:cent_x+rad);

            % Set parameters to use as initial guess in the fits. For the 
            % lattice data, try NeighborhoodSize = 1000, MaxThreshold = 2000, 
            % WidthGuess = 500, OffsetGuess = 1000.

            if strcmp(microscope, 'LAT')

                NeighborhoodSize = 1000/pixelSize; %nm
                MaxThreshold = 2000; %intensity
                WidthGuess = 500 / pixelSize; %nm
                OffsetGuess = 1000; %intensity

                % For confocal data, try NeighborhoodSize = 1000, MaxThreshold = 20,
                % WidthGuess = 200, OffsetGuess = 10.
            else 
                NeighborhoodSize = 1000/pixelSize; %nm
                MaxThreshold = 30; %intensity
                WidthGuess = 200 / pixelSize; %nm
                OffsetGuess = 10; %intensity
            end

            [fits, rel_errors, ci, GaussianIntensity, f2, f4, AmpIntegral] =  ...
                fitGaussians(snip, NeighborhoodSize, MaxThreshold, ...
                WidthGuess, OffsetGuess, show_status);
            sigma_x = fits(3);
            sigma_y = fits(3);
            sigma_x2 = fits(7);
            sigma_y2 = fits(7);
            area = pi*(2*sigma_x^2)^2; %in pixels. this is two widths away from peak
            fixedAreaIntensity = 0;
            integration_radius = 5; %integrate 121 pixels around the spot
            c_x = fits(2) - rad + cent_x; %AR 7/14/16: same deal as line above
            c_y = fits(4) - rad + cent_y;
            int_x = [round(c_x - integration_radius), round(c_x + integration_radius)];
            int_y = [round(c_y - integration_radius), round(c_y + integration_radius)];    
    %         int_x = [round(c_x - fits(3)), round(c_x + fits(3))];
    %         int_y = [round(c_y - fits(5)), round(c_y + fits(5))];

            %disp(rel_errors);
            % Quality control.
            % TODO: make some quality control using the errors in
            % the fits. Using any(rel_errors > 0.3) (for
            % example) is not the best thing to do, because
            % sometimes the second gaussian doesn't get a good fit
            % but the first one does, and the second one is good
            % enough to position its center.
            snip_mask = snip*0;
            for i = 1:size(snip,1)
                for j = 1:size(snip,2)
                    if i > fits(4) - integration_radius && i < fits(4) + integration_radius && j > fits(2) - integration_radius && j < fits(2) + integration_radius 
                        snip_mask(i,j) = 1;
                    end
                end
            end

            if ~(sigma_x2 <= 0 || sigma_x <= 0 || sigma_x > 2000/pixelSize || sigma_y > 2000/pixelSize...
                    || sigma_x2 > 2000/pixelSize || sigma_y2 > 2000/pixelSize...
                    || GaussianIntensity == 0)

                if int_x(1) > 1 && int_y(1) > 1 && int_x(2) < size(im,2) && int_y(2) < size(im,1)
                    for w = int_x(1):int_x(2)
                        for v = int_y(1): int_y(2)
                            fixedAreaIntensity = fixedAreaIntensity + double(im(v,w)) - fits(end-1);
                        end
                    end
                end
                temp_particles = {{fixedAreaIntensity, c_x, c_y, fits(end-1), snip, ...
                    area, sigma_x, sigma_y, cent_y, cent_x, GaussianIntensity,inten,...
                    max_dog, snip_mask, sigma_x2, sigma_y2, fits(end), rel_errors, ci, f2, f4}};
            else
                temp_particles = {[]};
            end
       else 
           temp_particles = {[]};
       end
    else 
        temp_particles = {[]};
    end
end
