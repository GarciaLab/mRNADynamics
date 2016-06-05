function temp_particles = fit_single_spot(k, im, im_label, dog, neighb, rad, ...
    pixelSize, show_status, f)

[r,c] = find(im_label == k);

% Keep track of the maximum dog value

max_dog = max(max(dog(r,c)));

%Find spot centroids in the actual image by hunting for absolute maxima in
%neighborhoods around spots that were just located

possible_cent = [];
pcentloc = {};
cent = [];
cent_intensity = 0;
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
   if show_status
        set(0,'CurrentFigure', f);
        ellipse(neighb/2,neighb/2,0,cent_x,cent_y,'r');
   end

   if cent_y - rad > 1 && cent_x - rad > 1 && cent_y + rad < size(im, 1) && cent_x + rad < size(im,2)
       snip = im(cent_y-rad:cent_y+rad, cent_x-rad:cent_x+rad);
        
        % Set parameters to use as initial guess in the fitting. For the 
        % lattice data, try NeighborhoodSize = 1000, MaxThreshold = 2000, 
        % WidthGuess = 500, OffsetGuess = 1000.

        % For confocal data, try NeighborhoodSize = 1000, MaxThreshold = 20,
        % WidthGuess = 200, OffsetGuess = 10.
        
        NeighborhoodSize = 1000/pixelSize; %nm
        MaxThreshold = 30; %intensity
        WidthGuess = 200 / pixelSize; %nm
        OffsetGuess = 10; %intensity
        [fitting, rel_errors, GaussianIntensity] =  ...
            fitTwoGausses(snip, NeighborhoodSize, MaxThreshold, ...
            WidthGuess, OffsetGuess, show_status);

        disp(rel_errors);

        % Quality control.
        % TODO: make some quality control using the errors in
        % the fitting. Using any(rel_errors > 0.3) (for
        % example) is not the best thing to do, because
        % sometimes the second gaussian doesn't get a good fit
        % but the first one does, and the second one is good
        % enough to position its center.

        c_x = fitting(2) - rad + cent_x;
        c_y = fitting(4) - rad + cent_y;
%         int_x = [round(c_x - fitting(3)), round(c_x + fitting(3))];
%         int_y = [round(c_y - fitting(5)), round(c_y + fitting(5))];
        int_x = [round(c_x - 5), round(c_x + 5)];
        int_y = [round(c_y - 5), round(c_y + 5)];
        sigma_x = fitting(3);
        sigma_y = fitting(5);

        area = pi*sigma_x*sigma_y; %in pixels
        fixedAreaIntensity = 0;
        if int_x(1) > 1 && int_y(1) > 1 && int_x(2) < size(im,2) && int_y(2) < size(im,1)
            for w = int_x(1):int_x(2)
                for v = int_y(1): int_y(2)
                    fixedAreaIntensity = fixedAreaIntensity + double(im(v,w) - fitting(end));
                end
            end
            temp_particles = {{fixedAreaIntensity, c_x, c_y, fitting(end), snip, ...
                area, sigma_x, sigma_y, cent_y, cent_x, GaussianIntensity ,inten, max_dog}};
        else
            temp_particles = {[]};
        end
   else
       temp_particles = {[]};
   end
end
