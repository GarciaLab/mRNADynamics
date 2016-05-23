function temp_particles = fit_single_spot(k, im, im_label, neighb, rad, ...
    pixelSize, show_status)

[r,c] = find(im_label == k);

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
   temp_particles = {};
   if show_status
        ellipse(neighb/2,neighb/2,0,cent_y,cent_x,'r');
   end

   if cent_y - rad > 1 && cent_x - rad > 1 && cent_y + rad < size(im, 1) && cent_x + rad < size(im,2)
       snip = im(cent_y-rad:cent_y+rad, cent_x-rad:cent_x+rad);
%                      [f1, res1, f2, res2] = fitGausses(snip);
        
        % Set parameters to use as initial guess in the fitting. For the 
        % lattice data, try NeighborhoodSize = 1000, MaxThreshold = 2000, 
        % WidthGuess = 500, OffsetGuess = 1000.

        % For confocal data, try NeighborhoodSize = 500, MaxThreshold = 20,
        % WidthGuess = 100, OffsetGuess = 10.


        NeighborhoodSize = 1000/pixelSize; %nm
        MaxThreshold = 20; %intensity
        WidthGuess = 500 / pixelSize; %nm
        OffsetGuess = 1000; %intensity
        [f1, GaussianIntensity,res1, residual, exitflag, output, lambda, jacobian] =  ...
            fitTwoGausses(snip, NeighborhoodSize, MaxThreshold, ...
            WidthGuess, OffsetGuess, show_status);

        ci = nlparci(f1,residual,'jacobian',jacobian);
        errors = zeros(1, length(f1));
        for ndx = 1:length(ci)
            errors(ndx) = abs((abs(ci(ndx, 1)) - abs(ci(ndx, 2)))/2);
        end
        rel_errors = abs(errors./f1);
        disp(rel_errors);

        % Quality control.
        % TODO: make some quality control using the errors in
        % the fitting. Using any(rel_errors > 0.3) (for
        % example) is not the best thing to do, because
        % sometimes the second gaussian doesn't get a good fit
        % but the first one does, and the second one is good
        % enough to position its center.

%                     if f1(3) > rad+3 || f1(5) > rad+3
        if 1
            c_x = f1(2) - rad + cent_x;
            c_y = f1(4) - rad + cent_y;
%             int_x = [round(c_x - f1(3)), round(c_x + f1(3))];
%             int_y = [round(c_y - f1(5)), round(c_y + f1(5))];
            int_x = [round(c_x - 5), round(c_x + 5)];
            int_y = [round(c_y - 5), round(c_y + 5)];
            area = pi*f1(3)*f1(5); %in pixels
            fixedAreaIntensity = 0;
            if int_x(1) > 1 && int_y(1) > 1 && int_x(2) < size(im,2) && int_y(2) < size(im,1)
                for w = int_x(1):int_x(2)
                    for v = int_y(1): int_y(2)
                        fixedAreaIntensity = fixedAreaIntensity + double(im(v,w) - f1(end));
                    end
                end
                temp = {{fixedAreaIntensity, c_x, c_y, f1(6), snip, area, f1(3), f1(5),cent_y,cent_x,GaussianIntensity ,inten}};
                temp_particles = [temp_particles,temp];
            end
        end
    end
end
