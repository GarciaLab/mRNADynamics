% Fits 2nd order polynomials to fluorescence data and uses fit to calculate
% fluorescence. Must enter a threshold for the 'a' value of the fit, below
% which the mean of the middle 11 frames will be used to calculate maximum
% fluorescence instead of the polynomial maximum.
% Paul Marchando
% 11-28-2018

function [fake_dv_comb, fake_dv_check, polyfit_data, polyfit_out, ...
    polyfit_maximums] = DVpolyfit(working_dat_final,thresh)

% Fits a polynomial to nuclear fluorescence data. If the "a" value of the
% polynomial is less than a specified threshold or the maximum of the
% fitted polynomial isn't within the frames provided, the mean of the
% middle ten nuclear frames is computed as the fluorescence value.

polyfit_data = cell2mat(working_dat_final(:,1));
polyfit_out = zeros(length(polyfit_data),3);

tic
for n = 1:length(polyfit_data)
    
    polyfit_data(n,:) = transpose(smooth(polyfit_data(n,:)));
    polyfit_out(n,:) = polyfit(1:length(polyfit_data(n,:)),polyfit_data(n,:),2);
    
end
toc

polyfit_maximums_x = -(polyfit_out(:,2)) ./ (2 * polyfit_out(:,1));

polyfit_maximums_y = (4 * (polyfit_out(:,1)) .* (polyfit_out(:,3))...
    - (polyfit_out(:,2)).^2) ./ (4 * polyfit_out(:,1));

polyfit_maximums = [polyfit_maximums_x polyfit_maximums_y];

fake_dv_comb = [cell2mat(working_dat_final(:,2)) polyfit_maximums(:,2) ...
    cell2mat(working_dat_final(:,3))];

filter_max = polyfit_maximums(:,1) < 0 | polyfit_maximums(:,1) > 21;
filter_thresh = abs(polyfit_out(:,1)) < thresh;
filter_tot = filter_max | filter_thresh;
tot_indices = 1:length(polyfit_out);
tot_indices = tot_indices(filter_tot);

for j = 1:length(tot_indices)
    
    fake_dv_comb(tot_indices(j),2) = mean(polyfit_data(tot_indices(j),6:16));
    
end

fake_dv_check = fake_dv_comb;
fake_dv_comb = sortrows(fake_dv_comb,3);

end
