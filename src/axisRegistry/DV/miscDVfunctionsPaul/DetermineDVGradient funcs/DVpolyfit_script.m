% Fits 2nd order polynomials to fluorescence data and uses fit to calculate
% fluorescence. These values are then used to create gaussians for each
% frame of a movie.
% Paul Marchando
% 11-28-2018

polyfit_data = cell2mat(working_dat_final(:,1));
polyfit_out = zeros(length(polyfit_data),3);

thresh = input('Please input a threshold value for ''a'': ');

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

gauss_out = cell(fake_dv_comb(end,3),2);

tic
for m = 1:fake_dv_comb(end,3)
    
    temp_range = find(fake_dv_comb(:,3) == m);
    
    if isempty(temp_range) == 0 && length(temp_range) >= 4
        
        mygauss = fittype('a1*exp(-((x-b1)/c1)^2) + d1',...
            'dependent',{'y'},'independent',{'x'},...
            'coefficients',{'a1','b1','c1','d1'});
        options = fitoptions(mygauss);
        options.Lower = [0 -Inf 100 0];
        options.Upper = [Inf Inf 500 Inf];
        [temp_gauss, temp_gof] = fit(fake_dv_comb(temp_range(1):temp_range(end),1),...
            fake_dv_comb(temp_range(1):temp_range(end),2),mygauss,options);
        gauss_out{m,1} = temp_gauss;
        gauss_out{m,2} = temp_gof;
        
    end
    
end
toc

para_values = zeros(length(gauss_out),4);

for y = 1:length(gauss_out)
    
    para_values(y,:) = coeffvalues(gauss_out{y,1});
    
end
