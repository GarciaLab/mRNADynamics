% Values calculated from polynomial fitting are used to create gaussians
% for each frame of a movie.
% Paul Marchando
% 12-3-2018

function [gauss_out, para_values] = fitDVGaussian(fake_dv_comb)

gauss_out = cell(fake_dv_comb(end,3),2);

% Finds the set of data for an individual frame and fits a gaussian
% with offset (+ d) to the points, storing the parameters and gof
% information in a cell array.

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

% Converts cell array of parameters to a double array for use in plotting.

for y = 1:length(gauss_out)
    
    if isempty(gauss_out{y,1}) == 0
        
        para_values(y,:) = coeffvalues(gauss_out{y,1});
        
    end
    
end

para_values = [transpose(1:length(para_values)) para_values];
para_values = para_values(para_values(:,2) ~= 0,:);

end
