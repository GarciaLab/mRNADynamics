function [fits, relative_errors, residual, confidence_intervals, GaussianIntensity, gaussian, mesh] = ...
    fitSingleGaussian(snippet, ~, ~, widthGuess, offsetGuess, show_status, graphicsHandles)

% Fit Gaussians to the given locus within a snippet

%     persistent gh;
%     gh = memoize(@gaussianForSpot);

warning('off','MATLAB:singularMatrix')

snippet = double(snippet);
[y,x] = meshgrid(1:size(snippet,2), 1:size(snippet,1));

med = median(snippet(:));
%fits: [amplitude, x position, x width, y position, y width, offset, angle]

singleGaussian = gaussianForSpot(y, x, snippet);

%Define some more initial parameters for fitting

% initial_parameters = [max(snippet(:)), round(length(snippet)/2), widthGuess, round(length(snippet)/2), ...
%     widthGuess,median(snippet(:)), 0, 0, 0];


initial_parameters = [max(snippet(:)), round(size(snippet, 2)/2), round(size(snippet, 1)/2), ...
    0, widthGuess, widthGuess,median(snippet(:)), 0, 0];


lsqOptions=optimset('Display','none');

%fits: [amplitude, x position, x width, y position, y width, offset, angle]
lb_offset = 1/10; %this is empirical. corresponds to a weak background of 1 pixel per ten having a value of 1.
% lb = [0, 0, 0, 0, 0,lb_offset, 0, -med/2, -med/2];
% ub = [max(snippet(:))*1.5, size(snippet, 2), size(snippet, 1), size(snippet, 2), size(snippet, 1), max(snippet(:)), 2*pi, med/2, med/2];
% 
% [single_fit, ~, residual, ~, ~, ~, ~] = lsqnonlin(singleGaussian, ...
%     initial_parameters,lb,ub, lsqOptions);

% @(A, x0, y0, rho, sigma_x, sigma_y, offset, offset_x, offset_y)
lb = [0, 0, 0, -1, 0, 0,lb_offset, 0, 0];
ub = [max(snippet(:))*1.5, size(snippet, 2), size(snippet, 1), 1, size(snippet, 2), size(snippet, 1), max(snippet(:)),med/2, med/2];


[single_fit, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(singleGaussian, ...
    initial_parameters,lb,ub, lsqOptions);


%quality control. probably just noisy background so amplitude
%should be close to offset in this case.
% 
% if .1 > single_fit(3) | .1 > single_fit(5)...
%         | single_fit(3) > size(snippet)/2 | single_fit(5) > size(snippet)/2
%     initial_parameters = [lb_offset, round(length(snippet)/2), widthGuess, round(length(snippet)/2), ...
%         widthGuess,med, 0, 0, 0];
%     lb = [0, 0, 0, 0, 0,lb_offset, 0, -med/2, -med/2];
%     ub = [median(snippet(:))/2, size(snippet, 1), size(snippet, 1), size(snippet, 2), size(snippet, 2), max(snippet(:)), 2*pi, med/2, med/2];
%     
%     [single_fit, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(singleGaussian, ...
%         initial_parameters,lb,ub, lsqOptions);
% end

confidence_intervals = single(nlparci(single_fit,residual,'jacobian',jacobian));
errors = zeros(1, length(single_fit));
for i = 1:length(confidence_intervals)
    errors(i) = abs((abs(confidence_intervals(i, 1)) - abs(confidence_intervals(i, 2)))/2);
end
relative_errors = single(abs(errors./single_fit));

residual = single(residual);

fits = single(single_fit);
GaussianIntensity = single(sum(sum(singleGaussian(single_fit) + snippet - single_fit(6))));

%Display
gaussian = single(singleGaussian(single_fit));
mesh = {uint8(y), uint8(x)};

if show_status && ~isempty(graphicsHandles)
    gAx = graphicsHandles(4);
    snipAx = graphicsHandles(6);
    rawAx = graphicsHandles(8);
    surf(gAx, y, x, gaussian + snippet);
    title(gAx, 'Single Gaussian fit')
    snipBig = imresize(snippet,10);
    imshow(snipBig,[], 'Parent', snipAx);
    surf(rawAx,y, x, snippet);
    title(rawAx,'Raw data');
end

end