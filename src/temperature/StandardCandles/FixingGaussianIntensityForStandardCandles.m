Prefix = '2022-01-19-eVasax5-120mer-GFP-lamin-T20C_856x856_LA1_Zoom2_zstep025_50uW-E1';
outdir = 'S:/Gabriella/Dropbox\StandardCandles\GaussianIntensityTesting';
z4 = load([outdir '\Z4Test120220126_1421_j2064fitSingleGaussian.mat']);
z5 = load([outdir '\Z5Test120220126_1421_j2064fitSingleGaussian.mat']);
z6 = load([outdir '\Z6Test120220126_1421_j2064fitSingleGaussian.mat']);
z7 = load([outdir '\Z7Test120220126_1421_j2064fitSingleGaussian.mat']);
% function [fits, relative_errors, residual, confidence_intervals, GaussianIntensity, gaussian, mesh] = ...
%     fitSingleGaussian(snippet, ~, ~, widthGuess, offsetGuess, show_status, graphicsHandles)

% Fit Gaussians to the given locus within a snippet
%%
warning('off','MATLAB:singularMatrix')

z4snippet = double(z4.snippet(11:17,11:17));

yDim = size(z4snippet, 1);
xDim = size(z4snippet, 2);
m4 = max(z4snippet(:));

%let's cache this for efficiency
%realized caching does weird stuff in parpools
%gonna disable that til i figure out why. 
% persistent y x 
%if isempty(y), 
    [y,x] = meshgrid(1:size(z4snippet,2), 1:size(z4snippet,1));
%end

% [y,x] = meshgrid(1:xDim, 1:yDim);

med4 = median(z4snippet(:));
%fits: [amplitude, x position, x width, y position, y width, offset, angle]

singleGaussian4 = gaussianForSpot(y, x, z4snippet);
widthGuess = 1;

%Define some more initial parameters for fitting
initial_parameters = [m4, round(xDim/2), round(yDim/2), ...
    0, widthGuess, widthGuess,med4, 0, 0];

%fits: [amplitude, x position, x width, y position, y width, offset, angle]
lb_offset = 1/10; %this is empirical. corresponds to a weak background of 1 pixel per ten having a value of 1.
% lb = [0, 0, 0, 0, 0,lb_offset, 0, -med/2, -med/2];
% ub = [max(snippet(:))*1.5, size(snippet, 2), size(snippet, 1), size(snippet, 2), size(snippet, 1), max(snippet(:)), 2*pi, med/2, med/2];
%
% [single_fit, ~, residual, ~, ~, ~, ~] = lsqnonlin(singleGaussian, ...
%     initial_parameters,lb,ub, lsqOptions);

% @(A, x0, y0, rho, sigma_x, sigma_y, offset, offset_x, offset_y)
%lb = [0, 0, 0, -1, 0, 0,lb_offset, 0, 0];
% CHANGE
lb = [0, 0, 0, -1, 0, 0,0, 0, 0];
ub = [m4*1.5, xDim, yDim, 1,...
    3, 3, m4, m4/2, m4/2];

%let's cache this for efficiency
% persistent lsqOptions 
% CHANGE
%if isempty(lsqOptions), 
    lsqOptions=optimset('Display','none'); 
    %end

[single_fit4, ~, residual4, ~, ~, ~, jacobian4] = lsqnonlin(singleGaussian4, ...
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

confidence_intervals4 = single(nlparci(single_fit4,residual4,'jacobian',jacobian4));
errors4 = zeros(1, length(single_fit4));
for i = 1:length(confidence_intervals4)
    errors4(i) = (1/2)* abs(...
        (abs( confidence_intervals4(i, 1) )...
        - abs( confidence_intervals4(i, 2)) )...
        );
end

relative_errors4 = single(abs(errors4./single_fit4));

residual4 = single(residual4);

fits4 = single(single_fit4);
% Shouldn't this also subtract off any background associated with
% single_fit(8) and single_fit(9)? (GM: 1/16/22)
GaussianIntensityZ4 = single(sum(sum(singleGaussian4(single_fit4) + z4snippet - single_fit4(7)-...
    single_fit4(8)*x-single_fit4(9)*y)));

%Display
gaussian4 = single(singleGaussian4(single_fit4));
mesh = {uint16(y), uint16(x)};

% if show_status && ~isempty(graphicsHandles)
%     gAx = graphicsHandles(4);
%     snipAx = graphicsHandles(6);
%     rawAx = graphicsHandles(8);
%     surf(gAx, y, x, gaussian + snippet);
%     title(gAx, 'Single Gaussian fit')
%     snipBig = imresize(snippet,10);
%     imshow(snipBig,[], 'Parent', snipAx);
%     surf(rawAx,y, x, snippet);
%     title(rawAx,'Raw data');
% end

%%
z5snippet = double(z5.snippet(11:17,11:17));

yDim = size(z5snippet, 1);
xDim = size(z5snippet, 2);
m5 = max(z5snippet(:));

%let's cache this for efficiency
%realized caching does weird stuff in parpools
%gonna disable that til i figure out why. 
% persistent y x 
%if isempty(y), 
    [y,x] = meshgrid(1:size(z5snippet,2), 1:size(z5snippet,1));
%end

% [y,x] = meshgrid(1:xDim, 1:yDim);

med5 = median(z5snippet(:));
%fits: [amplitude, x position, x width, y position, y width, offset, angle]

singleGaussian5 = gaussianForSpot(y, x, z5snippet);
widthGuess = 1;

%Define some more initial parameters for fitting
initial_parameters = [m5, round(xDim/2), round(yDim/2), ...
    0, widthGuess, widthGuess,med5, 0, 0];

%fits: [amplitude, x position, x width, y position, y width, offset, angle]
lb_offset = 1/10; %this is empirical. corresponds to a weak background of 1 pixel per ten having a value of 1.
% lb = [0, 0, 0, 0, 0,lb_offset, 0, -med/2, -med/2];
% ub = [max(snippet(:))*1.5, size(snippet, 2), size(snippet, 1), size(snippet, 2), size(snippet, 1), max(snippet(:)), 2*pi, med/2, med/2];
%
% [single_fit, ~, residual, ~, ~, ~, ~] = lsqnonlin(singleGaussian, ...
%     initial_parameters,lb,ub, lsqOptions);

% @(A, x0, y0, rho, sigma_x, sigma_y, offset, offset_x, offset_y)
%lb = [0, 0, 0, -1, 0, 0,lb_offset, 0, 0];
% CHANGE
lb = [0, 0, 0, -1, 0, 0,0, 0, 0];
ub = [m5*1.5, xDim, yDim, 1,...
    xDim, yDim, m5, m5/2, m5/2];

%let's cache this for efficiency
% persistent lsqOptions 
% CHANGE
%if isempty(lsqOptions), 
    lsqOptions=optimset('Display','none'); 
    %end

[single_fit5, ~, residual5, ~, ~, ~, jacobian5] = lsqnonlin(singleGaussian5, ...
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

confidence_intervals5 = single(nlparci(single_fit5,residual5,'jacobian',jacobian5));
errors5 = zeros(1, length(single_fit5));
for i = 1:length(confidence_intervals5)
    errors5(i) = (1/2)* abs(...
        (abs( confidence_intervals5(i, 1) )...
        - abs( confidence_intervals5(i, 2)) )...
        );
end

relative_errors5 = single(abs(errors5./single_fit5));

residual5 = single(residual5);

fits5 = single(single_fit5);
% Shouldn't this also subtract off any background associated with
% single_fit(8) and single_fit(9)? (GM: 1/16/22)
GaussianIntensityZ5 = single(sum(sum(singleGaussian5(single_fit5) + z5snippet - single_fit5(7)-...
    single_fit5(8)*x-single_fit5(9)*y)));

%Display
gaussian5 = single(singleGaussian5(single_fit5));
mesh = {uint16(y), uint16(x)};


%%
z6snippet = double(z6.snippet(11:17,11:17));

yDim = size(z6snippet, 1);
xDim = size(z6snippet, 2);
m6 = max(z6snippet(:));

%let's cache this for efficiency
%realized caching does weird stuff in parpools
%gonna disable that til i figure out why. 
% persistent y x 
%if isempty(y), 
    [y,x] = meshgrid(1:size(z6snippet,2), 1:size(z6snippet,1));
%end

% [y,x] = meshgrid(1:xDim, 1:yDim);

med6 = median(z6snippet(:));
%fits: [amplitude, x position, x width, y position, y width, offset, angle]

singleGaussian6 = gaussianForSpot(y, x, z6snippet);
widthGuess = 1;

%Define some more initial parameters for fitting
initial_parameters = [m6, round(xDim/2), round(yDim/2), ...
    0, widthGuess, widthGuess,med6, 0, 0];

%fits: [amplitude, x position, x width, y position, y width, offset, angle]
lb_offset = 1/10; %this is empirical. corresponds to a weak background of 1 pixel per ten having a value of 1.
% lb = [0, 0, 0, 0, 0,lb_offset, 0, -med/2, -med/2];
% ub = [max(snippet(:))*1.5, size(snippet, 2), size(snippet, 1), size(snippet, 2), size(snippet, 1), max(snippet(:)), 2*pi, med/2, med/2];
%
% [single_fit, ~, residual, ~, ~, ~, ~] = lsqnonlin(singleGaussian, ...
%     initial_parameters,lb,ub, lsqOptions);

% @(A, x0, y0, rho, sigma_x, sigma_y, offset, offset_x, offset_y)
%lb = [0, 0, 0, -1, 0, 0,lb_offset, 0, 0];
% CHANGE
lb = [0, 0, 0, -1, 0, 0,0, 0, 0];
ub = [m6*1.5, xDim, yDim, 1,...
    xDim, yDim, m6, m6/2, m6/2];

%let's cache this for efficiency
% persistent lsqOptions 
% CHANGE
%if isempty(lsqOptions), 
    lsqOptions=optimset('Display','none'); 
    %end

[single_fit6, ~, residual6, ~, ~, ~, jacobian6] = lsqnonlin(singleGaussian6, ...
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

confidence_intervals6 = single(nlparci(single_fit6,residual6,'jacobian',jacobian6));
errors6 = zeros(1, length(single_fit6));
for i = 1:length(confidence_intervals6)
    errors6(i) = (1/2)* abs(...
        (abs( confidence_intervals6(i, 1) )...
        - abs( confidence_intervals6(i, 2)) )...
        );
end

relative_errors6 = single(abs(errors6./single_fit6));

residual6 = single(residual6);

fits6 = single(single_fit6);
% Shouldn't this also subtract off any background associated with
% single_fit(8) and single_fit(9)? (GM: 1/16/22)
GaussianIntensityZ6 = single(sum(sum(singleGaussian6(single_fit6) + z6snippet - single_fit6(7)-...
    single_fit6(8)*x-single_fit6(9)*y)));

%Display
gaussian6 = single(singleGaussian6(single_fit6));
mesh = {uint16(y), uint16(x)};

%%
z7snippet = double(z7.snippet(11:17,11:17));

yDim = size(z7snippet, 1);
xDim = size(z7snippet, 2);
m7 = max(z7snippet(:));

%let's cache this for efficiency
%realized caching does weird stuff in parpools
%gonna disable that til i figure out why. 
% persistent y x 
%if isempty(y), 
    [y,x] = meshgrid(1:size(z7snippet,2), 1:size(z7snippet,1));
%end

% [y,x] = meshgrid(1:xDim, 1:yDim);

med7 = median(z7snippet(:));
%fits: [amplitude, x position, x width, y position, y width, offset, angle]

singleGaussian7 = gaussianForSpot(y, x, z7snippet);
widthGuess = 1;

%Define some more initial parameters for fitting
initial_parameters = [m7, round(xDim/2), round(yDim/2), ...
    0, widthGuess, widthGuess,med7, 0, 0];

%fits: [amplitude, x position, x width, y position, y width, offset, angle]
lb_offset = 1/10; %this is empirical. corresponds to a weak background of 1 pixel per ten having a value of 1.
% lb = [0, 0, 0, 0, 0,lb_offset, 0, -med/2, -med/2];
% ub = [max(snippet(:))*1.5, size(snippet, 2), size(snippet, 1), size(snippet, 2), size(snippet, 1), max(snippet(:)), 2*pi, med/2, med/2];
%
% [single_fit, ~, residual, ~, ~, ~, ~] = lsqnonlin(singleGaussian, ...
%     initial_parameters,lb,ub, lsqOptions);

% @(A, x0, y0, rho, sigma_x, sigma_y, offset, offset_x, offset_y)
%lb = [0, 0, 0, -1, 0, 0,lb_offset, 0, 0];
% CHANGE
lb = [0, 0, 0, -1, 0, 0,0, 0, 0];
ub = [m7*1.5, xDim, yDim, 1,...
    xDim, yDim, m7, m7/2, m7/2];

%let's cache this for efficiency
% persistent lsqOptions 
% CHANGE
%if isempty(lsqOptions), 
    lsqOptions=optimset('Display','none'); 
    %end

[single_fit7, ~, residual7, ~, ~, ~, jacobian7] = lsqnonlin(singleGaussian7, ...
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

confidence_intervals7 = single(nlparci(single_fit7,residual7,'jacobian',jacobian7));
errors7 = zeros(1, length(single_fit7));
for i = 1:length(confidence_intervals7)
    errors7(i) = (1/2)* abs(...
        (abs( confidence_intervals7(i, 1) )...
        - abs( confidence_intervals7(i, 2)) )...
        );
end

relative_errors7 = single(abs(errors7./single_fit7));

residual7 = single(residual7);

fits7 = single(single_fit7);
% Shouldn't this also subtract off any background associated with
% single_fit(8) and single_fit(9)? (GM: 1/16/22)
GaussianIntensityZ7 = single(sum(sum(singleGaussian7(single_fit7) + z7snippet - single_fit7(7)-...
    single_fit7(8)*x-single_fit7(9)*y)));

%Display
gaussian7 = single(singleGaussian6(single_fit7));
mesh = {uint16(y), uint16(x)};




