function [fits, relative_errors, residual, confidence_intervals, GaussianIntensity, gaussian, mesh] = ...
    fitSingleGaussian(snippet, neighborhoodSize, threshold, widthGuess, offsetGuess, show_status)

% Fit Gaussians to the given locus within a snippet

warning('off','MATLAB:singularMatrix')

snippet = double(snippet);
[mesh_y,mesh_x] = meshgrid(1:size(snippet,2), 1:size(snippet,1));

%fits: [amplitude, x position, x width, y position, y width, offset, angle] 

singleGaussian = gaussianForSpot(snippet);

%Define some more initial parameters for fitting

    initial_parameters = [max(snippet(:)), round(length(snippet)/2), widthGuess, round(length(snippet)/2), ...
        widthGuess,offsetGuess, 0];

%Perform fitting
lsqOptions=optimset('Display','none',... %Inherited these options from Mikhail Tikhonov's FISH analysis
'maxfunevals',10000,...
'maxiter',10000); 

lb_offset = 1/10; %this is empirical. corresponds to a weak background of 1 pixel per ten having a value of 1. 
lb = [max(snippet(:))*.5, 0, 1, 0, 1,lb_offset, 0];
ub = [max(snippet(:))*2, size(snippet, 1), size(snippet, 1), size(snippet, 2), size(snippet, 2), max(snippet(:)), inf];

[single_fit, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(singleGaussian, ...
    initial_parameters,lb,ub, lsqOptions);

confidence_intervals = nlparci(single_fit,residual,'jacobian',jacobian);
errors = zeros(1, length(single_fit));
for i = 1:length(confidence_intervals)
    errors(i) = abs((abs(confidence_intervals(i, 1)) - abs(confidence_intervals(i, 2)))/2);
end
relative_errors = abs(errors./single_fit);
fits = single_fit; 
GaussianIntensity = sum(sum(singleGaussian(single_fit) + double(snippet) - single_fit(6)));

%Display
gaussian = singleGaussian(single_fit);
mesh = {mesh_y, mesh_x};

    if show_status
        figure(2);
        surf(mesh_y, mesh_x, gaussian + snippet);
        title('Single Gaussian fit')
        set(gcf,'units', 'normalized', 'position',[0.01, .55, .33, .33]);
        figure(3)
        snipBig = imresize(snippet,10);
        set(gcf,'units', 'normalized', 'position',[0.4, .2, .1, .1])
        imshow(snipBig,[]);
        figure(4);
        surf(mesh_y, mesh_x, snippet);
        title('Raw data');
        set(gcf,'units', 'normalized', 'position',[.01, .1, .33, .33]);
    end
end
