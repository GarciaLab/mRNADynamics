function [fits, relative_errors, residual, confidence_intervals, GaussianIntensity, gaussian, mesh] = ...
    fitSingleGaussian(snippet, neighborhoodSize, threshold, widthGuess, offsetGuess, show_status)

% Fit Gaussians to the given locus within a snippet

warning('off','MATLAB:singularMatrix')

snippet = double(snippet);
[mesh_y,mesh_x] = meshgrid(1:size(snippet,2), 1:size(snippet,1));

% Single gaussian function
singleGaussian = @(params) (params(1).*...
    exp(-(...
    (((cos(params(7)))^2 / (2*params(3)^2) ) + ((sin(params(7)))^2 / 2*params(5)^2))  .* (mesh_x-params(2)).^2 ...
    - 2*((-sin(2*params(7)) / (4*params(3)^2) ) + (sin(2*params(7)) / 4*params(5)^2)) .* (mesh_x-params(2)).*(mesh_y-params(4))...
    + (((sin(params(7)))^2 / (2*params(3)^2) ) + ((cos(params(7)))^2 / 2*params(5)^2)).* (mesh_y-params(4)).^2 ...
        )))...
    + params(6) - double(snippet);

%Define some more initial parameters for fitting
neighborhoodSize = 2*floor(neighborhoodSize/2) + 1; %notice that this is forced to be odd
    
hLocalMax = vision.LocalMaximaFinder;
hLocalMax.NeighborhoodSize = [neighborhoodSize, neighborhoodSize];
hLocalMax.Threshold = threshold;
centers = double(step(hLocalMax, snippet));    
if ~isempty(centers)
    initial_parameters = [max(max(snippet)), centers(1,2), widthGuess, centers(1,1), ...
        widthGuess,offsetGuess, 0];
else
    initial_parameters = [max(max(snippet)), round(length(snippet)/2), widthGuess, round(length(snippet)/2), ...
        widthGuess,offsetGuess, 0];
end

%Perform fitting
lsqOptions=optimset('Display','none',... %Inherited these options from Mikhail Tikhonov's FISH analysis
'maxfunevals',10000,...
'maxiter',10000); 

[single_fit, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(singleGaussian, ...
    initial_parameters,zeros(1,7),inf(1,7), lsqOptions);

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
