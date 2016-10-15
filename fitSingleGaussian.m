function [fits, relative_errors, confidence_intervals, GaussianIntensity, gaussian, mesh] = ...
    fitSingleGaussian(snippet, neighborhoodSize, threshold, widthGuess, offsetGuess, show_status)

% Fit Gaussians to the given locus within a snippet

snippet = double(snippet);
[mesh_y,mesh_x] = meshgrid(1:size(snippet,2), 1:size(snippet,1));

%Single and double Gaussians with the full range of parameter space,
%including an angular parameter

singleGaussian = @(params) params(1).*exp( -( ...
    -( (cos(params(7))^2 / 2*params(3)^2) + (sin(params(7))^2 / 2*params(5)^2) )*(x-params(2)).^2 ...
    -2*( (sin(2*params(7)) / 4*params(3)^2) + (sin(2*params(7)) / 4*params(5)^2)) * (x-params(2)).*(y-params(4))...
    + ( (sin(params(7))^2 / 2*params(3)^2) + (cos(params(7))^2 / 2*params(5)^2) ) *(y-params(4)).^2))...
    + params(6) - double(snip);
    
neighborhoodSize = 2*floor(neighborhoodSize/2) + 1;
    
hLocalMax = vision.LocalMaximaFinder;
hLocalMax.NeighborhoodSize = [neighborhoodSize, neighborhoodSize];
hLocalMax.Threshold = threshold;
centers = double(step(hLocalMax, snippet));    

initial_parameters = [max(max(snippet)), centers(1,2), widthGuess, centers(1,1), ...
        widthGuess,offsetGuess, 0];
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
relative_errors = abs(errors./double_fit);

snip_cent = size(snippet)/2;
fits = single_fit; 
GaussianIntensity = sum(sum(singleGaussian(single_fit) + double(snippet) - single_fit(end)));

% 
% if show_status
%     figure(2)
%     if size(centers,1) == 1
%         surf(y, x, singleGaussian(fits) + double(snip));
%         title('Single gaussian fits')
%     elseif size(centers,1) == 2
%         surf(y, x, doubleGaussian(double_fit) + double(snip));
%         title('Double gaussian fits')
%     end
%     figure(3)
%     imshow(imresize(snip,10),[]);
%     figure(4)
%     surf(y, x, double(snip));
%     title('Raw data');

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
