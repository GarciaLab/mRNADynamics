function [fits, relative_errors, confidence_intervals, GaussianIntensity, gaussian, mesh] = ...
    fitGaussians(snippet, neighborhoodSize, threshold, widthGuess, offsetGuess, show_status)

% Fit Gaussians to the given locus within a snippet

% snippet = CPsmooth(snippet,'Gaussian Filter',1.3,0); %AR 7/6/16: Not sure this
% is needed.

snippet = double(snippet);
[mesh_y,mesh_x] = meshgrid(1:size(snippet,2), 1:size(snippet,1));

singleGaussian = @(params) params(1).*exp((-1/2).*(((mesh_x-params(2))./params(3)).^2 ...
        + ((mesh_y-params(4))./params(5)).^2)) + params(6) - double(snippet);
    
%Single and double Gaussians with the full range of parameter space,
%including an angular parameter

% singleGaussian = @(params) params(1).*exp( -( ...
%     -( (cos(params(7))^2 / 2*params(3)^2) + (sin(params(7))^2 / 2*params(5)^2) )*(x-params(2)).^2 ...
%     -2*( (sin(2*params(7)) / 4*params(3)^2) + (sin(2*params(7)) / 4*params(5)^2)) * (x-params(2)).*(y-params(4))...
%     + ( (sin(params(7))^2 / 2*params(3)^2) + (cos(params(7))^2 / 2*params(5)^2) ) *(y-params(4)).^2))...
%     + params(6) - double(snip);
% doubleGaussian = @(params) params(1).*exp( -( ...
%     -( (cos(params(12))^2 / 2*params(3)^2) + (sin(params(12))^2 / 2*params(5)^2) )*(x-params(2)).^2 ...
% -2*( (sin(2*params(12)) / 4*params(3)^2) + (sin(2*params(12)) / 4*params(5)^2)) * (x-params(2)).*(y-params(4))...
%     + ( (sin(params(12))^2 / 2*params(3)^2) + (cos(params(12))^2 / 2*params(5)^2) ) *(y-params(4)).^2))...
%     + params(6).*exp( -( ...
%     -( (cos(params(13))^2 / 2*params(8)^2) + (sin(params(13))^2 / 2*params(10)^2) )*(x-params(7)).^2 ...
% -2*( (sin(2*params(13)) / 4*params(8)^2) + (sin(2*params(13)) / 4*params(10)^2)) * (x-params(7)).*(y-params(9))...
%     + ( (sin(params(13))^2 / 2*params(8)^2) + (cos(params(13))^2 / 2*params(10)^2) ) *(y-params(9)).^2))...
%     + params(11) - double(snip);

%Double elliptical Gaussian with the elliptical axes parallel to the x-axis

% doubleGaussian = @(params) params(1).*exp((-1/2).*(((x-params(2))./params(3)).^2 ... 
%         + ((y-params(4))./params(5)).^2)) ... 
%         + params(6).*exp((-1/2).*(((x-params(7))./params(8)).^2  ...
%         + ((y-params(9))./params(10)).^2))+ params(11) - double(snip);
    
%For now I want to assume circular Gaussians. Only 9 free parameters! so
% easy
% doubleGaussian = @(params) params(1).*exp((-1/2).*(((x-params(2))./params(3)).^2 ... 
%         + ((y-params(4))./params(3)).^2)) ... 
%         + params(6).*exp((-1/2).*(((x-params(7))./params(8)).^2  ...
%         + ((y-params(9))./params(8)).^2))+ params(11) - double(snip);

doubleGaussian = @(params) params(1).*exp((-1/2).*(((mesh_x-params(2))./params(3)).^2 ... 
        + ((mesh_y-params(4))./params(3)).^2)) ... 
        + params(5).*exp((-1/2).*(((mesh_x-params(6))./params(7)).^2  ...
        + ((mesh_y-params(8))./params(7)).^2))+ params(9) - double(snippet);
    
neighborhoodSize = 2*floor(neighborhoodSize/2) + 1;
    
hLocalMax = vision.LocalMaximaFinder;
hLocalMax.NeighborhoodSize = [neighborhoodSize, neighborhoodSize];
hLocalMax.Threshold = threshold;
centers = double(step(hLocalMax, snippet));    

if size(centers,1)== 2
    initial_parameters = [max(max(snippet)), centers(1,2), widthGuess, centers(1,1), ...
            max(max(snippet)), centers(2,2), widthGuess, centers(2,1), ...
            offsetGuess];
elseif size(centers, 1) == 1
    initial_parameters = [max(max(snippet)), centers(1,2), widthGuess, centers(1,1), ...
            max(max(snippet)), centers(1,2), widthGuess, centers(1,1), ...
            offsetGuess];
else
    initial_parameters = [max(max(snippet)), round(size(snippet,1)/2), widthGuess, round(size(snippet,1)/2), ...
            max(max(snippet)), round(size(snippet,1)/2), widthGuess, round(size(snippet,1)/2), ...
            offsetGuess];
end

    lsqOptions=optimset('Display','none',... %Inherited these options from Mikhail Tikhonov's FISH analysis
    'maxfunevals',10000,...
    'maxiter',10000);

    [double_fit, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(doubleGaussian, ...
        initial_parameters,zeros(1,9),inf(1,9), lsqOptions);
    
    confidence_intervals = nlparci(double_fit,residual,'jacobian',jacobian);
    errors = zeros(1, length(double_fit));
    for i = 1:length(confidence_intervals)
        errors(i) = abs((abs(confidence_intervals(i, 1)) - abs(confidence_intervals(i, 2)))/2);
    end
    relative_errors = abs(errors./double_fit);
        
    snip_cent = size(snippet)/2;
    gaussian1_center = [double_fit(2), double_fit(4)];
    gaussian2_center = [double_fit(6), double_fit(8)];
    dif1 = gaussian1_center - snip_cent;
    dif2 = gaussian2_center - snip_cent;
    distance1 = sqrt(sum(abs(dif1).^2,2));
    distance2 = sqrt(sum(abs(dif2).^2,2));
    fits = double_fit; 
    %find the distance between sister chromatids 
    sister_chromatid_distance = sqrt((fits(2)-fits(6))^2 + (fits(4) - fits(8))^2); % in pixels
    fits(end+1) = sister_chromatid_distance;

%AR 7/6/2016: Why is this necessary? We already found local maxima in the
%identifySpot script. And it appears the algorithm is the same.

% Round to the nearest odd integer to be able to feed it to
% LocalMaximaFinder

% NeighborhoodSize = 2*floor(NeighborhoodSize/2) + 1;
%     
% hLocalMax = vision.LocalMaximaFinder;
% hLocalMax.NeighborhoodSize = [NeighborhoodSize, NeighborhoodSize];
% hLocalMax.Threshold = Threshold;
% centers = double(step(hLocalMax, snip));
% 
% if size(centers,1) == 1
%     init_params = [max(max(snip)), centers(1,2), WidthGuess, ... 
%         centers(1,1), WidthGuess, OffsetGuess, 0];
% 
%     [fits, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(singleGaussian, ...
%         init_params,zeros(1,6),inf(1,6));
%     
%     ci = nlparci(fits,residual,'jacobian',jacobian);
%     errors = zeros(1, length(fits));
%     for ndx = 1:length(ci)
%         errors(ndx) = abs((abs(ci(ndx, 1)) - abs(ci(ndx, 2)))/2);
%     end
%     rel_errors = abs(errors./fits);
%     
% elseif size(centers,1) == 2
%     init_params = [max(max(snip)), centers(1,2), WidthGuess, centers(1,1), ...
%         WidthGuess, max(max(snip)), centers(2,2), WidthGuess, centers(2,1), ...
%         WidthGuess, OffsetGuess];
% 
%     [double_fit, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(doubleGaussian, ...
%         init_params,zeros(1,11),inf(1,11));
%     
%     ci = nlparci(double_fit,residual,'jacobian',jacobian);
%     errors = zeros(1, length(double_fit));
%     for ndx = 1:length(ci)
%         errors(ndx) = abs((abs(ci(ndx, 1)) - abs(ci(ndx, 2)))/2);
%     end
%     rel_errors = abs(errors./double_fit);
%     
%     % Keep only the gaussian closer to the center of the snippet.
%     
%     snip_cent = size(snip)./2;
%     gaussian1_cent = [double_fit(2), double_fit(4)];
%     gaussian2_cent = [double_fit(7), double_fit(9)];
%     dif1 = gaussian1_cent - snip_cent;
%     dif2 = gaussian2_cent - snip_cent;
%     distance1 = sqrt(sum(abs(dif1).^2,2));
%     distance2 = sqrt(sum(abs(dif2).^2,2));
%     
%     if distance1 < distance2
%         fits = double_fit(1:5);
%         fits(6) = double_fit(end);
%     else
%         fits = double_fit(6:end);
%     end
% else
%     % TODO: I'm just making this so the script doesn't crash when it can't
%     % find any maxima. But it's ugly and should be replaced with something
%     % like assigning NaNs to all the return values.
%     init_params = [2000, 10, 5, 10, 5,1000, 0];
%     [fits, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(singleGaussian, ...
%         init_params,zeros(1,6),inf(1,6));
%     
%     ci = nlparci(fits,residual,'jacobian',jacobian);
%     errors = zeros(1, length(fits));
%     for ndx = 1:length(ci)
%         errors(ndx) = abs((abs(ci(ndx, 1)) - abs(ci(ndx, 2)))/2);
%     end
%     rel_errors = abs(errors./fits);
% end

%Compute intensities by integrating over the Gaussian fit. Offset
%subtracted

GaussianIntensity = sum(sum(doubleGaussian(double_fit) + double(snippet) - double_fit(end)));

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

gaussian = doubleGaussian(double_fit);
mesh = {mesh_y, mesh_x};

    if show_status
        figure(2);
        surf(mesh_y, mesh_x, gaussian + snippet);
        title('Double Gaussian fits')
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
