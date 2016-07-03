function [fits, rel_errors, GaussianIntensity] = ...
    fitTwoGausses(snip, NeighborhoodSize, Threshold, WidthGuess, OffsetGuess, show_status)

% Find local maxima in snip and use that information to decide if fits
% one or two gaussians. Also, use that information to define a reasonable 
% initial guess for the fits initial parameters.

snip = CPsmooth(snip,'Gaussian Filter',1.3,0);
snip = double(snip);
[y,x] = meshgrid(1:size(snip,2), 1:size(snip,1));

singleGaussian = @(params) params(1).*exp((-1/2).*(((x-params(2))./params(3)).^2 ...
        + ((y-params(4))./params(5)).^2)) + params(6) - double(snip);

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
    
doubleGaussian = @(params) params(1).*exp((-1/2).*(((x-params(2))./params(3)).^2 ... 
        + ((y-params(4))./params(5)).^2)) ... 
        + params(6).*exp((-1/2).*(((x-params(7))./params(8)).^2  ...
        + ((y-params(9))./params(10)).^2))+ params(11) - double(snip);

% Round to the nearest odd integer to be able to feed it to
% LocalMaximaFinder

NeighborhoodSize = 2*floor(NeighborhoodSize/2) + 1;
    
hLocalMax = vision.LocalMaximaFinder;
hLocalMax.NeighborhoodSize = [NeighborhoodSize, NeighborhoodSize];
hLocalMax.Threshold = Threshold;
centers = double(step(hLocalMax, snip));

if size(centers,1) == 1
    init_params = [max(max(snip)), centers(1,2), WidthGuess, ... 
        centers(1,1), WidthGuess, OffsetGuess, 0];

    [fits, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(singleGaussian, ...
        init_params,zeros(1,6),inf(1,6));
    
    ci = nlparci(fits,residual,'jacobian',jacobian);
    errors = zeros(1, length(fits));
    for ndx = 1:length(ci)
        errors(ndx) = abs((abs(ci(ndx, 1)) - abs(ci(ndx, 2)))/2);
    end
    rel_errors = abs(errors./fits);
    
elseif size(centers,1) == 2
    init_params = [max(max(snip)), centers(1,2), WidthGuess, centers(1,1), ...
        WidthGuess, max(max(snip)), centers(2,2), WidthGuess, centers(2,1), ...
        WidthGuess, OffsetGuess];

    [double_fit, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(doubleGaussian, ...
        init_params,zeros(1,11),inf(1,11));
    
    ci = nlparci(double_fit,residual,'jacobian',jacobian);
    errors = zeros(1, length(double_fit));
    for ndx = 1:length(ci)
        errors(ndx) = abs((abs(ci(ndx, 1)) - abs(ci(ndx, 2)))/2);
    end
    rel_errors = abs(errors./double_fit);
    
    % Keep only the gaussian closer to the center of the snippet.
    
    snip_cent = size(snip)./2;
    gaussian1_cent = [double_fit(2), double_fit(4)];
    gaussian2_cent = [double_fit(7), double_fit(9)];
    dif1 = gaussian1_cent - snip_cent;
    dif2 = gaussian2_cent - snip_cent;
    distance1 = sqrt(sum(abs(dif1).^2,2));
    distance2 = sqrt(sum(abs(dif2).^2,2));
    
    if distance1 < distance2
        fits = double_fit(1:5);
        fits(6) = double_fit(end);
    else
        fits = double_fit(6:end);
    end
else
    % TODO: I'm just making this so the script doesn't crash when it can't
    % find any maxima. But it's ugly and should be replaced with something
    % like assigning NaNs to all the return values.
    init_params = [2000, 10, 5, 10, 5,1000, 0];
    [fits, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(singleGaussian, ...
        init_params,zeros(1,6),inf(1,6));
    
    ci = nlparci(fits,residual,'jacobian',jacobian);
    errors = zeros(1, length(fits));
    for ndx = 1:length(ci)
        errors(ndx) = abs((abs(ci(ndx, 1)) - abs(ci(ndx, 2)))/2);
    end
    rel_errors = abs(errors./fits);
end

%Compute intensities by integrating over the Gaussian fit. Offset
%subtracted
% TODO: consider the case with two gaussians here

GaussianIntensity = sum(sum(singleGaussian(fits) + double(snip) - fits(6)));

if show_status
    figure(2)
    if size(centers,1) == 1
        surf(y, x, singleGaussian(fits) + double(snip));
        title('Single gaussian fits')
    elseif size(centers,1) == 2
        surf(y, x, doubleGaussian(double_fit) + double(snip));
        title('Double gaussian fits')
    end
    figure(3)
    imshow(imresize(snip,10),[]);
    figure(4)
    surf(y, x, double(snip));
    title('Raw data');
end
