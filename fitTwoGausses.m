function [f1, res1, residual, exitflag, output, lambda, jacobian, GaussianIntensity] = ...
    fitTwoGausses(snip, NeighborhoodSize, Threshold, WidthGuess, OffsetGuess, show)

% Find local maxima in snip and use that information to decide if fitting
% one or two gaussians. Also, use that information to define a reasonable 
% initial guess for the fitting initial parameters.

snip = CPsmooth(snip,'Gaussian Filter',1.3,0);
snip = double(snip);
[y,x] = meshgrid(1:size(snip,2), 1:size(snip,1));

singleGaussian = @(params) params(1).*exp((-1/2).*(((x-params(2))./params(3)).^2 ...
        + ((y-params(4))./params(5)).^2)) + params(6) - double(snip);

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
        centers(1,1), WidthGuess, OffsetGuess];

    [f1, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(singleGaussian, ...
        init_params,[0,0,0,0,0,0],[inf,inf,inf,inf,inf,inf]);
elseif size(centers,1) == 2
    init_params = [max(max(snip)), centers(1,2), WidthGuess, centers(1,1), ...
        WidthGuess, max(max(snip)), centers(2,2), WidthGuess, centers(2,1), ...
        WidthGuess, OffsetGuess];

    [f1, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(doubleGaussian, ...
        init_params,[0,0,0,0,0,0,0,0,0,0,0],[inf,inf,inf,inf,inf,inf,inf,inf,inf,inf,inf]);
else
    % TODO: I'm just making this so the script doesn't crash when it can't
    % find any maxima. But it's ugly and should be replaced with something
    % like assigning NaNs to all the return values.
    init_params = [2000, 10, 5, 10, 5,1000];
    [f1, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(singleGaussian, ...
        init_params,[0,0,0,0,0,0],[inf,inf,inf,inf,inf,inf]);
end

%Compute intensities by integrating over the Gaussian fit
% TODO: consider the case with two gaussians here

GaussianIntensity = sum(sum(singleGaussian(f1)+double(snip)-f1(end)));

if show
    figure(2)
    if size(centers,1) == 1
        surf(y, x, singleGaussian(f1) + double(snip));
        title('single gaussian')
    elseif size(centers,1) == 2
        surf(y, x, doubleGaussian(f1) + double(snip));
        title('double gaussian')
    end
    figure(3)
    imshow(imresize(snip,10),[]);
    figure(4)
    surf(y, x, double(snip));
end
