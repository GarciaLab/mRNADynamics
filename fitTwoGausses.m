function [f1, res1, residual, exitflag, output, lambda, jacobian] = ...
    fitTwoGausses(snip, NeighborhoodSize, Threshold, WidthGuess, OffsetGuess, show)

snip = double(snip);
[y,x] = meshgrid(1:size(snip,2), 1:size(snip,1));

singleGaussian = @(params) params(1).*exp((-1/2).*(((x-params(2))./params(3)).^2 ...
        + ((y-params(4))./params(5)).^2)) + params(6) - double(snip);

doubleGaussian = @(params) params(1).*exp((-1/2).*(((x-params(2))./params(3)).^2 ... 
        + ((y-params(4))./params(5)).^2)) ... 
        + params(6).*exp((-1/2).*(((x-params(7))./params(8)).^2  ...
        + ((y-params(9))./params(10)).^2))+ params(11) - double(snip);

hLocalMax = vision.LocalMaximaFinder;
hLocalMax.NeighborhoodSize = NeighborhoodSize;
hLocalMax.Threshold = Threshold;
centers = double(step(hLocalMax, snip));

if size(centers,1) == 1
    init_params = [max(max(snip)), centers(1,2), WidthGuess, ... 
        centers(1,1), WidthGuess, OffsetGuess];

    [f1, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(singleGaussian, ...
        init_params,[0,0,0,0,0,0],[inf,inf,inf,inf,inf,inf]);
else
    init_params = [max(max(snip)), centers(1,2), WidthGuess, centers(1,1), ...
        WidthGuess, max(max(snip)), centers(2,2), WidthGuess, centers(2,1), ...
        WidthGuess, OffsetGuess];

    [f1, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(doubleGaussian, ...
        init_params,[0,0,0,0,0,0,0,0,0,0,0],[inf,inf,inf,inf,inf,inf,inf,inf,inf,inf,inf]);
end

if show
    figure(2)
    if size(centers,1) == 1
        surf(y, x, singleGaussian(f1) + double(snip));
        title('single gaussian')
    else
        surf(y, x, doubleGaussian(f1) + double(snip));
        title('double gaussian')
    end
    figure(3)
    imshow(imresize(snip,10),[]);
    figure(4)
    surf(y, x, double(snip));
end
