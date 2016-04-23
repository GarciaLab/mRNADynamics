%function [f1, res1, f2, res2] = fitGausses(snip)
function [f1, res1] = fitGausses(snip, show)

snip = double(snip);
[y,x] = meshgrid(1:size(snip,1), 1:size(snip,2));
gfit1 = @(params) params(1).*exp((-1/2).*(((x-params(2))./params(3)).^2 + ((y-params(4))./params(5)).^2)) + params(6) - double(snip);
init_params = [2000, 10, 5, 10, 5,1000];
[f1,res1] = lsqnonlin(gfit1, init_params,[0,0,0,0,0,0],[inf,inf,inf,inf,inf,inf]);
if show
    figure(2)
    surf(y, x, gfit1(f1) + snip);
    title('single gaussian')
    figure(3)
    imshow(imresize(snip,10),[]);
end

% %Fit two 2D Gaussians
% options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective',...
%     'MaxFunctionEvaluations',10000,'MaxIterations',10000,'UseParallel',1);
% gfit2 = @(params) params(1).*exp((-1/2).*(((x-params(2))./params(3)).^2 + ((y-params(4))./params(5)).^2)) + params(6).*exp((-1/2).*(((x-params(7))./params(8)).^2 + ((y-params(9))./params(10)).^2))+ params(11) - double(snip);
% init_params = [1000, 30, 5, 30, 5,1000, 30, 5, 30, 5, 1000];
% [f2,res2] = lsqnonlin(gfit2, init_params,[500,0,1,0,1,500,0,1,0,1,500],[5000,size(snip,1),size(snip,1),size(snip,2),size(snip,2),5000,size(snip,1),size(snip,1),size(snip,2),size(snip,2),5000],options);
% figure(4)
% surf(x, y, gfit2(f2) + double(snip));
% title('double gaussian')
