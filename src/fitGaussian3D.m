function [fits, intensity] = fitGaussian3D(snip3D, initial_params)

    [mesh_y,mesh_x, mesh_z] = meshgrid(1:size(snip3D,2), 1:size(snip3D,1), 1:size(snip3D, 3));

%     % Single 2D generalized gaussian function
%     single3DGaussian = @(params) params(1).*...
%         exp(-( params(2)*(mesh_x-params(3)).^2 + 2*params(4)*(mesh_x-params(3))*(mesh_y-params(5)) + params(6)*(mesh_y-params(5)).^2 ))...
%         + params(7) - double(snippet);


%     % Single 3D spherical gaussian function
%     single3DGaussian = @(params) params(1).*...
%         exp(-( ((mesh_x-params(2)).^2 + (mesh_y-params(3)).^2 + (mesh_z-params(4)).^2 ) / 2*params(5).^2) )...
%         + params(6) - snip3D;

    % Single 3D generalized gaussian function
    single3DGaussian = @(params) params(1).*...
        exp(-( params(5).*(mesh_x-params(2)).^2 + params(6).*(mesh_x-params(2)).*(mesh_y-params(3)) + params(7).*(mesh_x-params(2)).*(mesh_z-params(4))+...
        params(8).*(mesh_y-params(3)).^2 + params(9).*(mesh_y-params(3)).*(mesh_z-params(4)) + params(10).*(mesh_z-params(4)).^2 ))...
        + params(11) - snip3D;
    
centroid_guess = [size(snip3D, 1)/2, size(snip3D, 2)/2, initial_params(4)];

% initial_parameters = [initial_params(1), centroid_guess(1),centroid_guess(2), centroid_guess(3), ...
%         initial_params(5),initial_params(6)];   
initial_parameters = [initial_params(1), centroid_guess(1),centroid_guess(2), centroid_guess(3), ...
        1, 1, 1, 1, 1, 1,initial_params(6)];    

    
%fitting options
lsqOptions=optimset('Display','none',... %Inherited these options from Mikhail Tikhonov's FISH analysis
'maxfunevals',10000,...
'maxiter',10000); 

%constraints
% lb = [min(min(min(snip3D))), 1, 1, 1, 0, 0];
% ub = [inf, size(snip3D, 1), size(snip3D, 2), size(snip3D, 3), size(snip3D, 1)/2, max(max(max(snip3D)))];

lb_offset = 1/10; %this is empirical. corresponds to a weak background of 1 pixel per ten having a value of 1. 
lb = [max(max(max(snip3D))) / 2, 1, 1, 1, -inf, -inf, -inf, -inf, -inf, -inf, lb_offset];
ub = [inf, size(snip3D, 1), size(snip3D, 2), size(snip3D, 3), inf, inf, inf, inf, inf, inf,max(max(max(snip3D)))];

[fits, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(single3DGaussian, ...
    initial_parameters,lb,ub, lsqOptions);


% %Display
% m = max(snip3D(:));
% figure(1)
% [y,x] = meshgrid(1:size(snip3D,2), 1:size(snip3D,1));
gaussian = single3DGaussian(fits);
% proj = max(gaussian + snip3D, [], 3);
% surf(y,x, proj)
% zlim([0, m]);
% figure(2)
% surf(y, x, max(snip3D, [], 3))
% zlim([0, m]);
% 
vol = size(snip3D, 1)*size(snip3D,2)*size(snip3D,3);
intensity = sum(gaussian(:) + snip3D(:)) - vol*fits(end);
% 
% figure(2);
% surf(mesh_y, mesh_x, gaussian + snippet);
% title('Single Gaussian fit')
% set(gcf,'units', 'normalized', 'position',[0.01, .55, .33, .33]);
% figure(3)
% snipBig = imresize(snippet,10);
% set(gcf,'units', 'normalized', 'position',[0.4, .2, .1, .1])
% imshow(snipBig,[]);
% figure(4);
% surf(mesh_y, mesh_x, snippet);
% title('Raw data');
% set(gcf,'units', 'normalized', 'position',[.01, .1, .33, .33]);

    
end


