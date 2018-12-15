function [fits, intensity] = fitGaussian3D(snip3D, initial_params, zstep, varargin)

%%Fitting
displayFigures = 0;
for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'displayFigures')
        displayFigures = 1;
    end
end

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
    
%%% params and fits: %%%
%(1)amplitude (2) x (3) y (4) z (5)sigma x (6)sigma xy
%(7)sigma xz (8)sigma y (9)sigma yz (10) sigma z (11) offset 

   
%fitting options
lsqOptions=optimset('Display','none',... %Inherited these options from Mikhail Tikhonov's FISH analysis
'maxfunevals',10000,...
'maxiter',10000); 

%constraints
% lb = [min(min(min(snip3D))), 1, 1, 1, 0, 0];
% ub = [inf, size(snip3D, 1), size(snip3D, 2), size(snip3D, 3), size(snip3D, 1)/2, max(max(max(snip3D)))];

lb_offset = 1/10; %this is empirical. corresponds to a weak background of 1 pixel per ten having a value of 1. 
lb = [max(max(max(snip3D))) / 2, 1 - size(snip3D, 1)*1.5, 1- size(snip3D, 2)*1.5, 1- size(snip3D, 3)*1.5, -inf, -inf, -inf, -inf, -inf, -inf, lb_offset];
% ub = [inf, size(snip3D, 1)*1.5, size(snip3D, 2)*1.5, size(snip3D, 3)*1.5, inf, inf, inf, inf, inf, inf,max(max(max(snip3D)))];
ub = [inf, size(snip3D, 1)*1.5, size(snip3D, 2)*1.5, size(snip3D, 3)*1.5, inf, inf, inf, ceil(0.8/zstep), inf, inf,max(max(max(snip3D)))];
[fits, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(single3DGaussian, ...
    initial_parameters,lb,ub, lsqOptions);

%Display

close all
vol = size(snip3D, 1)*size(snip3D,2)*size(snip3D,3);
gaussian = single3DGaussian(fits) + snip3D;
intensity = sum(gaussian(:)) - vol*fits(end);

%% All Plotting
if displayFigures
    figure(1)
    isosurface(gaussian, 3)
%     p = patch(isosurface(gaussian, 3));
%     p.FaceColor = 'none';
%     p.EdgeColor = 'green';
    axis tight

    % Generate arrays and grids/axis to plot %
    m = max(snip3D(:));
    f2 = figure(2);
    [yz,xz] = meshgrid(1:size(snip3D,2), 1:size(snip3D,1));
    [yx,zx] = meshgrid(1:size(snip3D,2), 1:size(snip3D,3));
    [xy,zy] = meshgrid(1:size(snip3D,1), 1:size(snip3D,3));
    projz = max(gaussian, [], 3);
    projxFit = [];
    projyFit = [];
    projxRaw = [];
    projyRaw = [];
    for z = 1:size(zx,1)
        maxxFit = max(gaussian, [], 1);
        maxyFit = max(gaussian, [], 2);
        maxxRaw = max(snip3D, [], 1);
        maxyRaw = max(snip3D, [], 2);
        projxFit(z,:) = maxxFit(:,:,z);
        projyFit(z,:) = maxyFit(:,:,z);
        projxRaw(z,:) = maxxRaw(:,:,z);
        projyRaw(z,:) = maxyRaw(:,:,z);
    end
    axesxyFit = axes(f2);
    axesyzFit = axes(f2);
    axesxzFit = axes(f2);
    axesxyRaw = axes(f2);
    axesyzRaw = axes(f2);
    axesxzRaw = axes(f2);
    
    %%%% Gaussian Fit Plotting %%%%
    subplot(3,2,2,axesxyFit)
    surf(yz,xz, projz)
    title(axesxyFit,{'Gaussian fit in X,Y';'(Z projected)'})
    xlabel(axesxyFit,'x-axis')
    ylabel(axesxyFit,'y-axis')
    
    subplot(3,2,4,axesyzFit)
    surf(yx,zx, projxFit)
    title(axesyzFit,{'Gaussian fit in Y,Z';'(X projected)'})
    xlabel(axesyzFit,'y-axis')
    ylabel(axesyzFit,'z-axis')
    
    subplot(3,2,6,axesxzFit)
    surf(xy,zy, projyFit)
    title(axesxzFit,{'Gaussian fit in X,Z';'(Y projected)'})
    xlabel(axesxzFit,'x-axis')
    ylabel(axesxzFit,'z-axis')

    %%%% Raw Data Plotting %%%%
    subplot(3,2,1,axesxyRaw)
    surf(yz, xz, max(snip3D, [], 3))
    title(axesxyRaw,{'Raw spot data in X,Y';'(Z projected)'})
    xlabel(axesxyRaw,'x-axis')
    ylabel(axesxyRaw,'y-axis')
    
    subplot(3,2,3,axesyzRaw)
    surf(yx, zx, projxRaw)
    title(axesyzRaw,{'Raw spot data in Y,Z';'(X projected)'})
    xlabel(axesyzRaw,'y-axis')
    ylabel(axesyzRaw,'z-axis')
    
    subplot(3,2,5,axesxzRaw)
    surf(xy, zy, projyRaw)
    title(axesxzRaw,{'Raw spot data in X,Z';'(Y projected)'})
    xlabel(axesxzRaw,'x-axis')
    ylabel(axesxzRaw,'z-axis')
    
    zlim([0, m]);
    pause(1)
end
    
%     figure(2);
%     surf(mesh_y, mesh_x, gaussian + snippet);
%     title('Single Gaussian fit')
%     set(gcf,'units', 'normalized', 'position',[0.01, .55, .33, .33]);
%     figure(3)
%     snipBig = imresize(snippet,10);
%     set(gcf,'units', 'normalized', 'position',[0.4, .2, .1, .1])
%     imshow(snipBig,[]);
%     figure(4);
%     surf(mesh_y, mesh_x, snippet);
%     title('Raw data');
%     set(gcf,'units', 'normalized', 'position',[.01, .1, .33, .33]); s

    
end


