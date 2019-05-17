function [fits, intensity, ci95, intensityError95] = fitGaussian3D(snip3D, initial_params, zStep, pixelSize, varargin)

%%Fitting
displayFigures = 0;
for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'displayFigures')
        displayFigures = 1;
    end
end


    
% NL: Is this right?    
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

single3DGaussian = gaussian3DForSpot(mesh_y,mesh_x, mesh_z, snip3D);

centroid_guess = [size(snip3D, 1)/2, size(snip3D, 2)/2, initial_params(2)];


initial_parameters = [initial_params(1), centroid_guess(1),centroid_guess(2), centroid_guess(3), ...
    initial_params(3)^(-2), 0, 0, initial_params(3)^(-2), 0, initial_params(3)^(-2),initial_params(4),...
    0, 0, 0,...
    0, 0, 0, 0, 0, 0];

%%% params and fits: %%%
%(1)amplitude (2) x (3) y (4) z (5)sigma x (6)sigma xy
%(7)sigma xz (8)sigma y (9)sigma yz (10) sigma z (11) offset
%(12)linear x offset (13) etc 14 etc


%fitting options
lsqOptions=optimset('Display','none','maxfunevals',10000,...
    'maxiter',10000);

%constraints
covxybound = 30/pixelSize; %nm. empirically determined. 
covzbound = 200/pixelSize;
xybound = 100/pixelSize; %nm. empirically determined.
zbound = 800/pixelSize;%nm. empirically determined.
lb_offset = 1/10; %this is empirical. corresponds to a weak background of 1 pixel per ten having a value of .
lb = [max(snip3D(:)) / 10, 1 - size(snip3D, 1)*1.5, 1- size(snip3D, 2)*1.5, 1- size(snip3D, 3)*1.5, 0, -covxybound, -covzbound, 0, -covzbound, 0, lb_offset,...
    -1, -1, -1,...
    -1, -1, -1, -1, -1, -1];
ub = [max(snip3D(:))*2, size(snip3D, 1)*1.5, size(snip3D, 2)*1.5, size(snip3D, 3)*1.5, xybound, covxybound, covzbound, xybound, covzbound, zbound, max(snip3D(:))/3,...
    1, 1, 1,...
    1, 1, 1, 1, 1, 1];

[fits, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(single3DGaussian, ...
    initial_parameters,lb,ub, lsqOptions);

ci63 = nlparci(fits, residual, 'jacobian', jacobian, 'alpha', .37); %63% confidence intervals for the fit. 
ci95 = nlparci(fits,residual,'jacobian',jacobian); %95% confidence intervals for the fit.




%%
%Intensity calculations
% vol = size(snip3D, 1)*size(snip3D,2)*size(snip3D,3);
gaussian = single3DGaussian(fits) + snip3D;

% intensity = sum(gaussian(:)) - vol*fits(end);


%%% fits and fits: %%%
%(1)amplitude (2) x (3) y (4) z (5)sigma x (6)sigma xy
%(7)sigma xz (8)sigma y (9)sigma yz (10) sigma z (11) offset
% formula: amplitude * sqrt(pi^3) / sqrt(det(cov))
%det cov = A*B*C - A*(F/2)^2 - B*(E/2)^2 -
%C*(D/2)^2 + D*E*F/4 
%where A = fits(5).  B = fits(8). C = fits(10). 
% D = fits(6). E = fits(7). F = fits(9)

syms amplitude A B C D E F
intCalc = amplitude * sqrt(pi^3) / sqrt( A*B*C - A*(F/2)^2 - B*(E/2)^2 -...
C*(D/2)^2 + D*E*F/4);


vals = [fits(1),fits(5), fits(8), fits(10), fits(6), fits(7), fits(9)]; 

errs = [ci95(1,2) - fits(1), ci95(5,2) - fits(5), ci95(8,2) - fits(8), ...
    ci95(10,2) - fits(10), ci95(6,2) - fits(6),...
    ci95(7,2) - fits(7), ci95(9,2) - fits(9)]; %get errors from ci95 or ci63

[intensity, intensityError95] = PropError(intCalc, [amplitude A B C D E F], vals, errs);

if ~isreal(intensity)
    disp('uh oh complex');
end

%% All Plotting
if displayFigures
    
    close all
    
    figure(1)
    isosurface(gaussian, 3)
    axis tight
    
    % Generate arrays and grids/axis to plot
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
    
    %%%% Gaussian Fit Surface Plotting %%%%
    
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
    pause(.1)
end

if intensityError95 > intensity*2 || ~isreal(intensity)
    intensity = NaN;
    intensityError95 = NaN;
end

end

function single3DGaussian = gaussian3DForSpot(mesh_y, mesh_x, mesh_z, snip3D)

    %%% params and fits: %%%
    %(1)amplitude (2) x (3) y (4) z (5)sigma x (6)sigma xy
    %(7)sigma xz (8)sigma y (9)sigma yz (10) sigma z (11) offset
    %(12)linear x offset (13) etc 14 etc

    
    single3DGaussian = @(params) params(1).*...
    exp(-( params(5).*(mesh_x-params(2)).^2 + params(6).*(mesh_x-params(2)).*(mesh_y-params(3)) + params(7).*(mesh_x-params(2)).*(mesh_z-params(4))+...
    params(8).*(mesh_y-params(3)).^2 + params(9).*(mesh_y-params(3)).*(mesh_z-params(4)) + params(10).*(mesh_z-params(4)).^2 ))...
    + params(11) +...
    params(12).*mesh_x + params(13).*mesh_y + params(14).*mesh_z...
    + params(15).*(mesh_x.^2) + params(16).*(mesh_y.^2) + params(17).*(mesh_z.^2) +...
    + params(18).*(mesh_x.*mesh_z) + params(19).*mesh_x.*mesh_y + params(20).*mesh_y.*mesh_z...
    - double(snip3D);
     
end


