function [fits, intensity, CI, intensityError95] = fitGaussian3D(snip3D, initial_p, zStep, pixelSize, varargin)

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
%     single3DGaussian = @(p) p(1).*...
%         exp(-( p(2)*(mesh_x-p(3)).^2 + 2*p(4)*(mesh_x-p(3))*(mesh_y-p(5)) + p(6)*(mesh_y-p(5)).^2 ))...
%         + p(7) - double(snippet);


%     % Single 3D spherical gaussian function
%     single3DGaussian = @(p) p(1).*...
%         exp(-( ((mesh_x-p(2)).^2 + (mesh_y-p(3)).^2 + (mesh_z-p(4)).^2 ) / 2*p(5).^2) )...
%         + p(6) - snip3D;

% Single 3D generalized gaussian function

initial_p = double(initial_p);
single3DGaussian = gaussian3DForSpot(mesh_y,mesh_x, mesh_z, snip3D);

centroid_guess = [size(snip3D, 1)/2, size(snip3D, 2)/2, initial_p(2)];

initial_parameters =[initial_p(1), centroid_guess(1),centroid_guess(2), centroid_guess(3), ...
    initial_p(3)^(-2), 0, 0, initial_p(3)^(-2), 0, initial_p(3)^(-2),initial_p(4),...
    0, 0, 0,...
    ];

%%% p and fits: %%%
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
    ];
ub = [max(snip3D(:))*2, size(snip3D, 1)*1.5, size(snip3D, 2)*1.5, size(snip3D, 3)*1.5, xybound, covxybound, covzbound, xybound, covzbound, zbound, max(snip3D(:))/3,...
    1, 1, 1,...
    ];

[fits, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(single3DGaussian, ...
    initial_parameters,lb,ub, lsqOptions);

% ci63 = nlparci(fits, residual, 'jacobian', jacobian, 'alpha', .37); %63% confidence intervals for the fit. 
CI = nlparci(fits,residual,'jacobian',jacobian); %95% confidence intervals for the fit.



%%
%Intensity calculations
% vol = size(snip3D, 1)*size(snip3D,2)*size(snip3D,3);
gaussian = single3DGaussian(fits) + snip3D;

% intensity = sum(gaussian(:)) - vol*fits(end);


%%% fits and fits: %%%
%(1)amplitude (2) x (3) y (4) z (5)sigma x (6)sigma xy
%(7)sigma xz (8)sigma y (9)sigma yz (10) sigma z (11) offset
% covariance matrix = { {A D E} {D B F} {E F C} }
% formula: amplitude * sqrt(pi^3) / sqrt(det(cov))
%det cov = A*B*C - A*(F/2)^2 - B*(E/2)^2 -
%C*(D/2)^2 + D*E*F/4 
%where A = fits(5).  B = fits(8). C = fits(10). 
% D = fits(6). E = fits(7). F = fits(9)

syms amplitude A B C D E F
intCalc = amplitude * sqrt(pi^3) / sqrt( A*B*C - A*(F/2)^2 - B*(E/2)^2 -...
C*(D/2)^2 + D*E*F/4);


vals = [fits(1),fits(5), fits(8), fits(10), fits(6), fits(7), fits(9)]; 

errs = [CI(1,2) - fits(1), CI(5,2) - fits(5), CI(8,2) - fits(8), ...
    CI(10,2) - fits(10), CI(6,2) - fits(6),...
    CI(7,2) - fits(7), CI(9,2) - fits(9)]; %get errors from ci95 or ci63

[intensity, intensityError95] = PropError(intCalc, [amplitude A B C D E F], vals, errs);

%and because these values will never reach double precision:
intensity = single(intensity);
intensityError95 = single(intensityError95);
fits = single(fits);
CI = single(CI);

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
%     pause(.1)
end

if intensityError95 > intensity*2 || ~isreal(intensity)
    intensity = NaN;
    intensityError95 = NaN;
end
% 
% if ~isreal(intensity)
%     intensity = NaN;
%     intensityError95 = NaN;
% end

end

function [single3DGaussian, single3DGaussian2] = gaussian3DForSpot(y, x, z, snip3D)

    %%% p and fits: %%%
    %(1)amplitude (2) x (3) y (4) z (5)sigma x (6)sigma xy
    %(7)sigma xz (8)sigma y (9)sigma yz (10) sigma z (11) offset
    %(12)linear x offset (13) etc 14 etc

    
    single3DGaussian = @(p) p(1).*...
    exp(-( p(5).*(x-p(2)).^2 + p(6).*(x-p(2)).*(y-p(3)) + p(7).*(x-p(2)).*(z-p(4))+...
    p(8).*(y-p(3)).^2 + p(9).*(y-p(3)).*(z-p(4)) + p(10).*(z-p(4)).^2 ))...
    + p(11) +...
    p(12).*x + p(13).*y + p(14).*z...
    - double(snip3D);
% 
%     % %%% parameters and fits: %%%
%     %(1)amplitude (2) x (3) y (4) z 
%     %5. rhoxy 6. rhoxz 7. rhozy
%     %8. sigma x 9. sigma y  10. sigma z 
%     %(11) offset (12)linear x offset (13) linear y offset (14) linear z offset
% 
%     single3DGaussian2 = @(p) p(1).*exp(-(1/2)*(...
%      ( (p(8)^2*p(9)^2*p(10)^2) *...
%     ( 1 - 2*p(5)*p(6)*p(7)+ p(5) + p(6) + p(7) ) )^(-1) ...
%     .*...
%     ...
%     (x-p(2))*( (p(9)*p(10)*(p(7)-1))*(x-p(2))+...
%     ((p(8)*p(9)*p(10)^2*(p(5)-p(6)))*(y-p(3))+...
%     ((p(8)*p(9)^2*p(10)*(p(6)-p(5)))*(z-p(4)))+...
%     (y-p(3))*(((p(8)*p(9)*p(10)^2*(p(5)-p(6)))*(x-p(2))+...
%     (p(8)*p(10)*(p(6)-1))*(y-p(3))+...
%     K*(z-p(4)))+...
%     (z-p(4))*(I*(x-p(2))+...
%     K*(y-p(3))+...
%     (p(8)*p(9)*(p(5)-1))*(z-p(4)))...  
%     )...
%     + p(11) + p(12).*x + p(13).*y + p(14).*z...
%     ...
%     - double(snip3D);

        single3DGaussian2= [];

    
     
end


