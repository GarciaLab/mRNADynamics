function out = fitEllipse(data_points)
%wrapper function  for ellipse fitting from 
%https://github.com/zygmuntszpak/guaranteed-ellipse-fitting-with-a-confidence-region-and-an-uncertainty-measure
%output is [xCenter,yCenter, majAxis, minAxis, orientation (radians)]

%fprintf('Algebraic ellipse parameters of our method: \n')
[theta_fastguaranteed] = fast_guaranteed_ellipse_estimate(data_points);

% fprintf('Geometric ellipse parameters \n')
% fprintf('(majAxis, minAxis, xCenter,yCenter, orientation (radians)): \n')
geometricEllipseParameters = ...
            fromAlgebraicToGeometricParameters(theta_fastguaranteed);
 
%rearrange the output a bit
out = [];
out(1) = geometricEllipseParameters(3);
out(2) = geometricEllipseParameters(4);
out(3) = geometricEllipseParameters(1);
out(4) = geometricEllipseParameters(2);
out(5) = geometricEllipseParameters(5);


% fprintf('Covariance matrix of geometric parameters: \n')
% geoCov =  compute_covariance_of_geometric_parameters(...
%                                theta_fastguaranteed, data_points);
%  fprintf('Standard deviation of geometric parameters: \n')                          
%  stds = sqrt(diag(geoCov)) 