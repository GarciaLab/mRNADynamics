function geometricEllipseParameters = fitEllipse(data_points)
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
geometricEllipseParameters = permute(geometricEllipseParameters, [3 4 1 2 5]);


% fprintf('Covariance matrix of geometric parameters: \n')
% geoCov =  compute_covariance_of_geometric_parameters(...
%                                theta_fastguaranteed, data_points);
%  fprintf('Standard deviation of geometric parameters: \n')                          
%  stds = sqrt(diag(geoCov)) 