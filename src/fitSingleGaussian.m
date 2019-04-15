function [fits, relative_errors, residual, confidence_intervals, GaussianIntensity, gaussian, mesh] = ...
    fitSingleGaussian(snippet, ~, ~, widthGuess, offsetGuess, show_status, graphicsHandles)

    % Fit Gaussians to the given locus within a snippet

    warning('off','MATLAB:singularMatrix')

    snippet = double(snippet);
    [mesh_y,mesh_x] = meshgrid(1:size(snippet,2), 1:size(snippet,1));

    %fits: [amplitude, x position, x width, y position, y width, offset, angle] 

    singleGaussian = gaussianForSpot(snippet);
%     [logGaussian, singleGaussian] = gaussianForSpot(snippet);
    
    %Define some more initial parameters for fitting

        initial_parameters = [max(snippet(:)), round(length(snippet)/2), widthGuess, round(length(snippet)/2), ...
            widthGuess,offsetGuess, 0];

    %Perform fitting
%     lsqOptions=optimset('Display','none',... %Inherited these options from Mikhail Tikhonov's FISH analysis
%     'maxfunevals',1000,...
%     'maxiter',1000); 
%     lsqOptions=optimset('Display','none', 'UseParallel', true);
lsqOptions=optimset('Display','none');

    lb_offset = 1/10; %this is empirical. corresponds to a weak background of 1 pixel per ten having a value of 1. 
    lb = [max(snippet(:))*.5, 0, 1, 0, 1,lb_offset, 0];
    ub = [max(snippet(:))*2, size(snippet, 1), size(snippet, 1), size(snippet, 2), size(snippet, 2), max(snippet(:)), 2*pi];

    
        [single_fit, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(singleGaussian, ...
            initial_parameters,lb,ub, lsqOptions);
        

        confidence_intervals = nlparci(single_fit,residual,'jacobian',jacobian);
        errors = zeros(1, length(single_fit));
        for i = 1:length(confidence_intervals)
            errors(i) = abs((abs(confidence_intervals(i, 1)) - abs(confidence_intervals(i, 2)))/2);
        end
        relative_errors = abs(errors./single_fit);
   
        
    fits = single_fit; 
    GaussianIntensity = sum(sum(singleGaussian(single_fit) + double(snippet) - single_fit(6)));

    %Display
    gaussian = singleGaussian(single_fit);
    mesh = {mesh_y, mesh_x};

    if show_status && ~isempty(graphicsHandles)
        gAx = graphicsHandles(4);
        snipAx = graphicsHandles(6);
        rawAx = graphicsHandles(8);
        surf(gAx, mesh_y, mesh_x, gaussian + snippet);
        title(gAx, 'Single Gaussian fit')
        snipBig = imresize(snippet,10);
        imshow(snipBig,[], 'Parent', snipAx);
        surf(rawAx,mesh_y, mesh_x, snippet);
        title(rawAx,'Raw data');
    end
    
end