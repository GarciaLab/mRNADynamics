function [fits, relative_errors, residual, confidence_intervals, GaussianIntensity, dgaussian, mesh] = ...
    fitGaussians(snippet, neighborhoodSize, threshold, widthGuess, offsetGuess, show_status)

% Fit Gaussians to the given locus within a snippet

    snippet = double(snippet);
    [mesh_y,mesh_x] = meshgrid(1:size(snippet,2), 1:size(snippet,1));
    [smesh_y, smesh_x] = meshgrid(1:3, 1:3);
    singleGaussian = @(params) (params(1).*...
        exp(-(...
        (((cos(params(7)))^2 / (2*params(3)^2) ) + ((sin(params(7)))^2 / 2*params(5)^2))  .* (mesh_x-params(2)).^2 ...
        - 2*((-sin(2*params(7)) / (4*params(3)^2) ) + (sin(2*params(7)) / 4*params(5)^2)) .* (mesh_x-params(2)).*(mesh_y-params(4))...
        + (((sin(params(7)))^2 / (2*params(3)^2) ) + ((cos(params(7)))^2 / 2*params(5)^2)).* (mesh_y-params(4)).^2 ...
            )))...
        + params(6) - double(snippet);

    doubleGaussian = @(params) params(1).*exp((-1/2).*(((mesh_x-params(2))./params(3)).^2 ... 
            + ((mesh_y-params(4))./params(3)).^2)) ... 
            + params(5).*exp((-1/2).*(((mesh_x-params(6))./params(7)).^2  ...
            + ((mesh_y-params(8))./params(7)).^2))+ params(9) - double(snippet);
    
    neighborhoodSize = 2*floor(neighborhoodSize/2) + 1; %force odd neighborhood

    hLocalMax = vision.LocalMaximaFinder;
    hLocalMax.NeighborhoodSize = [neighborhoodSize, neighborhoodSize];
    hLocalMax.Threshold = threshold;
    centers = double(step(hLocalMax, snippet));    

    %initial parameters for single gaussian fit
    if ~isempty(centers)
        initial_parameters_s = [max(max(snippet)), centers(1,2), widthGuess, centers(1,1), ...
            widthGuess,offsetGuess, 0];
    else
        initial_parameters_s = [max(max(snippet)), round(length(snippet)/2), widthGuess, round(length(snippet)/2), ...
            widthGuess,offsetGuess, 0];
    end

    %initial parameters for double gaussian fit
    if size(centers,1)== 2
        initial_parameters = [max(max(snippet)), centers(1,2), widthGuess, centers(1,1), ...
                max(max(snippet)), centers(2,2), widthGuess, centers(2,1), ...
                offsetGuess];
    elseif size(centers, 1) == 1
        initial_parameters = [max(max(snippet)), centers(1,2), widthGuess, centers(1,1), ...
                max(max(snippet)), centers(1,2), widthGuess, centers(1,1), ...
                offsetGuess];
    else
        initial_parameters = [max(max(snippet)), round(size(snippet,1)/2), widthGuess, round(size(snippet,1)/2), ...
                max(max(snippet)), round(size(snippet,1)/2), widthGuess, round(size(snippet,1)/2), ...
                offsetGuess];
    end

    lsqOptions=optimset('Display','none',... %Inherited these options from Mikhail Tikhonov's FISH analysis
    'maxfunevals',10000,...
    'maxiter',10000);
    [single_fit, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(singleGaussian, ...
        initial_parameters_s,zeros(1,7),inf(1,7), lsqOptions);

    [double_fit, res1, residual, exitflag, output, lambda, jacobian] = lsqnonlin(doubleGaussian, ...
        initial_parameters,zeros(1,9),inf(1,9), lsqOptions);
    
    confidence_intervals = nlparci(double_fit,residual,'jacobian',jacobian);
    errors = zeros(1, length(double_fit));
    for i = 1:length(confidence_intervals)
        errors(i) = abs((abs(confidence_intervals(i, 1)) - abs(confidence_intervals(i, 2)))/2);
    end
    relative_errors = abs(errors./double_fit);
          
    %The idea behind this structure is that the first row is a single
    %gaussian fit and the following two rows are the sister chromatids of a
    %double gaussian fit
    
    fits = struct;
    fits(1).Amplitude = single_fit(1);
    fits(1).XPosition = single_fit(2);
    fits(1).YPosition = single_fit(4);
    fits(1).sigmaX = single_fit(3);
    fits(1).sigmaY = single_fit(5);
    fits(1).Offset = single_fit(6);
    fits(2).Amplitude = double_fit(1);
    fits(2).XPosition = double_fit(2);
    fits(2).YPosition = double_fit(4);
    fits(2).sigmaX = double_fit(3);
    fits(2).sigmaY = double_fit(3);
    fits(2).Offset = double_fit(9);
    fits(3).Amplitude = double_fit(5);
    fits(3).XPosition = double_fit(6);
    fits(3).YPosition = double_fit(8);
    fits(3).sigmaX = double_fit(7);
    fits(3).sigmaY = double_fit(7);
    fits(3).Offset = double_fit(9);
    fitfields = fieldnames(fits);
    for i = 1:length(fitfields)
        fits(4).(fitfields{i}) = [];
        fits(5).(fitfields{i}) = [];
    end
    
    im2 = zeros(size(snippet));
    for i = 2:size(snippet,1) - 1
        for j = 2:size(snippet,2) - 1
            if snippet(i,j) > snippet(i,j+1) &&...
               snippet(i,j) > snippet(i, j-1) &&...
               snippet(i,j) > snippet(i+1, j+1) &&...
               snippet(i,j) > snippet(i+1, j) &&...
               snippet(i,j) > snippet(i+1, j-1) &&...
               snippet(i,j) > snippet(i-1, j+1) &&...
               snippet(i,j) > snippet(i-1, j) &&...
               snippet(i,j) > snippet(i-1, j-1) &&...           
               mean([snippet(i,j+1), snippet(i,j-1),...
                   snippet(i-1,j+1), snippet(i-1,j), snippet(i-1, j+1),...
                   snippet(i+1,j+1), snippet(i+1,j), snippet(i+1, j+1)])...
                   >fits(1).Offset;
               im2(i,j) = 1;
            end
        end
    end
    poss = {};
    dist = [];
    im2 = ((im2.*snippet)>100);
    for i =1:size(snippet,1)
        for j = 1:size(snippet,2)
            if i==1 || j ==1 ||i ==size(snippet,1) || j==size(snippet,2)
                im2(i,j)=0;
            end
            if im2(i,j) ~= 0
                poss = [poss, [i,j]];
                dist = [dist,sqrt((i-ceil(size(snippet,1)/2))^2 + (j-ceil(size(snippet,2)/2))^2)];
            end
        end
    end    
 
%     while length(dist)>2
%         [~,i] = max(dist); 
%         im2(poss{i}(1), poss{i}(2)) = 0;
%         dist(i) = [];
%         poss(i) = [];
%     end

    center = 2; %middle pixel of sub-snippet
    if length(dist) > 0
        peak1 = snippet(poss{1}(1), poss{1}(2));
        d1snip = snippet(poss{1}(1)-1:poss{1}(1)+1, poss{1}(2)-1:poss{1}(2)+1);
        initial_parameters_d1 = [peak1-fits(1).Offset,widthGuess];
        doubleGaussian1 = @(params) params(1)*exp((-1/2).*(((smesh_x-center)./params(2)).^2 ... 
        + ((smesh_y-center)./params(2)).^2)) ... 
        + fits(1).Offset - double(d1snip); 
        d1_fit = lsqnonlin(doubleGaussian1, ...
        initial_parameters_d1,zeros(1,2),Inf(1,2),lsqOptions);
        fits(4).Amplitude = peak1;
        fits(4).XPosition = poss{1}(1);
        fits(4).YPosition = poss{1}(2);
        fits(4).sigmaX = d1_fit(2);
        fits(4).sigmaY = d1_fit(2);
        fits(4).Offset = fits(1).Offset;
        d1gaussian = doubleGaussian1(d1_fit);
        
        if length(dist) == 2
            peak2 = snippet(poss{2}(1), poss{2}(2));
            d2snip = snippet(poss{2}(1)-1:poss{2}(1)+1, poss{2}(2)-1:poss{2}(2)+1);
            initial_parameters_d2 = [peak2-fits(2).Offset,widthGuess];    
            doubleGaussian2 = @(params) params(1)*exp((-1/2).*(((smesh_x-center)./params(2)).^2 ... 
            + ((smesh_y-center)./params(2)).^2)) ... 
            + fits(1).Offset - double(d2snip);
            d2_fit = lsqnonlin(doubleGaussian2, ...
            initial_parameters_d2,zeros(1,2),inf(1,2), lsqOptions);
            fits(5).Amplitude = peak2;
            fits(5).XPosition = poss{2}(1);
            fits(5).YPosition = poss{2}(2);
            fits(5).sigmaX = d2_fit(2);
            fits(5).sigmaY = d2_fit(2);
            fits(5).Offset = fits(1).Offset;
            d2gaussian = doubleGaussian2(d2_fit);
        end
    end
    
    figure(4)
    imshow(imresize(im2,15),[])
    set(gcf,'units', 'normalized', 'position',[0.6, .2, .2, .2])
    
   
   if length(dist)>3
       1
   end
    
        
%Compute intensities by integrating over the Gaussian fit. Offset
%subtracted

    GaussianIntensity = sum(sum(doubleGaussian(double_fit) + double(snippet) - double_fit(end)));

    dgaussian = doubleGaussian(double_fit);
    sgaussian = singleGaussian(single_fit);
    mesh = {mesh_y, mesh_x};
    smesh = {smesh_y, smesh_x};


        if show_status
            figure(2);
            subplot(2,2,4)
            surf(mesh_y, mesh_x, dgaussian + snippet);
            title('Double Gaussian fits')
            set(gcf,'units', 'normalized', 'position',[0.01, .4, .5, .5]);
            zlim([0, max(max(snippet))])
            subplot(2,2,2)
            surf(mesh_y, mesh_x, sgaussian + snippet);
            title('Single Gaussian fits')
            zlim([0, max(max(snippet))])
            subplot(2,2,[1,3])
            surf(mesh_y, mesh_x, snippet);
            zlim([0, max(max(snippet))])
            title('Raw data');
%             set(gcf,'units', 'normalized', 'position',[.01, .1, .33, .33]);
            figure(3)
            snipBig = imresize(snippet,10);
            set(gcf,'units', 'normalized', 'position',[0.4, .2, .1, .1])
            imshow(snipBig,[]);
            figure(5)
            set(gcf,'units', 'normalized', 'position',[0.01, .1, .5, .5]);
            

            if length(dist) >0
                subplot(2,2,1)
                surf(smesh_y, smesh_x, d1snip)
                zlim([0, max(max(snippet))])
                subplot(2,2,2)
                surf(smesh_y, smesh_x, d1gaussian + d1snip);
                zlim([0, max(max(snippet))])

                if length(dist) == 2
                    subplot(2,2,3)
                    surf(smesh_y, smesh_x, d2snip)
                    subplot(2,2,4)
                    surf(smesh_y, smesh_x, d2gaussian + d2snip);
                    zlim([0, max(max(snippet))])
                else
                    subplot(2,2,3)
                    imshow(zeros(3))
                    subplot(2,2,4)
                    imshow(zeros(3))
                end
            else
                subplot(2,2,1)
                imshow(zeros(3))
                subplot(2,2,2)
                imshow(zeros(3))
                subplot(2,2,3)
                imshow(zeros(3))
                subplot(2,2,4)
                imshow(zeros(3))
            end
        end

end
