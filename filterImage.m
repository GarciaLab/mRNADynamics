function f = filterImage(im, filterType, sigmas)
    
    rad = 3; %rule of thumb is kernel size is 3x the Gaussian sigma
    f = [];
    
    q = length(sigmas);
    switch q        
        case 0
            %do nothing
        case 1
            s = str2double(sigmas{1});
        case 2
            s1 = str2double(sigmas{1});
            s2 = str2double(sigmas{2});
        otherwise
            s = str2double(sigmas{end});            
    end

    switch filterType        
        case 'Gaussian_blur'
            f = imgaussfilt(im, s);
        case 'Edges'
            f = edge(im, 'canny', s);
        case 'Difference_of_Gaussian'
            filterSize = rad*s2;            
            imDoG{s1, s2} = conv2(single(im), single(fspecial('gaussian',filterSize, s1) - fspecial('gaussian',filterSize, s2)),'same');
            f = padarray(imDoG{s1, s2}(filterSize:end-filterSize-1, filterSize:end-filterSize-1), [filterSize,filterSize]);
        case 'Laplacian'
            filterSize = rad*s;            
            h = fspecial('log', filterSize, s);
            f = imfilter(im, h);
        case 'Mean'
            h = fspecial('average', s);
            f = imfilter(im, h);
        case 'Structure_smallest'
            G=fspecial('gauss',[round(rad*s1), round(rad*s1)], s1); 
            [Gx,Gy] = gradient(G);   
            imStructSmall{s1, s2} = im;
            %Compute Gaussian partial derivatives
            Dx = conv2(Gx,im);
            Dy = conv2(Gy, im);                        
            %Smooth elements of the structure tensor
            S11 = conv2(G, Dx.^2);
            S12 = conv2(G, Dx.*Dy);
            S21 = S12;
            S22 = conv2(G, Dy.^2);
            %Make eigenimages from the structure tensors
            for p = 1:size(im, 1)
                for q = 1:size(im, 2)
                    S = [S11(p, q), S12(p, q), S21(p, q), S22(p,q)];
                    lambdas = eig(S);
                    lambdaMin = min(lambdas);
                    f(p,q) = lambdaMin;
                end
            end
        case 'Structure_largest'
            G=fspecial('gauss',[round(rad*s1), round(rad*s1)], s1); 
            [Gx,Gy] = gradient(G);              
            %Compute Gaussian partial derivatives
            Dx = conv2(Gx,im);
            Dy = conv2(Gy, im);                        
            %Smooth elements of the structure tensor
            S11 = conv2(G, Dx.^2);
            S12 = conv2(G, Dx.*Dy);
            S21 = S12;
            S22 = conv2(G, Dy.^2);
            %Make eigenimages from the structure tensors
            for p = 1:size(im, 1)
                for q = 1:size(im, 2)
                    S = [S11(p, q), S12(p, q), S21(p, q), S22(p,q)];
                    lambdas = eig(S);
                    lambdaMax = max(lambdas);
                    f(p,q) = lambdaMax;
                end
            end
        case 'Median'
            imMedian{s} = colfilt(im,[s s],'sliding',@median);
        case 'Maximum'
            imMax{s} = colfilt(im,[s s],'sliding',@max);
        case 'Maximum'
            imVariance{s} = colfilt(im,[s s],'sliding',@var);
        case 'Maximum'
            imMin{s} = colfilt(im,[s s],'sliding',@min);
        case 'Hessian_smallest'
            G=fspecial('gauss',[round(rad*s), round(rad*s)], s); 
            [Gx,Gy] = gradient(G);
            [Gxx, Gxy] = gradient(Gx);
            [Gyy, ~] = gradient(Gy);
            %Compute elements of the Hessian matrix
            H11 = conv2(Gxx, im);
            H12 = conv2(Gxy, im);
            H21 = H12;
            H22 = conv2(Gyy, im);
            %Make eigenimages from the Hessian
            for p = 1:size(im, 1)
                for q = 1:size(im, 2)
                    H = [H11(p, q), H12(p, q), H21(p, q), H22(p,q)];
                    lambdas = eig(H);
                    lambdaMin = min(lambdas);
                    f(p,q) = lambdaMin;
                end
            end
        case 'Hessian_largest'
            G=fspecial('gauss',[round(rad*s), round(rad*s)], s); 
            [Gx,Gy] = gradient(G);
            [Gxx, Gxy] = gradient(Gx);
            [Gyy, ~] = gradient(Gy);
            %Compute elements of the Hessian matrix
            H11 = conv2(Gxx, im);
            H12 = conv2(Gxy, im);
            H21 = H12;
            H22 = conv2(Gyy, im);
            %Make eigenimages from the Hessian
            for p = 1:size(im, 1)
                for q = 1:size(im, 2)
                    H = [H11(p, q), H12(p, q), H21(p, q), H22(p,q)];
                    lambdas = eig(H);
                    lambdaMax = max(lambdas);
                    f(p,q) = lambdaMax;
                end
            end
        otherwise
            %do nothing
    end
end