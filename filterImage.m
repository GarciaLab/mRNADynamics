function fint = filterImage(im, filterType, sigmas)

    dim = length(size(im));
    rad = 3; %rule of thumb is kernel size is 3x the Gaussian sigma
    f = [];
    fint = [];
    
    im = double(im);
    
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
    if exist('s','var')
        filterSize = rad*s;
        if ~mod(filterSize,2)
            filterSize = filterSize + 1;
        end
    end

    switch filterType        
        case 'Gaussian_blur'
            if dim == 2
                f = imgaussfilt(im,s);
            elseif dim == 3
                f = imgaussfilt3(im, s);
            end
        case 'Edges'
            if dim == 2
                f = edge(im, 'canny', s);
            elseif dim == 3
                f = canny(im, s);
            end
        case 'Difference_of_Gaussian'
            filterSize2 = rad*s2;
            if dim==2
                f = imgaussfilt(im, s1) - imgaussfilt(im, s2);
                f = padarray(f(filterSize2:end-filterSize2-1, filterSize2:end-filterSize2-1), [filterSize2,filterSize2]);
            elseif dim == 3
                f = imgaussfilt3(im, s1) - imgaussfilt3(im, s2);
            end
        case 'Laplacian'
            if dim==2
                h = fspecial('log', filterSize,s);
                f = imfilter(im, h);
            elseif dim==3
                G = gauss3D(s);
                [Gx,Gy] = gradient(G);
                [Gxx,~] = gradient(Gx);
                [~,Gyy] = gradient(Gy);
                LoG = Gxx + Gyy;
                f = fastConv3D(LoG, im);
            end
        case 'Mean'
            if dim == 2
                h = fspecial('average', s);
                f = imfilter(im, h);
            elseif dim == 3
                h = ones(filterSize, filterSize, filterSize);
                h = h./sum(sum(sum(h)));
                f = fastConv3D(h, im);
            end
        case 'Structure_smallest'
            if dim==2
                filterSize2 = rad*s1;
                G=fspecial('gauss',[round(filterSize2), round(filterSize2)], s1); 
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
                        S = [S11(p, q), S12(p, q); S21(p, q), S22(p,q)];
                        lambdas = eig(S);
                        lambdaMin = min(lambdas);
                        f(p,q) = lambdaMin;
                    end
                end
            elseif dim==3
                filterSize2 = rad*s1;
                G1 = gauss3D(s1);
                G = fastConv3D(G1, im);
                [Gx,Gy,Gz] = gradient(G);   
                %Compute Gaussian partial derivatives
                Dx = fastConv3D(Gx,im);
                Dy = fastConv3D(Gy, im);  
                Dz = fastConv3D(Gz, im);
                %Smooth elements of the structure tensor
                S11 = fastConv3D(G, Dx.^2);
                S12 = fastConv3D(G, Dx.*Dy);
                S13 = fastConv3D(G, Dx.*Dz);
                S21 = S12;
                S22 = fastConv3D(G, Dy.^2);
                S23 = fastConv3D(G, Dy.*Dz);
                S31 = S13;
                S32 = S23;
                S33 = fastConv3D(G, Dz.^2);
                %Make eigenimages from the structure tensors
                l = size(im, 1);
                m = size(im, 2);
                n = size(im, 3);
                f = zeros(l,m,n);
                parfor p = 1:l
                    for q = 1:m
                        for r = 1:n
                            S = [S11(p,q,r), S12(p,q,r), S13(p,q,r);...
                                S21(p,q,r), S22(p,q,r), S23(p,q,r);...
                                S31(p,q,r), S32(p,q,r), S33(p,q,r)];
                            f(p,q,r) = min(eig(S));
                        end
                    end
                end
            end
        case 'Structure_largest'
            if dim==2
                filterSize2 = rad*s1;
                G=fspecial('gauss',[round(filterSize2), round(filterSize2)], s1); 
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
                        S = [S11(p, q), S12(p, q); S21(p, q), S22(p,q)];
                        lambdas = eig(S);
                        lambdaMax = max(lambdas);
                        f(p,q) = lambdaMax;
                    end
                end
            elseif dim==3
                filterSize2 = rad*s1;
                G1 = gauss3D(s1);
                G = fastConv3D(G1, im);
                [Gx,Gy,Gz] = gradient(G);   
                %Compute Gaussian partial derivatives
                Dx = fastConv3D(Gx,im);
                Dy = fastConv3D(Gy, im);  
                Dz = fastConv3D(Gz, im);
                %Smooth elements of the structure tensor
                S11 = fastConv3D(G, Dx.^2);
                S12 = fastConv3D(G, Dx.*Dy);
                S13 = fastConv3D(G, Dx.*Dz);
                S21 = S12;
                S22 = fastConv3D(G, Dy.^2);
                S23 = fastConv3D(G, Dy.*Dz);
                S31 = S13;
                S32 = S23;
                S33 = fastConv3D(G, Dz.^2);
                %Make eigenimages from the structure tensors
                l = size(im, 1);
                m = size(im, 2);
                n = size(im, 3);
                f = zeros(l,m,n);
                parfor p = 1:l
                    for q = 1:m
                        for r = 1:n
                            S = [S11(p,q,r), S12(p,q,r), S13(p,q,r);...
                                S21(p,q,r), S22(p,q,r), S23(p,q,r);...
                                S31(p,q,r), S32(p,q,r), S33(p,q,r)];
                            f(p,q,r) = max(eig(S));
                        end
                    end
                end
            end
        case 'Median'
            if dim==2
                f = colfilt(im,[s s],'sliding',@median);
            elseif dim==3
%                 f = ordfilt3(im, 'med', filterSize); %i need to rewrite
%                 this algorithm because it doesn't work
            end
        case 'Maximum'
            if dim==2
                f = colfilt(im,[s s],'sliding',@max);
            elseif dim==3
%                 f = ordfilt3(im, 'max', filterSize); %i need to rewrite
%                 this algorithm because it doesn't work
            end
        case 'Variance'
            if dim==2
                f = colfilt(im,[s s],'sliding',@var);
            elseif dim==3
%                  f = ordfilt3(im, 'var', filterSize); %i need to rewrite
%                 this algorithm because it doesn't work
            end
        case 'Minimum'
            if dim==2
                f = colfilt(im,[s s],'sliding',@min);
            elseif dim==3
%                  f = ordfilt3(im, 'min', filterSize); %i need to rewrite
%                 this algorithm because it doesn't work
            end
        case 'Hessian_smallest'
            if dim==2
                G=fspecial('gauss',[round(filterSize), round(filterSize)], s); 
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
                        H = [H11(p, q), H12(p, q); H21(p, q), H22(p,q)];
                        lambdas = eig(H);
                        lambdaMin = min(lambdas);
                        f(p,q) = lambdaMin;
                    end
                end
            elseif dim ==3
                G1 = gauss3D(s);
                [Gx,Gy,Gz] = gradient(G1); 
                [Gxx, Gxy, Gxz] = gradient(Gx);
                [Gyy, Gyx, Gyz] = gradient(Gy);
                [Gzx, Gzy, Gzz] = gradient(Gz);
                %Compute elements of the Hessian matrix
                H11 = fastConv3D(Gxx,im);
                H12 = fastConv3D(Gxy,im);
                H13 = fastConv3D(Gxz, im);
                H21 = H12;
                H22 = fastConv3D(Gyy, im);
                H23 = fastConv3D(Gyz, im);
                H31 = H13;
                H32 = H23;
                H33 = fastConv3D(Gzz, im);
                %Make eigenimages from the Hessian
                l = size(im, 1);
                m = size(im, 2);
                n = size(im, 3);
                f = zeros(l,m,n);
                parfor p = 1:l
                    for q = 1:m
                        for r = 1:n
                            H = [H11(p, q, r), H12(p, q, r), H13(p, q, r);...
                                H21(p, q, r), H22(p,q, r), H23(p,q,r);...
                                H31(p,q,r), H32(p,q,r), H33(p,q,r)];
                            f(p,q,r) = min(eig(H));
                        end
                    end
                end              
            end
        case 'Hessian_largest'
            if dim==2 
                G=fspecial('gauss',[round(filterSize), round(filterSize)], s); 
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
                        H = [H11(p, q), H12(p, q); H21(p, q), H22(p,q)];
                        lambdas = eig(H);
                        lambdaMax = max(lambdas);
                        f(p,q) = lambdaMax;
                    end
                end
            elseif dim==3
                G1 = gauss3D(s);
                [Gx,Gy,Gz] = gradient(G1); 
                [Gxx, Gxy, Gxz] = gradient(Gx);
                [Gyy, Gyx, Gyz] = gradient(Gy);
                [Gzx, Gzy, Gzz] = gradient(Gz);
                %Compute elements of the Hessian matrix
                H11 = fastConv3D(Gxx, im);
                H12 = fastConv3D(Gxy, im);
                H13 = fastConv3D(Gxz, im);
                H21 = H12;
                H22 = fastConv3D(Gyy, im);
                H23 = fastConv3D(Gyz, im);
                H31 = H13;
                H32 = H23;
                H33 = fastConv3D(Gzz, im);
                %Make eigenimages from the Hessian
                l = size(im, 1);
                m = size(im, 2);
                n = size(im, 3);
                f = zeros(l,m,n);
                parfor p = 1:l
                    for q = 1:m
                        for r = 1:n
                            H = [H11(p, q, r), H12(p, q, r), H13(p, q, r);...
                                H21(p, q, r), H22(p,q, r), H23(p,q,r);...
                                H31(p,q,r), H32(p,q,r), H33(p,q,r)];
                            f(p,q,r) = max(eig(H));
                        end
                    end
                end              
            end
        otherwise
            %do nothing
    end
    
    if ~isempty(f)
        %feature rescaling
        fmin = min(min(min(f)));
        fmax = max(max(max(f)));
        fprime = (f - fmin)./(fmax - fmin);

        %reduce precision and recast for memory and speed
        ndigits = ceil(abs(log(std(fprime(:)))));
        fround = round(fprime, ndigits);
        fint = int16(fround*10^ndigits);
    end
    
end