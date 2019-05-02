function fint = filterImage(im, filterType, sigmas, varargin)

    dim = length(size(im));
    rad = 3; %rule of thumb is kernel size is 3x the Gaussian sigma
    f = [];
    fint = [];
        
    im = double(im);
    
    %convert string sigmas to doubles
    if ~isempty(sigmas)
        if ischar(sigmas{1})
            for i=1:length(sigmas)
                sigmas{i}=str2double(sigmas{i});
            end
        end
    end
    
    q = length(sigmas);
    switch q        
        case 0
            %do nothing
        case 1
            s1 = sigmas{1};
        case 2
            s1 = sigmas{1};
            s2 = sigmas{2};
        otherwise
            s1 = str2double(sigmas{end});            
    end
    if exist('s1','var')
        filterSize = round(rad*s1);
    end
    
    if ~isempty(varargin)
        customSize = varargin{1};
        filterSize = round(customSize);
    end
    
    if ~mod(filterSize,2)
        filterSize = filterSize + 1;
    end

    switch filterType  
        case 'Identity'
            f=im;
        case 'Gaussian_blur'
            if dim == 2
                f = imgaussfilt(im,s1);
            elseif dim == 3
                f = imgaussfilt3(im, s1);
            end
        %AR 10/3/18 these next two have not yet been implemented. 
        case 'bright_spot_psf'
              %21 x 21
              spotfilt = [0,0.150000000000000,0.0500000000000000,0,0.0500000000000000,0,0.100000000000000,0.100000000000000,0,0.100000000000000,0.0500000000000000,0.0500000000000000,0,0.100000000000000,0,0,0,0,0.0500000000000000,0.100000000000000,0.100000000000000;0.0500000000000000,0.0500000000000000,0,0,0.100000000000000,0,0.0500000000000000,0,0,0,0.0500000000000000,0.0500000000000000,0,0,0.100000000000000,0.150000000000000,0,0,0.200000000000000,0.100000000000000,0.100000000000000;0.0500000000000000,0.0500000000000000,0,0.0500000000000000,0.0500000000000000,0,0,0,0.100000000000000,0.0500000000000000,0,0.0500000000000000,0,0.0500000000000000,0,0,0.200000000000000,0,0.100000000000000,0.100000000000000,0;0.0500000000000000,0.0500000000000000,0,0.0500000000000000,0,0.100000000000000,0,0.150000000000000,0,0.0500000000000000,0.0500000000000000,0,0.0500000000000000,0.0500000000000000,0,0.0500000000000000,0,0.100000000000000,0,0.100000000000000,0.0500000000000000;0,0,0.100000000000000,0.150000000000000,0.0500000000000000,0.100000000000000,0,0.150000000000000,0.0500000000000000,0.100000000000000,0.100000000000000,0.0500000000000000,0.100000000000000,0.150000000000000,0,0.150000000000000,0.0500000000000000,0,0,0.0500000000000000,0;0,0.100000000000000,0.200000000000000,0.100000000000000,0,0,0.100000000000000,0.100000000000000,0.200000000000000,0.300000000000000,0.150000000000000,0,0.200000000000000,0,0.150000000000000,0,0.0500000000000000,0,0.0500000000000000,0,0.100000000000000;0,0.0500000000000000,0,0,0,0.0500000000000000,0,0.150000000000000,0.100000000000000,0.150000000000000,0.150000000000000,0.200000000000000,0.250000000000000,0.200000000000000,0.250000000000000,0.0500000000000000,0.150000000000000,0,0.0500000000000000,0.0500000000000000,0.0500000000000000;0.0500000000000000,0.0500000000000000,0.100000000000000,0.0500000000000000,0.0500000000000000,0.0500000000000000,0.200000000000000,0.0500000000000000,0.200000000000000,0.250000000000000,0.300000000000000,0.300000000000000,0.0500000000000000,0.350000000000000,0.200000000000000,0.150000000000000,0.150000000000000,0.0500000000000000,0.0500000000000000,0.0500000000000000,0.0500000000000000;0.0500000000000000,0.100000000000000,0.0500000000000000,0.0500000000000000,0.150000000000000,0.150000000000000,0.550000000000000,0.200000000000000,0.250000000000000,0.650000000000000,0.550000000000000,0.500000000000000,0.550000000000000,0.250000000000000,0.350000000000000,0.0500000000000000,0.100000000000000,0.200000000000000,0.150000000000000,0.150000000000000,0;0.0500000000000000,0,0.100000000000000,0.200000000000000,0.200000000000000,0.250000000000000,0.0500000000000000,0.150000000000000,0.300000000000000,1,0.500000000000000,0.650000000000000,0.750000000000000,0.450000000000000,0.350000000000000,0.150000000000000,0.150000000000000,0.100000000000000,0.100000000000000,0,0.0500000000000000;0,0,0.0500000000000000,0,0.250000000000000,0.100000000000000,0.400000000000000,0.200000000000000,0.600000000000000,0.500000000000000,0.750000000000000,0.700000000000000,0.550000000000000,0.400000000000000,0.250000000000000,0.0500000000000000,0,0.150000000000000,0.200000000000000,0.0500000000000000,0.0500000000000000;0.0500000000000000,0.0500000000000000,0.0500000000000000,0.100000000000000,0.150000000000000,0.150000000000000,0.500000000000000,0.500000000000000,0.400000000000000,0.550000000000000,0.600000000000000,0.600000000000000,0.650000000000000,0.500000000000000,0.150000000000000,0.150000000000000,0.150000000000000,0.0500000000000000,0.100000000000000,0,0.100000000000000;0,0.200000000000000,0,0.0500000000000000,0.0500000000000000,0.250000000000000,0.150000000000000,0.300000000000000,0.450000000000000,0.400000000000000,0.600000000000000,0.550000000000000,0.350000000000000,0.300000000000000,0.200000000000000,0.100000000000000,0.0500000000000000,0.0500000000000000,0.100000000000000,0.0500000000000000,0;0.100000000000000,0.100000000000000,0.0500000000000000,0.100000000000000,0.100000000000000,0.300000000000000,0.0500000000000000,0.350000000000000,0.250000000000000,0.750000000000000,0.450000000000000,0.500000000000000,0.550000000000000,0.250000000000000,0.150000000000000,0.0500000000000000,0.0500000000000000,0.200000000000000,0,0.100000000000000,0.0500000000000000;0.0500000000000000,0.0500000000000000,0.100000000000000,0.0500000000000000,0.250000000000000,0.150000000000000,0.250000000000000,0.150000000000000,0.200000000000000,0.350000000000000,0.250000000000000,0.300000000000000,0.150000000000000,0.100000000000000,0,0.0500000000000000,0.100000000000000,0.150000000000000,0.100000000000000,0,0;0,0.100000000000000,0.0500000000000000,0.0500000000000000,0.0500000000000000,0,0.0500000000000000,0.200000000000000,0.200000000000000,0.100000000000000,0.300000000000000,0.100000000000000,0.200000000000000,0,0.150000000000000,0.0500000000000000,0.100000000000000,0,0.100000000000000,0.150000000000000,0;0,0.100000000000000,0.150000000000000,0,0.0500000000000000,0.0500000000000000,0.0500000000000000,0.250000000000000,0.0500000000000000,0.100000000000000,0.100000000000000,0.0500000000000000,0.0500000000000000,0,0,0,0.0500000000000000,0.0500000000000000,0.0500000000000000,0.100000000000000,0.0500000000000000;0,0,0.0500000000000000,0.0500000000000000,0.100000000000000,0.100000000000000,0.100000000000000,0.150000000000000,0.0500000000000000,0,0.100000000000000,0.0500000000000000,0.0500000000000000,0.150000000000000,0.0500000000000000,0,0.150000000000000,0.100000000000000,0,0.200000000000000,0.150000000000000;0.0500000000000000,0.150000000000000,0.100000000000000,0.0500000000000000,0.100000000000000,0.100000000000000,0,0.0500000000000000,0.150000000000000,0.100000000000000,0.100000000000000,0.0500000000000000,0,0.100000000000000,0.0500000000000000,0.0500000000000000,0.0500000000000000,0,0.0500000000000000,0,0.150000000000000;0,0,0.0500000000000000,0,0.0500000000000000,0.150000000000000,0.100000000000000,0.0500000000000000,0.100000000000000,0.0500000000000000,0.150000000000000,0.0500000000000000,0,0.150000000000000,0.0500000000000000,0.100000000000000,0,0.0500000000000000,0.0500000000000000,0.0500000000000000,0.0500000000000000;0.0500000000000000,0.0500000000000000,0.0500000000000000,0,0,0,0.0500000000000000,0.150000000000000,0.0500000000000000,0,0.100000000000000,0.100000000000000,0.0500000000000000,0.100000000000000,0,0,0.0500000000000000,0.0500000000000000,0.100000000000000,0.100000000000000,0.0500000000000000];
              f = conv2(im, spotfilt, 'same');
        case 'dim_spot_psf'
%             f = spotfilt(im,s1, dim, 'small');
              f = im;
        case 'Edges'
            if dim == 2
                f = edge(im, 'canny', s1);
            elseif dim == 3
                f = canny(im, s1);
            end
        case 'Difference_of_Gaussian'
            if s2 < s1 
                error('DoG filter requires sigma 1 < sigma 2')
            end
            %assumes sigma 2 > sigma 1
            if dim==2
                f = conv2(im, fspecial('gaussian',filterSize, s1) - fspecial('gaussian',filterSize, s2),'same');
            elseif dim == 3
                 f = fastConv3D(gauss3D(s1,3*s2+1) - gauss3D(s2), im);
            end
        case 'Laplacian'
            if dim==2
                h = fspecial('log', filterSize,s1);
                f = imfilter(im, h);
            elseif dim==3
                G = gauss3D(s1);
                [Gx,Gy] = gradient(G);
                [Gxx,~] = gradient(Gx);
                [~,Gyy] = gradient(Gy);
                LoG = Gxx + Gyy;
                f = fastConv3D(LoG, im);
            end
        case 'Mean'
            if dim == 2
                h = fspecial('average', s1);
                f = imfilter(im, h);
            elseif dim == 3
                h = ones(filterSize, filterSize, filterSize);
                h = h./sum(sum(sum(h)));
                f = fastConv3D(h, im);
            end
        case 'Structure_smallest'
            if dim==2
                G=fspecial('gauss',[filterSize, filterSize], s1); 
                [Gx,Gy] = gradient(G);   
                %Compute Gaussian partial derivatives
                Dx = conv2(im, Gx,'same');
                Dy = conv2(im, Gy, 'same');
                %Smooth elements of the structure tensor
                S11 = conv2(Dx.^2,G,'same');
                S12 = conv2(Dx.*Dy,G,'same');
                S21 = S12;
                S22 = conv2(Dy.^2,G,'same');
                %Make eigenimages from the structure tensors
                for p = 1:size(im, 1)
                    for q = 1:size(im, 2)
                        S = [S11(p, q), S12(p, q); S21(p, q), S22(p,q)]; 
                        f(p,q) = min(eig(S));
                    end
                end
            elseif dim==3
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
                G=fspecial('gauss',[filterSize, filterSize], s1); 
                [Gx,Gy] = gradient(G);              
                %Compute Gaussian partial derivatives
                Dx = conv2(im,Gx,'same');
                Dy = conv2(im,Gy,'same');                        
                %Smooth elements of the structure tensor
                S11 = conv2(Dx.^2,G,'same');
                S12 = conv2(Dx.*Dy,G,'same');
                S21 = S12;
                S22 = conv2(Dy.^2,'same');
                %Make eigenimages from the structure tensors
                for p = 1:size(im, 1)
                    for q = 1:size(im, 2)
                        S = [S11(p, q), S12(p, q); S21(p, q), S22(p,q)];
                        f(p,q) = max(eig(S));
                    end
                end
            elseif dim==3
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
                f = imgaussfilt(im,s1);
                f = ordfilt2(f,ceil(filterSize*filterSize/2),ones(filterSize,filterSize));
            elseif dim==3
%                 f = ordfilt3(im, 'med', filterSize); %i need to rewrite
%                 this algorithm because it doesn't work
                f = im;
            end
        case 'Maximum'
            if dim==2
                f = imgaussfilt(im,s1);
                se = strel('disk',ceil(filterSize/2));
                f = imdilate(f,se);
%                 f = ordfilt2(f,(filterSize*filterSize),ones(filterSize,filterSize));
%                 f = imgaussfilt(f,s);
            elseif dim==3
%                 f = ordfilt3(im, 'max', filterSize); %i need to rewrite
%                 this algorithm because it doesn't work
                f = im;
            end
        case 'Variance'
            if dim==2
                f = imgaussfilt(im,s1);
                f = colfilt(f,[filterSize filterSize],'sliding',@variance);
            elseif dim==3
%                  f = ordfilt3(im, 'var', filterSize); %i need to rewrite
%                 this algorithm because it doesn't work
                f = im;
            end
        case 'Minimum'
            if dim==2
                f = imgaussfilt(im,s1);
                se = strel('disk',ceil(filterSize/2));
                f = imerode(f,se);
            elseif dim==3
%                  f = ordfilt3(im, 'min', filterSize); %i need to rewrite
%                 this algorithm because it doesn't work
                f = im;                       
            end
        case 'Std'
            if dim==2
                f = imgaussfilt(im,s1);
                f = stdfilt(f);
            elseif dim==3
%                  f = ordfilt3(im, 'min', filterSize); %i need to rewrite
%                 this algorithm because it doesn't work
                f = im;                       
            end
        case 'Hessian_smallest'
            if dim==2
                G=fspecial('gauss',[filterSize, filterSize], s1); 
                [Gx,Gy] = gradient(G);
                [Gxx, Gxy] = gradient(Gx);
                [Gyy, ~] = gradient(Gy);
                %Compute elements of the Hessian matrix
                H11 = conv2(im,Gxx,'same');
                H12 = conv2(im,Gxy,'same');
                H21 = H12;
                H22 = conv2(im,Gyy,'same');
                %Make eigenimages from the Hessian
                for p = 1:size(im, 1)
                    for q = 1:size(im, 2)
                        H = [H11(p, q), H12(p, q); H21(p, q), H22(p,q)];
                        f(p,q) = min(eig(H));
                    end
                end
            elseif dim ==3
                G1 = gauss3D(s1);
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
                G=fspecial('gauss',[filterSize, filterSize], s1); 
                [Gx,Gy] = gradient(G);
                [Gxx, Gxy] = gradient(Gx);
                [Gyy, ~] = gradient(Gy);
                %Compute elements of the Hessian matrix
                H11 = conv2(im,Gxx,'same');
                H12 = conv2(im,Gxy,'same');
                H21 = H12;
                H22 = conv2(im,Gyy,'same');
                %Make eigenimages from the Hessian
                for p = 1:size(im, 1)
                    for q = 1:size(im, 2)
                        H = [H11(p, q), H12(p, q); H21(p, q), H22(p,q)];
                        f(p,q) = max(eig(H));
                    end
                end
            elseif dim==3
                G1 = gauss3D(s1);
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
  
    %this ended up altering the image in some circumstances, so it's
    %disabled for now. 
    
%     if ~isempty(f) && sum(f(:)) ~= 0
%         %feature rescaling
%         fmin = min(min(min(f)));
%         fmax = max(max(max(f)));
%         fprime = (f - fmin)./(fmax - fmin);
% 
%         %reduce precision and recast for memory and speed
%         ndigits = ceil(abs(log(std(fprime(:)))));
%         fround = round(fprime, ndigits);
%         fint = int16(fround*10^ndigits);
%     else
        fint = f;
%     end
    
end