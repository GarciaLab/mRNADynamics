function f = filterImage(im, filterType, sigmas, varargin)

persistent dh; persistent dh3;
% dh = memoize(@DoG); dh3 = memoize(@DoG3);
dh3 = @DoG3;
persistent g3;
% g3 = memoize(@gauss3D);
persistent Gx; persistent Gy; persistent Gz;
% grad = memoize(@gradient);
persistent Gxx; persistent Gxy; persistent Gxz; persistent Gyy; persistent Gyz; persistent Gzz;

zStep = 400; %nm. default.

%initialize
filterSizeZ = NaN; s1 = NaN; s2 = NaN; sigmaZ = NaN; filterSizeXY = NaN;

gpu = strcmpi(class(im), 'gpuArray');

padding = 'symmetric';
for args = 1:length(varargin)
    if strcmpi(varargin{args}, 'filterSize')
        filterSizeXY = round(varargin{args+1});
    elseif strcmpi(varargin{args}, 'zStep')
        zStep = varargin{args+1};
    elseif strcmpi(varargin{args}, 'padding')
        padding = varargin{args+1};
    end
end

getSigmasAndFilterSizes;

dim = length(size(im));
f = [];

switch filterType
    case 'Identity'
        f=im;
    case 'Gaussian_blur'
        if dim == 2
            f = imgaussfilt(im,s1);
        elseif dim == 3
            f = imgaussfilt3(im, [s1, s1, sigmaZ]);
        end
    case 'bright_spot_psf'
        if dim == 2
            %21 x 21
            spotfilt = [0,0.150000000000000,0.0500000000000000,0,0.0500000000000000,0,0.100000000000000,0.100000000000000,0,0.100000000000000,0.0500000000000000,0.0500000000000000,0,0.100000000000000,0,0,0,0,0.0500000000000000,0.100000000000000,0.100000000000000;0.0500000000000000,0.0500000000000000,0,0,0.100000000000000,0,0.0500000000000000,0,0,0,0.0500000000000000,0.0500000000000000,0,0,0.100000000000000,0.150000000000000,0,0,0.200000000000000,0.100000000000000,0.100000000000000;0.0500000000000000,0.0500000000000000,0,0.0500000000000000,0.0500000000000000,0,0,0,0.100000000000000,0.0500000000000000,0,0.0500000000000000,0,0.0500000000000000,0,0,0.200000000000000,0,0.100000000000000,0.100000000000000,0;0.0500000000000000,0.0500000000000000,0,0.0500000000000000,0,0.100000000000000,0,0.150000000000000,0,0.0500000000000000,0.0500000000000000,0,0.0500000000000000,0.0500000000000000,0,0.0500000000000000,0,0.100000000000000,0,0.100000000000000,0.0500000000000000;0,0,0.100000000000000,0.150000000000000,0.0500000000000000,0.100000000000000,0,0.150000000000000,0.0500000000000000,0.100000000000000,0.100000000000000,0.0500000000000000,0.100000000000000,0.150000000000000,0,0.150000000000000,0.0500000000000000,0,0,0.0500000000000000,0;0,0.100000000000000,0.200000000000000,0.100000000000000,0,0,0.100000000000000,0.100000000000000,0.200000000000000,0.300000000000000,0.150000000000000,0,0.200000000000000,0,0.150000000000000,0,0.0500000000000000,0,0.0500000000000000,0,0.100000000000000;0,0.0500000000000000,0,0,0,0.0500000000000000,0,0.150000000000000,0.100000000000000,0.150000000000000,0.150000000000000,0.200000000000000,0.250000000000000,0.200000000000000,0.250000000000000,0.0500000000000000,0.150000000000000,0,0.0500000000000000,0.0500000000000000,0.0500000000000000;0.0500000000000000,0.0500000000000000,0.100000000000000,0.0500000000000000,0.0500000000000000,0.0500000000000000,0.200000000000000,0.0500000000000000,0.200000000000000,0.250000000000000,0.300000000000000,0.300000000000000,0.0500000000000000,0.350000000000000,0.200000000000000,0.150000000000000,0.150000000000000,0.0500000000000000,0.0500000000000000,0.0500000000000000,0.0500000000000000;0.0500000000000000,0.100000000000000,0.0500000000000000,0.0500000000000000,0.150000000000000,0.150000000000000,0.550000000000000,0.200000000000000,0.250000000000000,0.650000000000000,0.550000000000000,0.500000000000000,0.550000000000000,0.250000000000000,0.350000000000000,0.0500000000000000,0.100000000000000,0.200000000000000,0.150000000000000,0.150000000000000,0;0.0500000000000000,0,0.100000000000000,0.200000000000000,0.200000000000000,0.250000000000000,0.0500000000000000,0.150000000000000,0.300000000000000,1,0.500000000000000,0.650000000000000,0.750000000000000,0.450000000000000,0.350000000000000,0.150000000000000,0.150000000000000,0.100000000000000,0.100000000000000,0,0.0500000000000000;0,0,0.0500000000000000,0,0.250000000000000,0.100000000000000,0.400000000000000,0.200000000000000,0.600000000000000,0.500000000000000,0.750000000000000,0.700000000000000,0.550000000000000,0.400000000000000,0.250000000000000,0.0500000000000000,0,0.150000000000000,0.200000000000000,0.0500000000000000,0.0500000000000000;0.0500000000000000,0.0500000000000000,0.0500000000000000,0.100000000000000,0.150000000000000,0.150000000000000,0.500000000000000,0.500000000000000,0.400000000000000,0.550000000000000,0.600000000000000,0.600000000000000,0.650000000000000,0.500000000000000,0.150000000000000,0.150000000000000,0.150000000000000,0.0500000000000000,0.100000000000000,0,0.100000000000000;0,0.200000000000000,0,0.0500000000000000,0.0500000000000000,0.250000000000000,0.150000000000000,0.300000000000000,0.450000000000000,0.400000000000000,0.600000000000000,0.550000000000000,0.350000000000000,0.300000000000000,0.200000000000000,0.100000000000000,0.0500000000000000,0.0500000000000000,0.100000000000000,0.0500000000000000,0;0.100000000000000,0.100000000000000,0.0500000000000000,0.100000000000000,0.100000000000000,0.300000000000000,0.0500000000000000,0.350000000000000,0.250000000000000,0.750000000000000,0.450000000000000,0.500000000000000,0.550000000000000,0.250000000000000,0.150000000000000,0.0500000000000000,0.0500000000000000,0.200000000000000,0,0.100000000000000,0.0500000000000000;0.0500000000000000,0.0500000000000000,0.100000000000000,0.0500000000000000,0.250000000000000,0.150000000000000,0.250000000000000,0.150000000000000,0.200000000000000,0.350000000000000,0.250000000000000,0.300000000000000,0.150000000000000,0.100000000000000,0,0.0500000000000000,0.100000000000000,0.150000000000000,0.100000000000000,0,0;0,0.100000000000000,0.0500000000000000,0.0500000000000000,0.0500000000000000,0,0.0500000000000000,0.200000000000000,0.200000000000000,0.100000000000000,0.300000000000000,0.100000000000000,0.200000000000000,0,0.150000000000000,0.0500000000000000,0.100000000000000,0,0.100000000000000,0.150000000000000,0;0,0.100000000000000,0.150000000000000,0,0.0500000000000000,0.0500000000000000,0.0500000000000000,0.250000000000000,0.0500000000000000,0.100000000000000,0.100000000000000,0.0500000000000000,0.0500000000000000,0,0,0,0.0500000000000000,0.0500000000000000,0.0500000000000000,0.100000000000000,0.0500000000000000;0,0,0.0500000000000000,0.0500000000000000,0.100000000000000,0.100000000000000,0.100000000000000,0.150000000000000,0.0500000000000000,0,0.100000000000000,0.0500000000000000,0.0500000000000000,0.150000000000000,0.0500000000000000,0,0.150000000000000,0.100000000000000,0,0.200000000000000,0.150000000000000;0.0500000000000000,0.150000000000000,0.100000000000000,0.0500000000000000,0.100000000000000,0.100000000000000,0,0.0500000000000000,0.150000000000000,0.100000000000000,0.100000000000000,0.0500000000000000,0,0.100000000000000,0.0500000000000000,0.0500000000000000,0.0500000000000000,0,0.0500000000000000,0,0.150000000000000;0,0,0.0500000000000000,0,0.0500000000000000,0.150000000000000,0.100000000000000,0.0500000000000000,0.100000000000000,0.0500000000000000,0.150000000000000,0.0500000000000000,0,0.150000000000000,0.0500000000000000,0.100000000000000,0,0.0500000000000000,0.0500000000000000,0.0500000000000000,0.0500000000000000;0.0500000000000000,0.0500000000000000,0.0500000000000000,0,0,0,0.0500000000000000,0.150000000000000,0.0500000000000000,0,0.100000000000000,0.100000000000000,0.0500000000000000,0.100000000000000,0,0,0.0500000000000000,0.0500000000000000,0.100000000000000,0.100000000000000,0.0500000000000000];
            f = conv2(im, spotfilt, 'same');
        else
            f = im;
        end
    case 'dim_spot_psf'
        %             f = spotfilt(im,s1, dim, 'small');
        f = im;
    case 'Edges'
        if dim == 2
            f = edge(im, 'canny', s1);
        elseif dim == 3
            f = canny(im, [s1, s1, sigmaZ]);
        end
    case 'Difference_of_Gaussian'
        if s2 < s1
            error('DoG filter requires sigma 1 < sigma 2')
        end
        %assumes sigma 2 > sigma 1
        if dim==2
            d = DoG(filterSizeXY, s1, s2);
            %             try
            %                 gp = gpuDevice;
            %                 if gp.AvailableMemory < 1E9
            %                     gp = gpuDevice(1);
            %                 end
            %                 gim = gpuArray(im);
            %                 gd = gpuArray(d);
            %                 gpuf = conv2(gim, gd,'same');
            %                 f = gather(gpuf); clear gpuf; clear gim; clear gd;
            %             catch
            f = imfilter(im, d, 'same', padding);
            %         end
        elseif dim == 3
            if gpu
                d = gpuArray(dh3(s1, s2, sigmaZ));
            else
                disp('Defaulting to CPU.');
                d = dh3(s1, s2, sigmaZ);
            end
            f = imfilter(im, d, 'same', padding);
        end
    case 'Laplacian'
        if dim==2
            h = fspecial('log', filterSizeXY,s1);
            f = imfilter(im, h, 'corr', 'symmetric', 'same');
        elseif dim==3
            %             G = g3(s1, 'sigmaZ', sigmaZ);
            %             [Gx,Gy] = grad(G);
            %             [Gxx,~] = grad(Gx);
            %             [~,Gyy] = grad(Gy);
            %             LoG = Gxx + Gyy;
            %             f = -imfilter(im, LoG, 'corr', 'symmetric', 'same');
            %                 f = fastConv3D(LoG, im);
            h = fspecial3('log', [filterSizeXY, filterSizeXY, filterSizeZ], [s1, s1, sigmaZ]);
            f = -imfilter(im, h, 'corr', 'symmetric', 'same');
        end
    case 'Mean'
        if dim == 2
            h = fspecial('average', s1);
            f = imfilter(im, h);
        elseif dim == 3
            h = fspecial3('average',filterSizeXY);
            f = imfilter(im, h, 'symmetric', 'same');
        end
    case {'Structure_smallest', 'Structure_largest'}
        if dim==2
            G=fspecial('gauss',[filterSizeXY, filterSizeXY], s1);
            [Gx,Gy] = grad(G);
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
                    if strcmpi(filterType, 'Structure_smallest')
                        f(p,q) = min(eig(S));
                    elseif strcmpi(filterType, 'Structure_largest')
                        f(p,q) = max(eig(S));
                    end
                end
            end
        elseif dim==3
            G = imgaussfilt3(im, [s1, s1, sigmaZ]);
            [Gx,Gy,Gz] = grad(G);
            %Compute Gaussian partial derivatives
            Dx = imfilter(im, Gx, 'corr', 'same', 'symmetric');
            Dy = imfilter(im, Gy, 'corr', 'same', 'symmetric');
            Dz = imfilter(im, Gz, 'corr', 'same', 'symmetric');
            
            %Smooth elements of the structure tensor
            S11 = imgaussfilt3(Dx.^2, [s1, s1, sigmaZ]);
            S12 = imgaussfilt3(Dx.*Dy, [s1, s1, sigmaZ]);
            S13 = imgaussfilt3(Dx.*Dz, [s1, s1, sigmaZ]);
            S22 =  imgaussfilt3(Dy.^2, [s1, s1, sigmaZ]);
            S23 = imgaussfilt3(Dy.*Dz, [s1, s1, sigmaZ]);
            S33 = imgaussfilt3(Dz.^2, [s1, s1, sigmaZ]);
            %Make eigenimages from the structure tensors
            l = size(im, 1); m = size(im, 2); n = size(im, 3);
            f = zeros(l,m,n);
            parfor p = 1:l
                for q = 1:m
                    for r = 1:n
                        S = [S11(p,q,r), S12(p,q,r), S13(p,q,r);...
                            S12(p,q,r), S22(p,q,r), S23(p,q,r);...
                            S13(p,q,r), S23(p,q,r), S33(p,q,r)];
                        if strcmpi(filterType, 'Structure_smallest')
                            f(p,q,r) = min(eig(S));
                        elseif strcmpi(filterType, 'Structure_largest')
                            f(p,q,r) = max(eig(S));
                        end
                    end
                end
            end
        end
    case 'Median'
        if dim==2
            %             f = imgaussfilt(im,s1);
            %             f = ordfilt2(f,ceil(filterSizeXY*filterSizeXY/2),ones(filterSizeXY,filterSizeXY));
            f = medfilt2(im, [filterSizeXY, filterSizeXY]);
        elseif dim==3
            f = medfilt3(im, [filterSizeXY, filterSizeXY, filterSizeZ]);
            f = im;
        end
    case 'Maximum'
        if dim==2
            f = imgaussfilt(im,s1);
            se = strel('disk',ceil(filterSizeXY/2));
            f = imdilate(f,se);
            %                 f = ordfilt2(f,(filterSize*filterSize),ones(filterSize,filterSize));
            %                 f = imgaussfilt(f,s);
        elseif dim==3
            f = imgaussfilt3(im, [s1, s1, sigmaZ]);
            se = strel('cuboid',[ceil(filterSizeXY/2), ceil(filterSizeXY/2), ceil(filterSizeZ/2)]);
            f = imdilate(f, se);
        end
        
    case 'Minimum'
        if dim==2
            f = imgaussfilt(im,s1);
            se = strel('disk',ceil(filterSizeXY/2));
            f = imerode(f,se);
        elseif dim==3
            f = imgaussfilt3(im, [s1, s1, sigmaZ]);
            se = strel('cuboid',[ceil(filterSizeXY/2), ceil(filterSizeXY/2), ceil(filterSizeZ/2)]);
            f = imerode(f, se);
        end
    case {'Std', 'Variance'}
        if dim==2
            f = imgaussfilt(im,s1);
            f = stdfilt(f);
        elseif dim==3
            %this blurs with sigma and sigmaZ, then stdfilts with sizes
            %dictated by sigma and sigma z.
            f = stdfilt(imgaussfilt3(im, [s1, s1, sigmaZ]), ones(filterSizeXY, filterSizeXY, ceil(3*sigmaZ)));
        end
        
        if strcmpi(filterType, 'Variance')
            f = f.^2;
        end
    case {'Hessian_smallest', 'Hessian_largest'}
        if dim==2
            G = fspecial('gauss',[filterSizeXY, filterSizeXY], s1);
            [Gx,Gy] = grad(G);
            [Gxx, Gxy] = grad(Gx);
            [Gyy, ~] = grad(Gy);
            %Compute elements of the Hessian matrix
            H11 = conv2(im,Gxx,'same');
            H12 = conv2(im,Gxy,'same');
            H21 = H12;
            H22 = conv2(im,Gyy,'same');
            %Make eigenimages from the Hessian
            for p = 1:size(im, 1)
                for q = 1:size(im, 2)
                    H = [H11(p, q), H12(p, q); H21(p, q), H22(p,q)];
                    if strcmpi(filterType, 'Hessian_smallest')
                        f(p,q) = -min(eig(H));
                    elseif strcmpi(filterType, 'Hessian_largest')
                        f(p,q) = max(eig(H));
                    end
                end
            end
        elseif dim ==3
            G1 = g3(s1, 'sigmaZ', sigmaZ);
            [Gx,Gy,Gz] = grad(G1);
            [Gxx, Gxy, Gxz] = grad(Gx);
            [Gyy, ~, Gyz] = grad(Gy);
            [~, ~, Gzz] = grad(Gz);
            %Compute elements of the Hessian matrix
            H11 = imfilter(im, Gxx, 'corr', 'same', 'symmetric');
            H12 = imfilter(im, Gxy, 'corr', 'same', 'symmetric');
            H13 = imfilter(im, Gxz, 'corr', 'same', 'symmetric');
            H22 = imfilter(im, Gyy, 'corr', 'same', 'symmetric');
            H23 = imfilter(im, Gyz, 'corr', 'same', 'symmetric');
            H33 = imfilter(im, Gzz, 'corr', 'same', 'symmetric');
            %Make eigenimages from the Hessian
            l = size(im, 1); m = size(im, 2); n = size(im, 3);
            f = zeros(l,m,n);
            for p = 1:l
                for q = 1:m
                    for r = 1:n
                        H = [H11(p, q, r), H12(p, q, r), H13(p, q, r);...
                            H12(p, q, r), H22(p,q, r), H23(p,q,r);...
                            H13(p,q,r), H23(p,q,r), H33(p,q,r)];
                        if strcmpi(filterType, 'Hessian_smallest')
                            f(p,q,r) = -min(eig(H));
                        elseif strcmpi(filterType, 'Hessian_largest')
                            f(p,q,r) = max(eig(H));
                        end
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
% fint = f;
%     end

    function dog = DoG(filterSize, s1, s2)
        dog = fspecial('gaussian',filterSize, s1) - fspecial('gaussian',filterSize, s2);
    end

    function dog = DoG3(s1, s2, sigmaZ)
        filterSizeXY = round(s2*3);
        filterSizeZ = round(sigmaZ*4*3);
        dog = fspecial3('gaussian',[filterSizeXY, filterSizeXY, filterSizeZ], [s1,s1,ceil(sigmaZ*4)]) - fspecial3('gaussian',[filterSizeXY, filterSizeXY, filterSizeZ], [s2,s2,ceil(sigmaZ*4)]);
    end

    function getSigmasAndFilterSizes
        rad = 3; %rule of thumb is kernel size is 3x the Gaussian sigma
        zStep = 400; %nm. default.
        %convert string sigmas to doubles
        if ~isempty(sigmas)
            if ischar(sigmas{1})
                for ii=1:length(sigmas)
                    sigmas{ii}=str2double(sigmas{ii});
                end
            end
        end
        
        len = length(sigmas);
        switch len
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
            filterSizeXY = round(rad*s1);
        end
        
        sigmaZ = 280 / zStep;
        filterSizeZ = ceil(sigmaZ*3);
        
        if ~mod(filterSizeXY,2)
            filterSizeXY = filterSizeXY + 1;
        end
    end


end