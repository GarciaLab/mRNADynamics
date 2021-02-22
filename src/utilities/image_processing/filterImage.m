function [im, successFlag] = filterImage(im, filterType, sigmas, varargin)
%
%
%
% supportedFilters = {'Gaussian_blur', 'Identity', 'Anisotropic_diffusion', 'bilateral',...
%     'Edges', 'Sobel', 'Difference_of_Gaussian', 'Laplacian', 'Mean', 'Structure_smallest', 'Structure_largest',...
%     'Median',  'Maximum', 'Minimum', 'Std', 'Variance',...
%       'Hessian_smallest', 'Hessian_largest', 'Entropy', 'Range', 'imsegkmeans'};
%


zStep = 400; %nm/voxel
rad = 3; %rule of thumb is kernel size is 3x the Gaussian sigma
sigmaZ = 280 / zStep; %280 nm. don't remember why i picked this, but it works for my data -AR
filterSizeZ = ceil(sigmaZ*rad);

successFlag = true;
numType = 'double';
padding = 'symmetric';

for arg = 1:length(varargin)
    if strcmpi(varargin{arg}, 'filterSize')
        filterSizeXY = varargin{arg+1};
    elseif strcmpi(varargin{arg}, 'zStep')
        zStep  = varargin{arg+1};
    end
end

gpu = strcmpi(class(im), 'gpuArray');

[s1, s2] = ...
    getSigmas(sigmas);

if isnan(s1)
    successFlag = false;
%     warning('FilterImage: Filter not recognized. Returning original image.')
    return
end

if ~exist('filterSizeXY', 'var')
    filterSizeXY = round(rad*s1);
end

if ~mod(filterSizeXY,2)
    filterSizeXY = filterSizeXY + 1;
end

s1Mat = [s1 s1];
s2Mat = [s2 s2];

yDim=size(im, 1); xDim = size(im, 2);
dim = length(size(im));
if dim==3
    zDim = size(im, 3);
    s1Mat = [s1Mat, sigmaZ];
end

switch filterType
    
    case 'Identity'
        %return the image
    case 'Gaussian_blur'
        if dim == 2
            im = imgaussfilt(im,s1);
        elseif dim == 3
            d = og3(s1, sigmaZ);
            im = imfilter(im, d, 'same', padding);
        end
    case 'Anisotropic_diffusion'
        if dim == 2
            im = imdiffusefilt(im,'NumberOfIterations', 20);
            %this probably will not produce the same results as the fiji
            %implementation. 
        elseif dim == 3
           %not supported
           successFlag = false;
        end
    case 'bilateral'
        if dim == 2
            im = imbilatfilt(im, 'Padding','symmetric');
            %this probably will not produce the same results as the fiji
            %implementation. 
        elseif dim == 3
           %not supported
           successFlag = false;
        end
    case 'bright_spot_psf'
        if dim == 2
            im = filtBrightSpot;
        else
            successFlag = false;
        end
    case 'imsegkmeans'
         if dim==2
             k = 2;
             im = imcomplement(imsegkmeans(imgaussfilt(single(im),s1), k));
         else
             successFlag = false;
         end
    case 'Edges'
        im = canny(im, s1Mat);
    case 'Sobel'
        if dim == 2
            im = edge(imgaussfilt(im,s1), 'Sobel');
        elseif dim == 3
           %not supported
           successFlag = false;
        end    
    case 'Difference_of_Gaussian'
        if s2 < s1
            error('DoG filter requires sigma 1 < sigma 2')
        end
        %assumes sigma 2 > sigma 1
        if dim==2
            d = DoG(filterSizeXY, s1, s2);
            im = imfilter(im, d, 'same', padding);
        elseif dim == 3
            d = DoG3(s1, s2, sigmaZ);
            im = imfilter(im, d, 'same', padding);
        end
    case 'Laplacian'
        if dim==2
            h = fspecial('log', filterSizeXY,s1);
            im = imfilter(im, h, 'corr', 'symmetric', 'same');
        elseif dim==3
            h = fspecial3('log', [filterSizeXY, filterSizeXY, filterSizeZ], [s1, s1, sigmaZ]);
            im = -imfilter(im, h, 'corr', 'symmetric', 'same');
        end
    case 'Mean'
        if dim == 2
            h = fspecial('average', s1);
            im = imfilter(im, h,  'same', padding);
        elseif dim == 3
            h = fspecial3('average',filterSizeXY);
            im = imfilter(im, h, 'symmetric', 'same');
        end
    case {'Structure_smallest', 'Structure_largest'}
        if dim==2
            G=fspecial('gauss',[filterSizeXY, filterSizeXY], s1);
            [Gx,Gy] = gradient(G);
            %Compute Gaussian partial derivatives
            Dx = conv2(im, abs(Gx),'same');
            Dy = conv2(im, abs(Gy), 'same');
            %Smooth elements of the structure tensor
            S11 = conv2(Dx.^2,G,'same');
            S12 = conv2(Dx.*Dy,G,'same');
            S21 = S12;
            S22 = conv2(Dy.^2,G,'same');
            %Make eigenimages from the structure tensors
            for y = 1:yDim
                for x = 1:xDim
                    S = [S11(y, x), S12(y, x); S21(y, x), S22(y,x)];
                    if strcmpi(filterType, 'Structure_smallest')
                        im(y,x) = min(eig(S));
                    elseif strcmpi(filterType, 'Structure_largest')
                        im(y,x) = max(eig(S));
                    end
                end
            end
        elseif dim==3
            G = imgaussfilt3(im, [s1, s1, sigmaZ]);
            [Gx,Gy,Gz] = gradient(G);
            %Compute Gaussian partial derivatives
            Dx = imfilter(im, abs(Gx), 'corr', 'same', 'symmetric');
            Dy = imfilter(im, abs(Gy), 'corr', 'same', 'symmetric');
            Dz = imfilter(im, abs(Gz), 'corr', 'same', 'symmetric');
            
            %Smooth elements of the structure tensor
            S11 = imgaussfilt3(Dx.^2, [s1, s1, sigmaZ]);
            S12 = imgaussfilt3(Dx.*Dy, [s1, s1, sigmaZ]);
            S13 = imgaussfilt3(Dx.*Dz, [s1, s1, sigmaZ]);
            S22 =  imgaussfilt3(Dy.^2, [s1, s1, sigmaZ]);
            S23 = imgaussfilt3(Dy.*Dz, [s1, s1, sigmaZ]);
            S33 = imgaussfilt3(Dz.^2, [s1, s1, sigmaZ]);
            %Make eigenimages from the structure tensors
            im = zeros(yDim,xDim,zDim);
            parfor y = 1:yDim
                for x = 1:xDim
                    for z = 1:zDim
                        S = [S11(y,x,z), S12(y,x,z), S13(y,x,z);...
                            S12(y,x,z), S22(y,x,z), S23(y,x,z);...
                            S13(y,x,z), S23(y,x,z), S33(y,x,z)];
                        if strcmpi(filterType, 'Structure_smallest')
                            im(y,x,z) = min(eig(S));
                        elseif strcmpi(filterType, 'Structure_largest')
                            im(y,x,z) = max(eig(S));
                        end
                    end
                end
            end
        end
    case 'Median'
        im = gather(im);
        
        if dim==2
            %             f = imgaussfilt(im,s1);
            %             f = ordfilt2(f,ceil(filterSizeXY*filterSizeXY/2),ones(filterSizeXY,filterSizeXY));
            im = medfilt2(im, [filterSizeXY, filterSizeXY], padding);
        elseif dim==3
            im = medfilt3(im, [filterSizeXY, filterSizeXY, filterSizeZ], padding);
        end
        if gpu
            im = gpuArray(im);
        end
        
    case 'Maximum'
        if dim==2
            im = imgaussfilt(im,s1);
            se = strel('disk',ceil(filterSizeXY/2));
            im = imdilate(im,se);
            %                 f = ordfilt2(f,(filterSize*filterSize),ones(filterSize,filterSize));
            %                 f = imgaussfilt(f,s);
        elseif dim==3
            im = imgaussfilt3(im, [s1, s1, sigmaZ]);
            se = strel('cuboid',[ceil(filterSizeXY/2), ceil(filterSizeXY/2), ceil(filterSizeZ/2)]);
            im = gather(im);
            im = imdilate(im, se);
            if gpu
                im = gpuArray(im);
            end
        end
        
    case 'Minimum'
        if dim==2
            im = imgaussfilt(im,s1);
            se = strel('disk',ceil(filterSizeXY/2));
            im = imerode(im,se);
        elseif dim==3
            im = imgaussfilt3(im, [s1, s1, sigmaZ]);
            se = strel('cuboid',[ceil(filterSizeXY/2), ceil(filterSizeXY/2), ceil(filterSizeZ/2)]);
            im = gather(im);
            im = imerode(im, se);
            if gpu
                im = gpuArray(im);
            end
        end
    case {'Entropy'}
        if dim==2
            im = entropyfilt(imgaussfilt(im,s1));
        elseif dim==3
           successFlag = false;
        end
      case {'Range'}
        if dim==2
            im = rangefilt(imgaussfilt(im,s1));
        elseif dim==3
            successFlag = false;
        end
    case {'Std', 'Variance'}
        if dim==2
            im = imgaussfilt(im,s1);
            opts = {};
            %             opts  = [opts, 'gpuArray'];
            im = stdfilt(im,ones(filterSizeXY, filterSizeXY), opts{:});
        elseif dim==3
            %this blurs with sigma and sigmaZ, then stdfilts with sizes
            %dictated by sigma and sigma z.
            im = gather(im);
            im = stdfilt(imgaussfilt3(im, [s1, s1, sigmaZ]), ones(filterSizeXY, filterSizeXY, ceil(3*sigmaZ)));
            if gpu
                im = gpuArray(im);
            end
        end
        
        if strcmpi(filterType, 'Variance')
            im = im.^2;
        end
        
    case {'Hessian_smallest', 'Hessian_largest'}
        if dim==2
            G = fspecial('gauss',[filterSizeXY, filterSizeXY], s1);
            [Gx,Gy] =gradient(G);
            [Gxx, Gxy] = gradient(Gx);
            [Gyy, ~] = gradient(Gy);
            %Compute elements of the Hessian matrix
            H11 = conv2(im,abs(Gxx),'same');
            H12 = conv2(im,abs(Gxy),'same');
            H21 = H12;
            H22 = conv2(im,abs(Gyy),'same');
            %Make eigenimages from the Hessian
            for y = 1:yDim
                for x = 1:xDim
                    H = [H11(y, x), H12(y, x); H21(y, x), H22(y,x)];
                    if strcmpi(filterType, 'Hessian_smallest')
                        im(y,x) = -min(eig(H));
                    elseif strcmpi(filterType, 'Hessian_largest')
                        im(y,x) = max(eig(H));
                    end
                end
            end
        elseif dim ==3
            G1 = og3(s1, sigmaZ);
            [Gx,Gy,Gz] = gradient(G1);
            [Gxx, Gxy, Gxz] = gradient(Gx);
            [Gyy, ~, Gyz] = gradient(Gy);
            [~, ~, Gzz] = gradient(Gz);
            %Compute elements of the Hessian matrix
            H11 = imfilter(im, abs(Gxx), 'corr', 'same', 'symmetric');
            H12 = imfilter(im, abs(Gxy), 'corr', 'same', 'symmetric');
            H13 = imfilter(im, abs(Gxz), 'corr', 'same', 'symmetric');
            H22 = imfilter(im, abs(Gyy), 'corr', 'same', 'symmetric');
            H23 = imfilter(im, abs(Gyz), 'corr', 'same', 'symmetric');
            H33 = imfilter(im, abs(Gzz), 'corr', 'same', 'symmetric');
            %Make eigenimages from the Hessian
            im = zeros(yDim,xDim,zDim);
            parfor y = 1:yDim
                for x = 1:xDim
                    for z = 1:zDim
                        H = [H11(y, x, z), H12(y, x, z), H13(y, x, z);...
                            H12(y, x, z), H22(y,x, z), H23(y,x,z);...
                            H13(y,x,z), H23(y,x,z), H33(y,x,z)];
                        if strcmpi(filterType, 'Hessian_smallest')
                            im(y,x,z) = -min(eig(H));
                        elseif strcmpi(filterType, 'Hessian_largest')
                            im(y,x,z) = max(eig(H));
                        end
                    end
                end
            end
        end
        
    otherwise
%     warning('FilterImage: Filter not recognized. Returning original image.')
        successFlag = false;
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

    function og = og3(s1, sigmaZ)
        filterSizeXY = round(s1*3);
        filterSizeZ = round(sigmaZ*4*3);
        og = fspecial3('gaussian',[filterSizeXY, filterSizeXY, filterSizeZ], [s1,s1,ceil(sigmaZ*4)]);
    end

    function [s1, s2] = ...
            getSigmas(sigmas)
        s1 = NaN;
        s2 = NaN;
        if ~iscell(sigmas)
            sigmas = {sigmas};
        end
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
    end


end