function [ x,y, varargout ] = interpolatedShift(FrameInfo, img1, img2, maxRadius, localMaxRadius, maxShift, precisionFactor, XY)
%INTERPOLATEDSHIFT Summary of this function goes here
%   Detailed explanation goes here


if ~exist('precisionFactor','var') || isempty(precisionFactor)
    precisionFactor = 1;
end

space_resolution = getDefaultParameters(FrameInfo,'space resolution');
tileSizeX = round(1*128*0.22/space_resolution/precisionFactor);
tileSizeY = round(1*128*0.22/space_resolution/precisionFactor);
subTileSize = round(1*64*0.22/space_resolution/precisionFactor);



if ~exist('maxRadius','var') || isempty(maxRadius)
    mask = ones(tileSizeX,tileSizeY);
    maxRadius = min(tileSizeX,tileSizeY)/8;
    localMaxRadius = maxRadius;
    maxShift = maxRadius*10;%7.5;
else
    [X,Y] = meshgrid(1:tileSizeX,1:tileSizeY);
    mu = 0.5*[tileSizeX tileSizeY];
    Sigma = eye(2)*maxRadius;
    mask = mvnpdf([X(:) Y(:)],mu,Sigma);
    mask = reshape(mask,tileSizeX,tileSizeY);
end

if ~exist('localMaxRadius','var') || isempty(localMaxRadius)
    localMaxRadius = 0.3*maxRadius;
    maxShift = maxRadius*10;%7.5;
end

if ~exist('maxShift','var') || isempty(maxShift)
    
    maxShift = maxRadius*10;%7.5;
end



localMaxMask = fspecial('disk',localMaxRadius);
localMaxMask = im2bw(mat2gray(localMaxMask),graythresh(mat2gray(localMaxMask)));
localMaxMask(round(length(localMaxMask)*0.5),round(length(localMaxMask)*0.5)) = 0;

%nTilesX = size(img1,1)/tileSizeX;
%nTilesY = size(img1,2)/tileSizeY;
nTilesX = round(size(img1,1)/subTileSize);
nTilesY = round(size(img1,2)/subTileSize);

img1 = medfilt2(img1);
img2 = medfilt2(img2);

for j = 1:nTilesY
    
    for jj = 1:nTilesX
        centerX = (jj-0.5)*subTileSize;
        centerY = (j-0.5)*subTileSize;
        spanX = max(1,round(centerX-0.5*tileSizeX)):min(size(img1,1),round(centerX+0.5*tileSizeX));
        spanY = max(1,round(centerY-0.5*tileSizeY)):min(size(img1,2),round(centerY+0.5*tileSizeY));
        fil = img1;%(spanX,spanY);%(img1((1+(jj-1)*tileSizeX):(tileSizeX+(jj-1)*tileSizeX),(1+(j-1)*tileSizeY):(tileSizeY+(j-1)*tileSizeY));
        
        [X,Y] = meshgrid(1:size(fil,1),1:size(fil,2));
        mu = [centerX centerY];
        Sigma = eye(2)*maxRadius*7;
        mask_tmp = mvnpdf([X(:) Y(:)],mu,Sigma);
        mask_tmp = reshape(mask_tmp,size(fil,2),size(fil,1))';
        
        fil = fil-mean(mean(fil(spanX,spanY)));
        
        cimg2 = img2-mean(mean(img2(spanX,spanY)));%img2-mean(mean(img2(((1+(jj-1)*tileSizeX):(tileSizeX+(jj-1)*tileSizeX)),(1+(j-1)*tileSizeY):(tileSizeY+(j-1)*tileSizeY))));
        f = fourierFilterInAGaussianWindow(fil.*mask_tmp,cimg2);
        minX = max(1,floor(size(f,1)*0.5-.5*tileSizeX));
        maxX = min(ceil(size(f,1)*0.5+.5*tileSizeX),size(f,1));
        minY = max(1,floor(size(f,2)*0.5-.5*tileSizeY));
        maxY = min(ceil(size(f,2)*0.5+.5*tileSizeY),size(f,2));
        f = f(minX:maxX,minY:maxY);%f(spanX,spanY);%f((1+(jj-1)*tileSizeX):(tileSizeX+(jj-1)*tileSizeX),(1+(j-1)*tileSizeY):(tileSizeY+(j-1)*tileSizeY));
        
        useGPU = 0;
        try
            gpuDevice;
            useGPU = 1;
        catch
            %do nothing;
        end
        
        if useGPU
            maxima = (f >= gather(imdilate(gpuArray(uint8(f)),localMaxMask)));
        else
            maxima = f >= imdilate(f,localMaxMask);
        end
        
        if any(size(f) ~= [tileSizeX,tileSizeY])
            [X,Y] = meshgrid(1:size(f,1),1:size(f,2));
            mu = 0.5*size(f)+0.5;
            Sigma = eye(2)*maxShift;
            mask_tmp = mvnpdf([X(:) Y(:)],mu,Sigma);
            mask_tmp = reshape(mask_tmp,size(f,2),size(f,1))';
            IND = find(f(maxima).*mask_tmp(maxima) == max(f(maxima).*mask_tmp(maxima)));
        else
            IND = find(f(maxima).*mask(maxima) == max(f(maxima).*mask(maxima)));
        end
        indMaxima = find(maxima);
        
        [tx,ty] = ind2sub(size(f),indMaxima(IND));
        
        x(j,jj) = mean((0.5*size(f,1)+0.5)-tx);
        
        y(j,jj) = mean((0.5*size(f,2)+0.5)-ty);
    end
end

%% Interpolate Values
if nargout > 2
    
    if exist('XY','var') && ~isempty(XY)
        
        method = 2;
        
        switch method
            case 1
                %% Using griddata : cubic spline interpolation
                %[xi,yi] = meshgrid(1:size(img2,1),1:size(img2,2));
                xi = XY(:,1); yi = XY(:,2);
                halfTileX = 0.5*subTileSize;
                halfTileY = 0.5*subTileSize;
                px = repmat( (halfTileX:subTileSize:(nTilesX*subTileSize-halfTileX))', 1, nTilesY );
                py = repmat( halfTileY:subTileSize:(nTilesY*subTileSize-halfTileY), nTilesX, 1 );% [0.5*tileSizeX*ones(nTilesY,1) (size(img2,1)-0.5*tileSizeX)*ones(nTilesY,1)];
                %py = [0.5*tileSizeY*ones(1,nTilesX); 192.5*ones(1,nTilesX); 320.5*ones(1,nTilesX); 448.5*ones(1,nTilesX)];
                Px = [0.5*ones(2,1) (size(img2,1)+0.5)*ones(2,1)];
                Py = [0.5*ones(1,2); (size(img2,2)+0.5)*ones(1,2)];
                vx = x([1,nTilesY,nTilesX*nTilesY-nTilesY+1,nTilesX*nTilesY]);
                vy = y([1,nTilesY,nTilesX*nTilesY-nTilesY+1,nTilesX*nTilesY]);
                
                warning('off','MATLAB:griddata:DuplicateDataPoints');
                gx = flipud(griddata([px(:);Px(:)],[py(:);Py(:)],[x(:); vx(:)],xi,yi,'cubic'));
                %gx = flipud(ksrmv([px(:) py(:);Px(:) Py(:)],[y(:);vy(:)],[],[xi(:),yi(:)]));
                gy = flipud(griddata([px(:);Px(:)],[py(:);Py(:)],[y(:); vy(:)],xi,yi,'cubic'));
                %gy = flipud(ksrmv([px(:) py(:);Px(:) Py(:)],[x(:);vx(:)],[],[xi(:),yi(:)]));
                %     gx = flipud(gaussianKernelRegression([px(:);Px(:)],[py(:);Py(:)],[y(:); vy(:)],xi,yi));
                %     gx = reshape(gx,size(img2,1),size(img2,2));
                %     gy = flipud(gaussianKernelRegression([px(:);Px(:)],[py(:);Py(:)],[x(:); vx(:)],xi,yi));
                %     gy = reshape(gy,size(img2,1),size(img2,2));
                warning('on','MATLAB:griddata:DuplicateDataPoints');
                varargout{1} = XY+[gx,gy];
                
            case 2
                
                %% using tpaps : thin plate spline interpolation
                
                halfTileX = 0.5*subTileSize;
                halfTileY = 0.5*subTileSize;
                px = repmat( (halfTileX:subTileSize:(nTilesX*subTileSize-halfTileX))', 1, nTilesY );
                py = repmat( halfTileY:subTileSize:(nTilesY*subTileSize-halfTileY), nTilesX, 1 );% [0.5*tileSizeX*ones(nTilesY,1) (size(img2,1)-0.5*tileSizeX)*ones(nTilesY,1)];
                Px = [0.5*ones(2,1) (size(img2,1)+0.5)*ones(2,1)];
                Py = [0.5*ones(1,2); (size(img2,2)+0.5)*ones(1,2)];
                vx = x([1,nTilesY,nTilesX*nTilesY-nTilesY+1,nTilesX*nTilesY]);
                vy = y([1,nTilesY,nTilesX*nTilesY-nTilesY+1,nTilesX*nTilesY]);
                
                X = [px(:);Px(:)];
                Y = [py(:);Py(:)];
                xy = [X';Y'];
                x = x'; y = y';
                z = [x(:)' vx;y(:)' vy];
                xyi = XY';
                
                st = tpaps(xy,z,0.6);
                interpolatedXY = fnval(st,xyi)';
                varargout{1} = interpolatedXY+XY;
                
        end
    else
        
        method = 2;
        
        switch method
            case 1
                %% Using griddata : cubic spline interpolation
                [xi,yi] = meshgrid(1:size(img2,1),1:size(img2,2));
                halfTileX = 0.5*subTileSize;
                halfTileY = 0.5*subTileSize;
                px = repmat( (halfTileX:subTileSize:(nTilesX*subTileSize-halfTileX))', 1, nTilesY );
                py = repmat( halfTileY:subTileSize:(nTilesY*subTileSize-halfTileY), nTilesX, 1 );% [0.5*tileSizeX*ones(nTilesY,1) (size(img2,1)-0.5*tileSizeX)*ones(nTilesY,1)];
                %py = [0.5*tileSizeY*ones(1,nTilesX); 192.5*ones(1,nTilesX); 320.5*ones(1,nTilesX); 448.5*ones(1,nTilesX)];
                Px = [0.5*ones(2,1) (size(img2,1)+0.5)*ones(2,1)];
                Py = [0.5*ones(1,2); (size(img2,2)+0.5)*ones(1,2)];
                vx = x([1,nTilesY,nTilesX*nTilesY-nTilesY+1,nTilesX*nTilesY]);
                vy = y([1,nTilesY,nTilesX*nTilesY-nTilesY+1,nTilesX*nTilesY]);
                
                warning('off','MATLAB:griddata:DuplicateDataPoints');
                gx = flipud(griddata([px(:);Px(:)],[py(:);Py(:)],[x(:); vx(:)],xi,yi,'cubic'));
                %gx = flipud(ksrmv([px(:) py(:);Px(:) Py(:)],[y(:);vy(:)],[],[xi(:),yi(:)]));
                gy = flipud(griddata([px(:);Px(:)],[py(:);Py(:)],[y(:); vy(:)],xi,yi,'cubic'));
                %gy = flipud(ksrmv([px(:) py(:);Px(:) Py(:)],[x(:);vx(:)],[],[xi(:),yi(:)]));
                %     gx = flipud(gaussianKernelRegression([px(:);Px(:)],[py(:);Py(:)],[y(:); vy(:)],xi,yi));
                %     gx = reshape(gx,size(img2,1),size(img2,2));
                %     gy = flipud(gaussianKernelRegression([px(:);Px(:)],[py(:);Py(:)],[x(:); vx(:)],xi,yi));
                %     gy = reshape(gy,size(img2,1),size(img2,2));
                warning('on','MATLAB:griddata:DuplicateDataPoints');
                varargout{1} = gx;
                varargout{2} = gy;
                
            case 2
                
                %% using tpaps : thin plate spline interpolation
                [xi,yi] = meshgrid(1:size(img2,1),1:size(img2,2));
                halfTileX = 0.5*subTileSize;
                halfTileY = 0.5*subTileSize;
                px = repmat( (halfTileX:subTileSize:(nTilesX*subTileSize-halfTileX))', 1, nTilesY );
                py = repmat( halfTileY:subTileSize:(nTilesY*subTileSize-halfTileY), nTilesX, 1 );% [0.5*tileSizeX*ones(nTilesY,1) (size(img2,1)-0.5*tileSizeX)*ones(nTilesY,1)];
                Px = [0.5*ones(2,1) (size(img2,1)+0.5)*ones(2,1)];
                Py = [0.5*ones(1,2); (size(img2,2)+0.5)*ones(1,2)];
                vx = x([1,nTilesY,nTilesX*nTilesY-nTilesY+1,nTilesX*nTilesY]);
                vy = y([1,nTilesY,nTilesX*nTilesY-nTilesY+1,nTilesX*nTilesY]);
                
                X = [px(:);Px(:)];
                Y = [py(:);Py(:)];
                xy = [X';Y'];
                x = x'; y = y';
                z = [x(:)' vx;y(:)' vy];
                xyi = [xi(:)';yi(:)'];
                
                st = tpaps(xy,z,0.9);
                vals = fnval(st,xyi);
                gx = vals(1,:);
                gx = reshape(gx,fliplr(size(img1)))';
                gy = vals(2,:);
                gy = reshape(gy,fliplr(size(img1)))';
                varargout{1} = gx;
                varargout{2} = gy;
                
        end
    end
end

end

