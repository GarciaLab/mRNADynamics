function [cMask, ellipseFrameWithEdges] = fitCirclesToNuclei(mask, varargin)

displayFigures = false;
doEllipse = true;
image = [];

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for k = 1:2:(numel(varargin)-1)
    if k ~= numel(varargin)
        eval([varargin{k} '=varargin{k+1};']);
    end
end

if displayFigures
    
    figure('Units', 'normalized', 'Position',[0.6441 0.0744 0.3184 0.3844]);
    tiledlayout('flow', 'TileSpacing', 'none', 'Padding', 'none')
    
end

% figure; imagesc(bw);
xDim = size(mask, 2);
yDim = size(mask, 1);

boundaryCell = bwboundaries(mask, 8, 'noholes');
stats = regionprops(~~mask, 'EquivDiameter', 'SubarrayIdx', 'Image');
averageEquivRadius = median([stats.EquivDiameter]/2);
%if an object is within borderThresh px of the image edge, let's not fit a circle to it and leave it be
borderThresh = averageEquivRadius; 
border = borderImage(mask);
borderDist = bwdist(border);

edgeMask = false(yDim, xDim);
cMask =  false(yDim, xDim);

ellipseFrame = [];
edgeEllipseFrame = [];

n = 0;
m = 0;
for k = 1:numel(boundaryCell)
    %     hold on
    xs = boundaryCell{k}(:, 1);
    ys = boundaryCell{k}(:, 2);
    %     plot(xs, ys, 'g.');
    
    if ~doEllipse
        %[x center y center R]
        [xfit,yfit, Rfit]= circfit(xs,ys);
    else
        ellipseParams = fitEllipse(boundaryCell{k});
        xfit = ellipseParams(1);
        yfit = ellipseParams(2);
        afit = ellipseParams(3);
        bfit = ellipseParams(4);
        thetafit = ellipseParams(5);
        Rfit = mean([afit, bfit]);
    end
    
    xSub = min(round(abs(xfit)), xDim);
    ySub = min(round(abs(yfit)), yDim);
    
    xSub = max(xSub, 1);
    ySub = max(ySub, 1);
    
    isFarFromBorder =  borderDist(ySub, xSub) > borderThresh;
    
    if isFarFromBorder
        
        n = n + 1;
        
        ellipseFrame(n, 2) = xfit;
        ellipseFrame(n, 1) = yfit;
        if ~doEllipse
            ellipseFrame(n, 3) = Rfit;
        else
            ellipseFrame(n, 3) = afit;
            ellipseFrame(n, 4) = bfit;
            ellipseFrame(n, 5) = thetafit + pi/2;
%             h = drawellipse('Center',[ellipseFrame(n, 1) ellipseFrame(n, 2)],'SemiAxes',[ellipseFrame(n, 3) ellipseFrame(n, 4)], ...
%             'RotationAngle',ellipseFrame(n, 5) * (360/(2*pi)),'StripeColor','m');
                 h = images.roi.Ellipse('Center',[ellipseFrame(n, 1) ellipseFrame(n, 2)],'SemiAxes',[ellipseFrame(n, 3) ellipseFrame(n, 4)], ...
              'RotationAngle',ellipseFrame(n, 5) * (360/(2*pi)),'StripeColor','m');
%             cMask = cMask + createMask(h);  
            cMask = cMask + poly2mask(h.Vertices(:, 1), h.Vertices(:, 2), size(cMask, 1), size(cMask, 2));
        end
        
    else
        
        m = m + 1;
        edgeEllipseFrame(m, 2) = xfit;
        edgeEllipseFrame(m, 1) = yfit;
        edgeEllipseFrame(m, 3) = Rfit;
        if doEllipse
            edgeEllipseFrame(m, 4) = Rfit;
            edgeEllipseFrame(m, 5) = 0;
        end
        
        reg = stats(k).Image;
        %         figure(87); imagesc(reg);
        reg = reg(:);
        % %             figure(88); imagesc(edgeMask);
        a = [stats(k).SubarrayIdx];
        sz = size(stats(k).Image);
        for xx = 1:numel(a{2})
            for yy = 1:numel(a{1})
                edgeMask(a{1}(yy), a{2}(xx)) = reg(sub2ind(sz,yy, xx));
            end
        end
        
    end
        
    %make circle with center
    %     rectangle('position',[xfit-Rfit,yfit-Rfit,Rfit*2,Rfit*2],...
    %     'curvature',[1,1],'linestyle','-','edgecolor','r');
end


if ~doEllipse
    cMask = makeNuclearMask(ellipseFrame, [size(mask, 1), size(mask, 2)], 'radiusScale', 1);
end

ellipseFrameWithEdges = cat(1, ellipseFrame, edgeEllipseFrame);

if ~isreal(ellipseFrameWithEdges)
    keyboard
end
cMask = cMask + edgeMask;
% figure; imshowpair(bw, cMask, 'montage');