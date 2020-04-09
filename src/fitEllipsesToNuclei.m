function [cMask, ellipseFrameWithEdges] = fitEllipsesToNuclei(mask, varargin)

displayFigures = false;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

if displayFigures
    
    figure('Units', 'normalized', 'Position',[0.6441 0.0744 0.3184 0.3844]);
    tiledlayout('flow', 'TileSpacing', 'none', 'Padding', 'none')
    
end

% figure; imagesc(bw);
xDim = size(mask, 1);
yDim = size(mask, 1);

boundaryCell = bwboundaries(mask, 8, 'noholes');
stats = regionprops(~~mask, 'EquivDiameter', 'SubarrayIdx', 'Image');
averageEquivRadius = median([stats.EquivDiameter]/2);
%if an object is within borderThresh px of the image edge, let's not fit a circle to it and leave it be
borderThresh = averageEquivRadius; 
border = borderImage(mask);
borderDist = bwdist(border);

edgeMask = false(xDim, yDim);

ellipseFrame = zeros(numel(boundaryCell), 3);
edgeEllipseFrame = [];

for i = 1:numel(boundaryCell)
    %     hold on
    xs = boundaryCell{i}(:, 1);
    ys = boundaryCell{i}(:, 2);
    %     plot(xs, ys, 'g.');
    
    %[x center y center R]
%     [xfit,yfit, Rfit]= circfit(xs,ys);
    ellipse_t = fit_ellipse(xs, ys);
    
    yFit = ellipse_t.Y0_in;
    xFit = ellipse_t.X0_in;
    
    xSub = min(round(abs(xFit)), xDim);
    ySub = min(round(abs(yFit)), yDim);
    phi = ellipse_t.phi;
    longAxis = ellipse_t.long_axis;
    shortAxis = ellipse_t.short_axis;
    
    
    isFarFromBorder =  borderDist(ySub, xSub) > borderThresh;
    
    if isFarFromBorder
        
        ellipseFrame(i, 1) = yfit;
        ellipseFrame(i, 2) = xfit;
        ellipseFrame(i, 3) = longAxis;
        ellipseFrame(i, 4) = shortAxis;
        ellipseFrame(i, 5) = phi;
        
    else
        edgeEllipseFrame(i, 1) = yfit;
        edgeEllipseFrame(i, 2) = xfit;
        edgeEllipseFrame(i, 3) = longAxis;
        edgeEllipseFrame(i, 4) = shortAxis;
        edgeEllipseFrame(i, 5) = phi;        
        
        reg = stats(i).Image;
        %         figure(87); imagesc(reg);
        reg = reg(:);
        % %             figure(88); imagesc(edgeMask);
        a = [stats(i).SubarrayIdx];
        sz = size(stats(i).Image);
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

cMask = makeNuclearMask(ellipseFrame, [size(mask, 1), size(mask, 2)], 'radiusScale', 1);

ellipseFrameWithEdges = cat(1, ellipseFrame, edgeEllipseFrame);

cMask = cMask + edgeMask;
% figure; imshowpair(bw, cMask, 'montage');