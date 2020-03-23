function [cMask, ellipseFrameWithEdges] = fitCirclesToNuclei(mask, varargin)

displayFigures = false;
image = [];

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

ellipseFrame = [];
edgeEllipseFrame = [];

n = 0;
m = 0;
for i = 1:numel(boundaryCell)
    %     hold on
    xs = boundaryCell{i}(:, 1);
    ys = boundaryCell{i}(:, 2);
    %     plot(xs, ys, 'g.');
    
    %[x center y center R]
    [xfit,yfit, Rfit]= circfit(xs,ys);
    xSub = min(round(abs(xfit)), xDim);
    ySub = min(round(abs(yfit)), yDim);
    
    xSub = max(xSub, 1);
    ySub = max(ySub, 1);
    
    isFarFromBorder =  borderDist(ySub, xSub) > borderThresh;
    
    if isFarFromBorder
        
        n = n + 1;
        
        ellipseFrame(n, 2) = xfit;
        ellipseFrame(n, 1) = yfit;
        ellipseFrame(n, 3) = Rfit;
        
    else
        
        m = m + 1;
        edgeEllipseFrame(m, 2) = xfit;
        edgeEllipseFrame(m, 1) = yfit;
        edgeEllipseFrame(m, 3) = Rfit;
        
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