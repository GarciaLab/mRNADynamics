function [innerMask, ellipseFrameWithEdges] = fitEllipsesToNuclei(mask, varargin)

warning('off','all'); %disabling all of the linear algebra spam

displayFigures = false;
doEllipse = true;
image = []; %#ok<NASGU>
maxAspectRatio = 4; %quality control to remove bad ellipses. major/minor
areaFilter = [];

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for k = 1:2:(numel(varargin)-1)
    if k ~= numel(varargin)
        eval([varargin{k} '=varargin{k+1};']);
    end
end

if displayFigures
    
    figure('Units', 'normalized', 'Position',[0.6441 0.0744 0.3184 0.3844]); %#ok<*UNRCH>
    tiledlayout('flow', 'TileSpacing', 'none', 'Padding', 'none')
    
end

% figure; imagesc(bw);
xDim = size(mask, 2);
yDim = size(mask, 1);

boundaryCell = bwboundaries(mask, 8, 'noholes');
stats = regionprops(~~mask, 'EquivDiameter', 'SubarrayIdx', 'Image');
averageEquivRadius = median([stats.EquivDiameter]/2);

%if an object is within borderThresh px of the image edge,
%let's not fit a circle to it and leave it be
borderThresh = averageEquivRadius;
border = borderImage(mask);
borderDist = bwdist(border);

edgeMask = false(yDim, xDim);
innerMask =  false(yDim, xDim);

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
        try
            ellipseParams = fitEllipse(boundaryCell{k});
        catch
            continue;
        end
        %not sure why complex parameters or NaNs happen. we'll skip this ellipse
        %in this case. also ensure good aspect ratios here
        if ~isreal(ellipseParams) ||...
                ellipseParams(3)/ellipseParams(4) > maxAspectRatio ||...
                ellipseParams(4)/ellipseParams(3) > maxAspectRatio ||...
                any(isnan(ellipseParams))
            
            
            continue;
            
        end
        
        xfit = ellipseParams(1);
        yfit = ellipseParams(2);
        afit = ellipseParams(3);
        bfit = ellipseParams(4);
        thetafit = ellipseParams(5);
        Rfit = mean([afit, bfit]);
        
    end
    
    xSub = max(min(round(abs(xfit)), xDim), 1);
    ySub = max(min(round(abs(yfit)), yDim), 1);
    
    isFarFromBorder =  borderDist(ySub, xSub) > borderThresh;
    
    if isFarFromBorder
        
        n = n + 1;
        
        ellipseFrame(n, 2) = xfit;
        ellipseFrame(n, 1) = yfit;
        if ~doEllipse
            ellipseFrame(n, 3) = Rfit;
            ellipseFrame(n, 5) = 0; %angle
        else
            ellipseFrame(n, 3) = afit;
            ellipseFrame(n, 4) = bfit;
            ellipseFrame(n, 5) = mod(thetafit + pi/2, 2*pi); %rotate theta to match
            %what roi.Ellipse expects
            
            h = images.roi.Ellipse('Center',[ellipseFrame(n, 1) ellipseFrame(n, 2)],...
                'SemiAxes',[ellipseFrame(n, 3) ellipseFrame(n, 4)], ...
                'RotationAngle',ellipseFrame(n, 5) * (360/(2*pi)),'StripeColor','m');
            innerMask = innerMask + poly2mask(h.Vertices(:, 1), h.Vertices(:, 2), size(innerMask, 1), size(innerMask, 2));
            %             h = drawellipse('Center',[ellipseFrame(n, 1) ellipseFrame(n, 2)],'SemiAxes',[ellipseFrame(n, 3) ellipseFrame(n, 4)], ...
            %             'RotationAngle',ellipseFrame(n, 5) * (360/(2*pi)),'StripeColor','m');
            %             cMask = cMask + createMask(h);
            
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
        
        region = stats(k).Image;
        %         figure(87); imagesc(reg);
        region = region(:);
        % %             figure(88); imagesc(edgeMask);
        a = [stats(k).SubarrayIdx];
        sz = size(stats(k).Image);
        for xx = 1:numel(a{2})
            for yy = 1:numel(a{1})
                edgeMask(a{1}(yy), a{2}(xx)) = region(sub2ind(sz,yy, xx));
            end
        end
        
    end
    
end

%this is here for backwards compatibility.
%it's being deprecated.
if ~doEllipse
    innerMask = makeNuclearMask(ellipseFrame,...
        [size(mask, 1), size(mask, 2)], 'radiusScale', 1);
end

% ellipseFrameWithEdges = cat(1, ellipseFrame, edgeEllipseFrame);
ellipseFrameWithEdges = ellipseFrame;

innerMask = innerMask + edgeMask;
% 
% ellipseFrameWithEdgesTemp = [];
% %quality control
% if ~isempty(areaFilter)
%     for n = 1:size(ellipseFrameWithEdges, 1)
%         ellipseArea = pi*ellipseFrameWithEdges(n, 3)*ellipseFrameWithEdges(n, 4);
%         ellipseAspectRatio = ellipseFrameWithEdges(n, 3) / ellipseFrameWithEdges(n, 4); 
%         
%         if ellipseArea > areaFilter(1) &&...
%                ellipseArea < areaFilter(2) &&...
%                 ellipseAspectRatio < maxAspectRatio &&...
%                 (1/ellipseAspectRatio) < maxAspectRatio
%                 
%             ellipseFrameWithEdgesTemp = [ellipseFrameWithEdgesTemp; ellipseFrameWithEdges(n, :)];
%             
%         end
%         
%         assert(size(ellipseFrameWithEdgesTemp, 1) <= size(ellipseFrameWithEdges, 1));
%         
%     end
% end

% ellipseFrameWithEdges = ellipseFrameWithEdgesTemp;

%double-check orientation angles are kosher
if ~isempty(ellipseFrameWithEdges)
    assert( all(ellipseFrameWithEdges(:, 5) <= 2*pi) );
end


% figure; imshowpair(bw, cMask, 'montage');


end