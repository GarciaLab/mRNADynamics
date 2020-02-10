function nuclearMask = makeNuclearMask(ellipsesFrame, dim, varargin)

%make a mask from the ellipses structure. this is currently used for
%segmenting loci in segmentSpots.

radScale = 1.3;

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'radScale')
        radScale = varargin{i+1};
    end
end

nuclearMask = false(dim(1), dim (2));

for e = 1:size(ellipsesFrame, 1)
    
    ceny = ellipsesFrame(e, 1);
    cenx = ellipsesFrame(e, 2);
    rad = ellipsesFrame(e,3)*radScale;
    xrange = max(ceil(cenx-rad), 1) : min(ceil(cenx+rad),dim(1));
    yrange = max(ceil(ceny-rad), 1) : min(ceil(ceny+rad),dim(2));
%     a = [xrange' yrange'];
%     v = [cenx, ceny];
%     d = a - v;
%     b = meshgrid(xrange, yrange);
%     nuclearMask(b) = 0;
    for x = xrange
        for y = yrange
            nuclearMask(x, y) =  sqrt(([x, y] - [cenx, ceny])*([x, y] - [cenx, ceny])') < rad;
%             nuclearMask(x, y) = norm([x, y] - [cenx, ceny]) < rad;
        end
    end
    
end




end


%%
% Different algorithm i might implement later 

% %Create the circle that we'll use as the mask
% IntegrationRadius=2;       %Radius of the integration region in um
% IntegrationRadius=floor(IntegrationRadius/FrameInfo(1).PixelSize); %Radius of the integration in pixels
% if ~mod(IntegrationRadius,2)
%     IntegrationRadius=IntegrationRadius+1;
% end
% 
% Circle=false(3*IntegrationRadius,3*IntegrationRadius);
% Circle=MidpointCircle(Circle,IntegrationRadius,1.5*IntegrationRadius+0.5,...
%     1.5*IntegrationRadius+0.5,1);
% 
% refFrame = ones(LinesPerFrame,PixelsPerLine, nSlices);
%     convRef = convn(refFrame, Circle, 'same');
%     edgeMask = convRef~=sum(Circle(:)); 
% convImage = imfilter(Image, double(Circle), 'same');
%             convImage(edgeMask) = NaN;
%              for j=1:length(tempSchnitz)
%                 CurrentIndex=find(tempSchnitz(j).frames==CurrentFrame);
%                 cenx=min(max(1,round(tempSchnitz(j).cenx(CurrentIndex))),PixelsPerLine);
%                 ceny=min(max(1,round(tempSchnitz(j).ceny(CurrentIndex))),LinesPerFrame);
%                 tempSchnitz(j).Fluo(CurrentIndex,1:nSlices,ChN) = single(convImage(ceny,cenx,:));
%             end