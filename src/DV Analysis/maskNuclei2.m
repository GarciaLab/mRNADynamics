function [centers, radii, nuclearMask] = maskNuclei2(raw, varargin)

displayFigures = true;

for i = 1:(numel(varargin)-1)
    if i ~= numel(varargin)
        if ~ischar(varargin{i+1})
            eval([varargin{i} '=varargin{i+1};']);
        end
    end
end

d = struct;
d.raw = raw;

% radiiRange = [40, 80];
% scaleRange = radiiRange(1):round(radiiRange(2));
% scaleRange = 10:15;
% n = 0;
% for scale = scaleRange
%     n = n + 1;

%the meat of the operation- the canny filter
scale = 12;
% d.imCanny = canny(rescale(d.raw), scale);
% 
% %the rest is just cleaning up the edges and filtering out non-toroidal objects
% d.fib = fibermetric(double(d.imCanny));
% d.fibSkel = bwmorph(d.fib, 'skel');
% imTemp = d.fibSkel;
% for k = 1:2
%     imTemp = bwmorph(imTemp, 'shrink');
%     imTemp = bwmorph(imTemp, 'clean');
%     imTemp = bwmorph(imTemp, 'bridge');
%     imTemp = bwmorph(imTemp, 'diag');
% end
% 
% 
% d.cleanedFibSkel= imTemp;
% se = strel('disk', 2);
% nCircles = [];
% nCircles = [nCircles, length(regionprops(logical(d.cleanedFibSkel)))];
% imTemp = d.cleanedFibSkel;
% chg = 0;
% tol = 1; %process stays constant for 5 steps
% while chg < tol
%     d.dil = imdilate(imTemp, se);
%     figure(3); imshow(d.dil);
%     d.fill= bwmorph(d.dil, 'fill');
%     figure(4); imshow(d.fill);
%     d.donut = bwpropfilt(d.fill,'EulerNumber',[0, 0]);
%     figure(5); imshow(d.donut);
%     nCircles = [nCircles, length(regionprops(logical(d.donut)))];
%     
%     if nCircles(end) == nCircles(end-1)
%         chg = chg + 1;
%     else
%         chg = 0;
%     end
%     if chg == tol
%         break;
%     end
%     
% end
% 
% 
% %%
% if displayFigures
%     ims = fieldnames(d);
%     imageCell = {};
%     for i = 1:length(ims)
%         imageCell{i} = d.(ims{i});
%     end
%     imageTile(imageCell, ims);
% end
% 
% 
% im2 = im2.*a;
% im2 = fibermetric(double(im2));
% im2 = bwmorph(im2, 'bridge');
% im2 = bwmorph(im2, 'diag');
% 
% 
% for k = 1:100
%     im2 = bwmorph(im2, 'skel');
%     im2 = bwmorph(im2, 'spur');
% end

% figure(1); imshow(im2, [])

%     a = imfill(im2, 'holes');
%     figure(4); imshow(a, [])
d.imk = imsegkmeans(imgaussfilt(d.raw,scale/2),2); %
figure(1); imagesc(d.imk)
a=rescale(d.imk);
cs= regionprops(logical(a), 'Centroid');
radii = regionprops(logical(a), 'EquivDiameter');
radii = [radii.EquivDiameter]/2;
centers=[];
for k = 1:length(cs)
    centers(k, :) = cs(k).Centroid;
end

nuclearMask = filledVisCircles(d.imk, centers, radii);

% nuclearMask = filledVisCircles(im2, centers, radii);
% 
viscircles(centers, radii)
figure(2); imshow(d.raw, []); viscircles(centers, radii);

end

% end