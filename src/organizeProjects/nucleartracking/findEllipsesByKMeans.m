function [centers, radii, fin] = findEllipsesByKMeans(im, varargin)

pixelSize = .1;
min_area = pi*(1/pixelSize)^2; 
max_area = pi*(6/pixelSize)^2;
displayFigures = false;

for i = 1:(numel(varargin)-1)
    if i ~= numel(varargin)
        if ~ischar(varargin{i+1})
            eval([varargin{i} '=varargin{i+1};']);
        end
    end
end

df = displayFigures;
yDim = size(im, 1);
xDim = size(im, 2);

se = strel('disk', 4);

top = @(X) imtophat(X, strel('disk',30));
dil = @(X) imdilate(X, se);
sk2 = @(X) imsegkmeans(dil(top(X)), 2);

[smths2, outs2] = smoothTile(im, [2, 16], 10, 'fun', sk2, 'displayFigures', df);
sk3 = @(X) imsegkmeans(dil(X), 3);

[smths3, outs3] = smoothTile(im, [2, 16], 10, 'fun', sk3,'displayFigures', df);


out2 = rescale(imsegkmeans(uint8(outs2), 2)); 
fin2 = wshed(out2,  'displayFigures', df); 

if df
    figure(4); imagesc(out2); figure(6); imagesc(fin2);
end

for n = 1:size(outs3, 3)
    
    a = squeeze(outs3(:, :, n));
%     b = rescale(label2rgb(a));
    ch1 = a==1;
    ch2 = a==2;
    ch3 = a==3;
%     figure; imagesc(ch1); figure; imagesc(ch2); figure; imagesc(ch3)
    stats1 = regionprops(ch1, 'all'); stats2 = regionprops(ch2, 'all'); stats3 = regionprops(ch3, 'all');
%     figure; histogram([stats1.Area]); figure; histogram([stats2.Area]); figure; histogram([stats3.Area]);
    a1 = max([stats1.Area]); a2 = max([stats2.Area]); a3 = max([stats3.Area]);
    e1 = min(([stats1.EulerNumber])); e2 = min([stats2.EulerNumber]); e3 = min([stats3.EulerNumber]);
    as = [a1, a2, a3];
    es = [e1, e2, e3];
    vals = [1 2 3];
    
    [m, bg] = min(es);
    
    outs3(:, :, n) = outs3(:, :, n) ~= bg;
    
end

out3 = rescale(imsegkmeans(uint8(outs3), 2));
fin3 = wshed(out3, 'displayFigures', df);

if df
    figure(7); imagesc(fin3);  figure(5); imagesc(out3)
end

stats2 = regionprops(imclearborder(fin2), 'all'); stats3 = regionprops(imclearborder(fin3), 'all');

u2 = std([stats2.EquivDiameter]);
u3 = std([stats3.EquivDiameter]);

fin = fin2; 
if u2>u3
    fin = fin3;
end
% if length(stats2)<length(stats3)
%     fin = fin3;
% end

cs= regionprops(logical(fin), 'Centroid');
radii = regionprops(logical(fin), 'EquivDiameter');
radii = [radii.EquivDiameter]/2;
centers=[];
for k = 1:length(cs)
    centers(k, :) = cs(k).Centroid;
end

fin = bwareafilt(fin, [min_area, max_area]);

%% processing done
cs= regionprops(logical(fin), 'Centroid');
radii = regionprops(logical(fin), 'EquivDiameter');
radii = [radii.EquivDiameter]/2;
centers=[];
for k = 1:length(cs)
    centers(k, :) = cs(k).Centroid;
end


if df
    viscircles(centers, radii)
    figure(2); imshow(im, []); viscircles(centers, radii);
end


end